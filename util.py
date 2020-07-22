"""
helper functions for Helmsman
"""

# system packages
from __future__ import print_function
import os
import sys
import warnings
import itertools
import collections
import csv
from joblib import Parallel, delayed
from logging import StreamHandler, getLogger as realGetLogger, Formatter
from colorama import Fore, Back, Style
# matrix+stats processing
import pandas as pd
import numpy as np
# vcf/fasta parsing
from cyvcf2 import VCF
from pyfaidx import Fasta
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
# PCA algorithms
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
# ignore nuisance warnings when loading nimfa package
warnings.filterwarnings("ignore", category=UserWarning)
# decomposition algorithms
import nimfa

sys.path.append(os.getcwd())


###############################################################################
# Configure color stream handler
# https://gist.github.com/jonaprieto/a61d9cade3ba19487f98
###############################################################################
class ColourStreamHandler(StreamHandler):
    """ A colorized output StreamHandler """

    # Some basic colour scheme defaults
    colours = {
        'DEBUG': Fore.CYAN,
        'INFO': Fore.GREEN,
        'WARN': Fore.YELLOW,
        'WARNING': Fore.YELLOW,
        'ERROR': Fore.RED,
        'CRIT': Back.RED + Fore.WHITE,
        'CRITICAL': Back.RED + Fore.WHITE
    }

    def emit(self, record):
        try:
            message = self.format(record)
            self.stream.write(self.colours[record.levelname] + message +
                              Style.RESET_ALL)
            self.stream.write(getattr(self, 'terminator', '\n'))
            self.flush()
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.handleError(record)


###############################################################################
# configure logger
###############################################################################
# class initLogger:
#     """ initialize logger """
#     def __init__(level):
#         self.level = level


def get_logger(name=None,
               fmt='[%(name)s::%(funcName)s] %(levelname)s %(message)s',
               level='INFO'):
    """ Get and initialize a colourised logging instance if the system supports
    it as defined by the log.has_colour
    :param name: Name of the logger
    :type name: str
    :param fmt: Message format to use
    :type fmt: str
    :return: Logger instance
    :rtype: Logger
    """
    log = realGetLogger(name)
    # Only enable colour if support was loaded properly
    handler = ColourStreamHandler()
    handler.setLevel(level)
    handler.setFormatter(Formatter(fmt))
    log.addHandler(handler)
    log.setLevel(level)
    log.propagate = 0  # Don't bubble up to the root logger
    return log


util_log = get_logger(__name__, level="DEBUG")


###############################################################################
# Manipulate sequence motifs etc.
###############################################################################
def getCategory(mu_type):
    """
    collapse mutation types per strand symmetry

    """

    # if re.match("^[ACGT]*$", mu_type):
    if mu_type in ('AC', 'TG'):
        category = "T_G"
    elif mu_type in ('AG', 'TC'):
        category = "T_C"
    elif mu_type in ('AT', 'TA'):
        category = "T_A"
    elif mu_type in ('CA', 'GT'):
        category = "C_A"
    elif mu_type in ('CG', 'GC'):
        category = "C_G"
    elif mu_type in ('CT', 'GA'):
        category = "C_T"
    else:
        category = "unknown"
    return category


def getMotif(sequence):
    """
    query reference genome for local sequence motif

    """

    motif = Seq(sequence, IUPAC.unambiguous_dna)
    altmotif = motif.reverse_complement()
    central_base_pos = (len(motif) - 1) // 2

    central_base = motif[central_base_pos]

    if central_base in ('C', 'T'):
        motif_a = motif
    else:
        motif_a = altmotif

    return motif_a


def indexSubtypes(motiflength):
    """
    define k-mer mutation subtypes

    """

    categories = ["T_G", "T_C", "T_A", "C_G", "C_T", "C_A"]
    bases = ["A", "C", "G", "T"]
    flank = (motiflength - 1) // 2

    if motiflength > 1:
        kmers = itertools.product(bases, repeat=motiflength - 1)

        subtypes_list = []

        for kmer in kmers:
            kmerstr = ''.join(kmer)

            for category in categories:
                ref = category[0]

                subtype = category + "." \
                    + kmerstr[0:flank] + ref + kmerstr[flank:(motiflength-1)]

                subtypes_list.append(subtype)
    else:
        ext = [".T", ".C"]
        extr = list(np.repeat(ext, 3))
        subtypes_list = [m + n for m, n in zip(categories, extr)]

    i = 0
    subtypes_dict = {}
    for subtype in sorted(subtypes_list):
        subtypes_dict[subtype] = i
        i += 1

    util_log.debug("%s %s-mer subtypes indexed", len(subtypes_dict.keys()),
                   motiflength)

    return subtypes_dict


def indexGroups(samplefile, groupvar):
    """
    Build dictionary with sample ID as key, group ID as value

    """

    sg_dict = {}

    f = open(samplefile, 'r', encoding="utf-8")
    reader = csv.DictReader(f, delimiter='\t')

    for row in reader:
        sg_dict[row['ID']] = row[groupvar]

    return sg_dict


def get_samples(sample_file):
    """
    get samples from input M matrix when using aggregation mode

    """

    samples = np.loadtxt(
        sample_file, dtype='S120', skiprows=1, delimiter='\t', usecols=(0, ))

    util_log.debug("%s contains %s samples", sample_file, len(samples))

    return samples


def parseSampleFile(samplefile):
    """
    get list of samples to keep if samplefile supplied

    """

    # f = open(args.input, 'r', encoding = "ISO-8859-1")
    f = open(samplefile, 'r', encoding="utf-8")
    reader = csv.DictReader(f, delimiter='\t')
    keep_samples = []
    for row in reader:
        keep_samples.append(row['ID'])

    return keep_samples


def get_samples_vcf(args, inputvcf):
    """
    get samples from VCF file

    """

    if args.samplefile:
        keep_samples = parseSampleFile(args.samplefile)
        vcf_reader = VCF(
            inputvcf, mode='rb', gts012=True, lazy=True, samples=keep_samples)
    else:
        vcf_reader = VCF(inputvcf, mode='rb', gts012=True, lazy=True)

    if (args.samplefile and args.groupvar):
        samples = indexGroups(args.samplefile, args.groupvar)
    else:
        samples = vcf_reader.samples

    return samples


class processInput:
    """
    Methods for parsing input data into sample x subtype count matrices:
    - MAF format
    - plain text format
    - Aggregation of existing subtype count matrices

    """

    def __init__(self, mode, args, subtypes_dict, par=False):
        self.mode = mode
        self.args = args
        self.subtypes_dict = subtypes_dict
        self.par = par

        if self.mode == "agg":
            self.data = self.process_agg()
        elif self.mode == "txt":
            self.data = self.process_txt()
        elif self.mode == "maf":
            self.data = self.process_maf()
        elif self.mode == "vcf":
            if (args.input.lower().endswith(('.vcf', '.vcf.gz', '.bcf'))
                    or args.input == "-"):
                par = False
                self.data = self.process_vcf(args.input)
            elif args.input.lower().endswith(('.txt')):
                self.par = True
                with open(args.input) as vcf_list_file:
                    vcf_list = vcf_list_file.read().splitlines()

                results = Parallel(n_jobs=args.cpus) \
                    (delayed(self.process_vcf)(vcf) \
                    for vcf in vcf_list)

                if args.rowwise:
                    count_matrix = np.vstack(results)

                    samples = np.array([])
                    for vcf in vcf_list:
                        samples = np.append(samples, get_samples_vcf(args, vcf))

                else:
                    nrow, ncol = results[1].shape
                    count_matrix = np.zeros((nrow, ncol))

                    for count_matrix_i in results:
                        count_matrix = np.add(count_matrix, count_matrix_i)
                    self.par
                    samples = np.array([get_samples_vcf(args, vcf_list[1])])
                self.data = collections.namedtuple('Out', ['M', 'samples'])(
                    count_matrix, samples)

    def process_vcf(self, inputfile):
        """
        Main function for parsing VCF

        """
        # initialize reference genome
        fasta_reader = Fasta(self.args.fastafile, read_ahead=1000000)

        # initialize vcf reader
        if self.args.samplefile:
            keep_samples = parseSampleFile(self.args.samplefile)

            vcf_reader = VCF(
                inputfile,
                mode='rb',
                gts012=True,
                lazy=True,
                samples=keep_samples)
        else:
            vcf_reader = VCF(inputfile, mode='rb', gts012=True, lazy=True)

        nbp = (self.args.length - 1) // 2

        # index samples
        if (self.args.samplefile and self.args.groupvar):
            all_samples = vcf_reader.samples

            sg_dict = indexGroups(self.args.samplefile, self.args.groupvar)
            samples = sorted(list(set(sg_dict.values())))

            # get boolean vector of samples that are in sample file
            samples_keep_match = np.isin(all_samples, list(sg_dict.keys()))

            # get indices of matching samples
            samples_keep_idx = np.where(samples_keep_match)

            # get list of individual sample ids to keep
            samples_keep = sorted(list(set(sg_dict.keys())))

            util_log.debug("%s samples will be pooled into %s groups: %s",
                           len(all_samples), len(samples), ",".join(samples))
        else:
            samples = vcf_reader.samples

        samples_dict = {}
        for i, sample in enumerate(samples):
            samples_dict[sample] = i

        # Query records in VCF and build matrix
        M = np.zeros((len(samples), len(self.subtypes_dict)))
        numsites_keep = 0
        numsites_skip = 0
        chrseq = '0'
        chr_check = "none"

        for record in vcf_reader:

            # Filter by SNP status, # alt alleles, and FILTER column
            if (not record.is_snp or len(record.ALT) != 1
                    or record.FILTER is not None):
                numsites_skip += 1
                continue

            # Filter by allele count
            if record.INFO['AC'] > self.args.maxac > 0:
                numsites_skip += 1
                continue

            row_chr = record.CHROM

            # check chromosome formatting matches between MAF and fasta files
            if numsites_keep == 0:
                if "chr1" in fasta_reader and "chr" not in row_chr:
                    chr_check = "add"
                    util_log.debug(
                        "formatting mismatch: 'chr' only in fasta file")
                elif "chr1" not in fasta_reader and "chr" in row_chr:
                    chr_check = "delete"
                    util_log.debug(
                        "formatting mismatch: 'chr' only in MAF file")
                else:
                    util_log.debug("chromosome formatting matches")

            if chr_check == "add":
                row_chr = "chr" + row_chr
            elif chr_check == "delete":
                row_chr = row_chr.replace('chr', '')

            if row_chr != chrseq:
                sequence = fasta_reader[row_chr]
                chrseq = row_chr

            # check and update chromosome sequence
            # if record.CHROM != chrseq:
            #     sequence = fasta_reader[record.CHROM]
            #     chrseq = record.CHROM

            lseq = sequence[record.POS - (nbp + 1):record.POS + nbp].seq

            mu_type = record.REF + str(record.ALT[0])
            category = getCategory(mu_type)
            motif_a = getMotif(lseq)
            subtype = str(category + "." + motif_a)

            if subtype not in self.subtypes_dict:
                numsites_skip += 1
                continue

            st = self.subtypes_dict[subtype]

            # currently only works with singletons--
            if (self.args.samplefile and self.args.groupvar):

                gt_new = record.gt_types

                if (self.args.impute and 3 in gt_new):
                    gt_complete = gt_new[gt_new != 3]
                    freq = sum(gt_complete) / len(gt_complete)
                    gt_new[gt_new == 3] = freq

                else:
                    gt_new[gt_new == 3] = 0

                # if not any("/" in b for b in record.gt_bases):
                if self.args.haploid:
                    gt_new = np.divide(gt_new, 2.)

                # get array of genotypes only for samples in samplefile
                gt_sub = gt_new[samples_keep_idx]

                if gt_sub.sum() == 0:
                    numsites_skip += 1
                    continue

                # initialize dict of group allele counts = 0
                sg_counts = {k: 0 for k in sorted(list(set(sg_dict.values())))}

                # initialize dict of allele counts per sample
                d2 = dict(zip(samples_keep, gt_sub))

                # iterate per-sample counts and update per-group counts
                for key, value in d2.items():
                    sg_counts[sg_dict[key]] += value

                # add to matrix
                M[:, st] = M[:, st] + list(sg_counts.values())
                numsites_keep += 1

            else:
                gt_new = record.gt_types
                if (self.args.impute and 3 in gt_new):
                    gt_complete = gt_new[gt_new != 3]
                    freq = sum(gt_complete) / len(gt_complete)
                    gt_new[gt_new == 3] = freq

                else:
                    gt_new[gt_new == 3] = 0

                # if not any("/" in b for b in record.gt_bases):
                if self.args.haploid:
                    gt_new = np.divide(gt_new, 2.)

                M[:, st] = M[:, st] + gt_new
                numsites_keep += 1
                # util_log.debug(gt_new)

            if numsites_keep % 100000 != 0:
                continue
            util_log.debug("%s : %s sites counted", inputfile, numsites_keep)

        util_log.debug("%s : %s sites counted", inputfile, numsites_keep)
        util_log.debug("%s : %s sites skipped", inputfile, numsites_skip)

        out = collections.namedtuple('Out', ['M', 'samples'])(M, samples)
        if self.par:
            out = M

        return out

    def process_maf(self):
        """
        process MAF files

        """

        fasta_reader = Fasta(self.args.fastafile, read_ahead=1000000)

        nbp = (self.args.length - 1) // 2
        samples_dict = {}

        # M = np.zeros((len(samples), len(subtypes_dict)))
        numsites_keep = 0
        numsites_skip = 0
        chrseq = '0'

        maf_file = open(self.args.input, 'r', encoding="ISO-8859-1")

        reader = csv.DictReader(
            filter(lambda row: row[0] != '#', maf_file), delimiter='\t')
        counter = 0
        chr_check = "none"
        for row in reader:

            if (row['Variant_Type'] not in ["SNP", "SNV"]):
                continue

            if 'Start_Position' in row:
                pos = int(row['Start_Position'])
            else:
                pos = int(row['Start_position'])
            ref = row['Reference_Allele']
            alt = row['Tumor_Seq_Allele2']
            row_chr = row['Chromosome']
            sample = row[self.args.groupvar]

            # check chromosome formatting matches between MAF and fasta files
            if counter == 0:
                if "chr1" in fasta_reader and "chr" not in row_chr:
                    chr_check = "add"
                    util_log.debug(
                        "formatting mismatch: 'chr' only in fasta file")
                elif "chr1" not in fasta_reader and "chr" in row_chr:
                    chr_check = "delete"
                    util_log.debug(
                        "formatting mismatch: 'chr' only in MAF file")
                else:
                    util_log.debug("chromosome formatting matches")

            if chr_check == "add":
                row_chr = "chr" + row_chr
            elif chr_check == "delete":
                row_chr = row_chr.replace('chr', '')

            if row_chr != chrseq:
                sequence = fasta_reader[row_chr]
                chrseq = row_chr

            # if row['Chromosome'] != chrseq:
            #     sequence = fasta_reader[row['Chromosome']]
            #     chrseq = row['Chromosome']

            counter += 1
            mu_type = ref + alt
            category = getCategory(mu_type)
            lseq = sequence[pos - (nbp + 1):pos + nbp].seq

            motif_a = getMotif(lseq)
            subtype = str(category + "." + motif_a)
            # st = subtypes_dict[subtype]

            if sample not in samples_dict:
                samples_dict[sample] = {}

            if subtype not in samples_dict[sample]:
                samples_dict[sample][subtype] = 1
            else:
                samples_dict[sample][subtype] += 1

            if counter % 1000 != 0:
                continue
            util_log.debug("%s : %s sites counted", self.args.input, counter)

        M = pd.DataFrame(samples_dict).T.fillna(0).values
        samples = sorted(samples_dict)

        out = collections.namedtuple('Out', ['M', 'samples'])(M, samples)
        return out

    def process_agg(self):
        """
        aggregate M matrices from list of input files

        """

        inputM = self.args.input
        colnames = ["ID"]
        M_colnames = colnames + list(sorted(self.subtypes_dict.keys()))
        colrange = range(1, len(M_colnames))

        if (inputM.lower().endswith('m_samples.txt')
                or inputM.lower().endswith('m_regions.txt')):
            with open(inputM) as f:
                file_list = f.read().splitlines()

            # M output by sample
            if inputM.lower().endswith('m_samples.txt'):

                M_out = np.array([M_colnames])

                samples = np.empty((0, 100))

                for mfile in file_list:
                    samples_it = get_samples(mfile)
                    samples = np.concatenate((samples, samples_it), axis=None)

                    M_it = np.loadtxt(mfile, skiprows=1, usecols=colrange)
                    M_it = np.concatenate((np.array([samples_it]).T, M_it),
                                          axis=1)
                    M_out = np.concatenate((M_out, M_it), axis=0)

                M = np.delete(M_out, 0, 0)
                M = np.delete(M, 0, 1)
                M = M.astype(np.float)

            # M output by region
            elif inputM.lower().endswith('m_regions.txt'):
                samples = get_samples(file_list[0])

                M_out = np.zeros((len(samples), len(M_colnames) - 1))
                for mfile in file_list:
                    M_it = np.loadtxt(mfile, skiprows=1, usecols=colrange)
                    M_out = np.add(M_out, M_it)

                M = M_out.astype(np.float)

        else:
            samples = get_samples(inputM)
            M = np.loadtxt(inputM, skiprows=1, usecols=colrange)
            M = M.astype(np.float)

        out = collections.namedtuple('Out', ['M', 'samples'])(M, samples)
        return out

    def process_txt(self):
        """
        process tab-delimited text file, containing the following columns:
        CHR    POS    REF    ALT    SAMPLE_ID

        """

        fasta_reader = Fasta(self.args.fastafile, read_ahead=1000000)

        nbp = (self.args.length - 1) // 2
        samples_dict = {}

        numsites_keep = 0
        numsites_skip = 0
        chrseq = '0'

        with open(self.args.input, 'r') as txt_file:
            reader = csv.reader(txt_file, delimiter='\t')

            for row in reader:
                chrom = row[0]
                pos = int(row[1])
                ref = row[2]
                alt = row[3]
                sample = row[4]

                if chrom != chrseq:
                    sequence = fasta_reader[chrom]
                    chrseq = chrom

                if (len(alt) == 1 and len(ref) == 1):
                    mu_type = ref + alt
                    category = getCategory(mu_type)
                    if nbp > 0:
                        lseq = sequence[pos - (nbp + 1):pos + nbp].seq
                    else:
                        lseq = sequence[pos - 1].seq
                        # eprint("lseq:", lseq)
                    motif_a = getMotif(lseq)
                    subtype = str(category + "." + motif_a)
                    
                    if subtype not in self.subtypes_dict:
                        continue

                    if sample not in samples_dict:
                        samples_dict[sample] = {}

                    if subtype not in samples_dict[sample]:
                        samples_dict[sample][subtype] = 1
                    else:
                        samples_dict[sample][subtype] += 1
            mdf = pd.DataFrame(samples_dict).T.fillna(0)
            samples = mdf.index.tolist() #instead of using samples_dict with sorted(), which leads to mismatching, simply retain the explicit ordering of the matrix dataframe.
            M = mdf.values 

        out = collections.namedtuple('Out', ['M', 'samples'])(M, samples)
        return out


class DecompModel:
    """
    Class for fitting PCA and NMF models

    """

    def __init__(self, M_run, rank, seed, decomp):
        self.M_run = M_run / (M_run.sum(axis=1) + 1e-8)[:, None]
        self.rank = rank
        self.seed = seed
        self.decomp = decomp

        self.evar_dict = {}

        if self.decomp == "pca":
            # standarize input matrix
            X_std = StandardScaler().fit_transform(self.M_run)

            # run PCA
            pca = PCA(n_components=self.M_run.shape[1])
            W = pca.fit_transform(X_std)
            H = pca.components_.T * np.sqrt(pca.explained_variance_)

            if self.rank > 0:
                self.modrank = self.rank
                evar = np.cumsum(pca.explained_variance_ratio_)[self.rank - 1]
                self.evar_dict[self.modrank] = evar

            elif self.rank == 0:
                util_log.debug("Finding optimal rank for %s decomposition",
                               decomp)
                evar_prev = 0
                i = 1
                for evar in np.cumsum(pca.explained_variance_ratio_):
                    self.modrank = i
                    # self.evar_list.append(evar)
                    self.evar_dict[self.modrank] = evar
                    if evar - evar_prev < 0.01:
                        self.modrank = i - 1
                        evar = evar_prev
                        break
                    evar_prev = evar
                    util_log.debug(
                        "Explained variance for first %s %s components: %s", i,
                        decomp.upper(), evar)
                    i += 1

            self.W = W[:, :self.modrank]
            self.H = H[:self.modrank, :]
        elif self.decomp == "nmf":

            if self.rank > 0:
                model = self.run_nmf_model(self.rank)
                self.modrank = self.rank

            elif self.rank == 0:
                util_log.debug("Finding optimal rank for %s decomposition",
                               decomp)
                self.evarprev = 0
                for i in range(1, self.M_run.shape[0]):
                    model = self.run_nmf_model(rank=i)
                    model_fit = model()
                    evar = model_fit.fit.evar()
                    self.modrank = i

                    if (i > 2 and evar - evarprev < 0.001):
                        model = self.run_nmf_model(rank=i - 1)
                        self.modrank = i - 1
                        break

                    self.evar_dict[self.modrank] = evar
                    evarprev = evar
                    util_log.debug(
                        "Explained variance for first %s %s components: %s", i,
                        decomp.upper(), evar)

            model_fit = model()
            self.evar_dict[self.modrank] = model_fit.fit.evar()
            self.W = model_fit.basis()
            self.H = model_fit.coef()

    # Specify NMF model
    # options can be added/modified per
    # http://nimfa.biolab.si/nimfa.methods.factorization.nmf.html
    def run_nmf_model(self, rank):
        """
        Run NMF model

        """

        prng = np.random.RandomState(self.seed)
        W_init = prng.rand(self.M_run.shape[0], rank)
        H_init = prng.rand(rank, self.M_run.shape[1])

        model = nimfa.Nmf(
            self.M_run,
            rank=rank,
            # seed=None,
            H=H_init,
            W=W_init,
            update="divergence",
            objective='div',
            n_run=1,
            max_iter=200)
        return model


class writeOutput:
    """
    Class of functions for writing the output of Helmsman.

    """

    def __init__(self, dat_paths, samples, subtypes_dict):
        self.dat_paths = dat_paths
        self.samples = samples
        self.subtypes_dict = subtypes_dict

    def writeW(self, decomp_data):
        """ write W matrix """
        num_sigs = decomp_data.W.shape[1]
        W_out = pd.DataFrame(
            data=decomp_data.W,
            index=self.samples[0],
            columns=["S" + str(i) for i in range(1, num_sigs + 1)])
        W_out.to_csv(self.dat_paths["W_path"], index_label="ID", sep="\t")

    def writeH(self, decomp_data):
        """ write H matrix """
        num_sigs = decomp_data.H.shape[0]
        H_out = pd.DataFrame(
            data=decomp_data.H,
            index=["S" + str(i) for i in range(1, num_sigs + 1)],
            columns=list(sorted(self.subtypes_dict.keys())))
        H_out.to_csv(self.dat_paths["H_path"], index_label="Sig", sep="\t")

    def writeM(self, count_matrix):
        """ write M matrix """

        count_matrix_df = pd.DataFrame(
            data=count_matrix,
            index=self.samples[0],
            columns=list(sorted(self.subtypes_dict.keys())))
        count_matrix_df.to_csv(
            self.dat_paths["M_path"], index_label="ID", sep="\t")

        freq_matrix = count_matrix / (count_matrix.sum(axis=1) + 1e-8)[:, None]

        freq_matrix_df = pd.DataFrame(
            data=freq_matrix,
            index=self.samples[0],
            columns=list(sorted(self.subtypes_dict.keys())))
        freq_matrix_df.to_csv(
            self.dat_paths["M_path_rates"], index_label="ID", sep="\t")


def writeR(package, projectdir, matrixname):
    """
    auto-generate R script

    """

    rscript_path = projectdir + "/" + "Helmsman_to_" + package + ".R"
    rscript = open(rscript_path, "w+")
    print("library(\"" + package + "\")", file=rscript)
    print("library(\"devtools\")", file=rscript)
    print("install_github(\"carjed/musigtools\")", file=rscript)
    print("library(\"musigtools\")", file=rscript)
    print(
        "mu_counts <- read.table(\"" + projectdir + "/" + matrixname +
        ".txt\", header=T, stringsAsFactors=F)",
        file=rscript)
    print("msm <- format_counts(mu_counts, \"" + package + "\")", file=rscript)
    print(
        "message(\"The mutation spectra matrix generated by Helmsman is " +
        "now formatted for use with the " + package + " package, and loaded " +
        "in a data frame named 'msm'. Please refer to the " + package +
        " documentation for help with analyzing this matrix\")",
        file=rscript)
