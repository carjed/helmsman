#!/usr/bin/python

# system packages
from __future__ import print_function
import os
import sys
import warnings

# ignore nuisance warnings when loading nimfa package
warnings.filterwarnings("ignore", category=UserWarning)

from logging import StreamHandler, DEBUG, getLogger as realGetLogger, Formatter
from colorama import Fore, Back, init, Style
import textwrap
import itertools
import timeit
import collections
import csv
import re

sys.path.append(os.getcwd())

# matrix+stats processing
from pandas import *
import numpy as np

# decomposition algorithms
import nimfa
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# vcf/fasta parsing
import cyvcf2 as vcf
from cyvcf2 import VCF
from pyfaidx import Fasta
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

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
            self.stream.write(
                self.colours[record.levelname] + 
                message + 
                Style.RESET_ALL
            )
            self.stream.write(getattr(self, 'terminator', '\n'))
            self.flush()
        except (KeyboardInterrupt, SystemExit):
            raise
        except:
            self.handleError(record)

###############################################################################
# configure logger
###############################################################################
class initLogger:
    def __init__(level):
        self.level = level

def getLogger(name=None, 
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

util_log = getLogger(__name__, level="DEBUG")

###############################################################################
# collapse mutation types per strand symmetry
###############################################################################
def getCategory(mu_type):
    # if re.match("^[ACGT]*$", mu_type):
    if (mu_type == "AC" or mu_type == "TG"):
        category = "T_G"
    elif (mu_type == "AG" or mu_type == "TC"):
        category = "T_C"
    elif (mu_type == "AT" or mu_type == "TA"):
        category = "T_A"
    elif (mu_type == "CA" or mu_type == "GT"):
        category = "C_A"
    elif (mu_type == "CG" or mu_type == "GC"):
        category = "C_G"
    elif (mu_type == "CT" or mu_type == "GA"):
        category = "C_T"
    else:
        category = "unknown"
    return category

###############################################################################
# query reference genome for local sequence motif
###############################################################################
def getMotif(pos, sequence):
    motif = Seq(sequence, IUPAC.unambiguous_dna)
    altmotif = motif.reverse_complement()
    central_base = (len(motif)-1)//2

    m1 = motif[central_base]
    m2 = altmotif[central_base]

    if (m1 == "C" or m1 == "T"):
        motif_a = motif
    else:
        motif_a = altmotif

    return motif_a

###############################################################################
# define k-mer mutation subtypes
###############################################################################
def indexSubtypes(motiflength):
    categories = ["T_G", "T_C", "T_A", "C_G", "C_T", "C_A"]
    bases = ["A", "C", "G", "T"]
    flank = (motiflength-1)//2

    if motiflength > 1:
        kmers = itertools.product(bases, repeat=motiflength-1)

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
        extr = list(np.repeat(ext,3))
        subtypes_list = [m+n for m,n in zip(categories,extr)]

    i = 0
    subtypes_dict = {}
    for subtype in sorted(subtypes_list):
        subtypes_dict[subtype] = i
        i += 1

    util_log.debug(str(len(subtypes_dict.keys())) + " " + 
        str(motiflength) + "-mer subtypes indexed")

    return subtypes_dict

###############################################################################
# Build dictionary with sample ID as key, group ID as value
###############################################################################
def indexGroups(samplefile, groupvar):
    sg_dict = {}
    
    f = open(samplefile, 'r', encoding = "utf-8")
    reader = csv.DictReader(f, delimiter='\t')

    for row in reader:
        sg_dict[row['ID']] = row[groupvar]
    
    return sg_dict

###############################################################################
# get list of samples to keep if samplefile supplied
###############################################################################
def parseSampleFile(samplefile):
    # f = open(args.input, 'r', encoding = "ISO-8859-1")
    f = open(samplefile, 'r', encoding = "utf-8")
    reader = csv.DictReader(f, delimiter='\t')
    keep_samples = []
    for row in reader:
        keep_samples.append(row['ID'])
        
    return keep_samples
    
###############################################################################
# get samples from VCF file
###############################################################################
def getSamplesVCF(args, inputvcf):
    
    if args.samplefile:
        keep_samples = parseSampleFile(args.samplefile)
        vcf_reader = VCF(inputvcf,
            mode='rb', gts012=True, lazy=True, samples=keep_samples)
    else:
        vcf_reader = VCF(inputvcf,
            mode='rb', gts012=True, lazy=True)

    if (args.samplefile and args.groupvar):
        all_samples = vcf_reader.samples
        samples = indexGroups(args.samplefile, args.groupvar)
    else:
        samples = vcf_reader.samples

    return samples

###############################################################################
# Main function for parsing VCF
###############################################################################
def processVCF(args, inputvcf, subtypes_dict, par):

    # initialize reference genome
    fasta_reader = Fasta(args.fastafile, read_ahead=1000000)

    # initialize vcf reader
    if args.samplefile:
        keep_samples = parseSampleFile(args.samplefile)

        vcf_reader = VCF(inputvcf,
            mode='rb', gts012=True, lazy=True, samples=keep_samples)
    else:
        vcf_reader = VCF(inputvcf,
            mode='rb', gts012=True, lazy=True)

    nbp = (args.length-1)//2

    # index samples
    if (args.samplefile and args.groupvar):
        all_samples = vcf_reader.samples

        sg_dict = indexGroups(args.samplefile, args.groupvar)
        samples = sorted(list(set(sg_dict.values())))

        util_log.debug(str(len(all_samples)) + " samples will be pooled into " +
            str(len(samples)) + " groups: " +  ",".join(samples))
    else:
        samples = vcf_reader.samples

    samples_dict = {}
    for i in range(len(samples)):
        samples_dict[samples[i]] = i

    # Query records in VCF and build matrix
    M = np.zeros((len(samples), len(subtypes_dict)))
    numsites_keep = 0
    numsites_skip = 0
    chrseq = '0'

    for record in vcf_reader:

        # Filter by SNP status, # alt alleles, and FILTER column
        if (not record.is_snp or 
            len(record.ALT) != 1 or 
            record.FILTER is not None):
            numsites_skip += 1
            continue

        # Filter by allele count
        if (record.INFO['AC'] > args.maxac > 0):
            numsites_skip += 1
            continue

        # check and update chromosome sequence
        if record.CHROM != chrseq:
            sequence = fasta_reader[record.CHROM]
            chrseq = record.CHROM

        lseq = sequence[record.POS-(nbp+1):record.POS+nbp].seq

        mu_type = record.REF + str(record.ALT[0])
        category = getCategory(mu_type)
        motif_a = getMotif(record.POS, lseq)
        subtype = str(category + "." + motif_a)

        if subtype not in subtypes_dict:
            numsites_skip += 1
            continue
        
        st = subtypes_dict[subtype]

        # currently only works with singletons--
        if (args.samplefile and args.groupvar):

            if record.gt_types.sum() == 0:
                numsites_skip += 1
                continue
            
            carrier = all_samples[record.gt_types.tolist().index(1)]
            if carrier not in sg_dict:
                numsites_skip += 1
                continue
            
            sample_gp = sg_dict[carrier]
            ind = samples.index(sample_gp)
            M[ind,st] += 1
            numsites_keep += 1

        else:
            gt_new = record.gt_types
            gt_new[gt_new == 3] = 0
            M[:,st] = M[:,st]+gt_new
            numsites_keep += 1

        if (numsites_keep%1000000 != 0): continue
        util_log.debug(inputvcf + ": " + str(numsites_keep) + " sites counted")

    util_log.info(inputvcf + ": " + str(numsites_keep) + " sites counted")
    util_log.info(inputvcf + ": " + str(numsites_skip) + " sites skipped")

    out = collections.namedtuple('Out', ['M', 'samples'])(M, samples)

    if par:
        return M
    else:
        return out

###############################################################################
# process MAF files
###############################################################################
def processMAF(args, subtypes_dict):
    
    fasta_reader = Fasta(args.fastafile, read_ahead=1000000)
    
    nbp = (args.length-1)//2
    samples_dict = {}

    # M = np.zeros((len(samples), len(subtypes_dict)))
    numsites_keep = 0
    numsites_skip = 0
    chrseq = '0'

    f = open(args.input, 'r', encoding = "ISO-8859-1")

    reader = csv.DictReader(filter(lambda row: row[0]!='#', f), delimiter='\t')
    counter = 0
    for row in reader:

        if(row['Variant_Type'] != "SNP"): continue
            
        pos = int(row['Start_position'])
        ref = row['Reference_Allele']
        alt = row['Tumor_Seq_Allele2']
        sample = row[args.groupvar]
        
        if row['Chromosome'] != chrseq:
            sequence = fasta_reader[row['Chromosome']]
            chrseq = row['Chromosome']
        
        counter += 1
        mu_type = ref + alt
        category = getCategory(mu_type)
        lseq = sequence[pos-(nbp+1):pos+nbp].seq
        
        motif_a = getMotif(pos, lseq)
        subtype = str(category + "." + motif_a)
        st = subtypes_dict[subtype]

        if sample not in samples_dict:
            samples_dict[sample] = {}

        if subtype not in samples_dict[sample]:
            samples_dict[sample][subtype] = 1
        else:
            samples_dict[sample][subtype] += 1

        if (counter%1000 != 0): continue
        util_log.debug(args.input + ": " + str(counter) + " sites counted")

    M = DataFrame(samples_dict).T.fillna(0).values
    samples = sorted(samples_dict)

    out = collections.namedtuple('Out', ['M', 'samples'])(M, samples)
    return out

###############################################################################
# process tab-delimited text file, containing the following columns:
# CHR    POS    REF    ALT    SAMPLE_ID
###############################################################################
def processTxt(args, subtypes_dict):

    fasta_reader = Fasta(args.fastafile, read_ahead=1000000)

    nbp = (args.length-1)//2
    samples_dict = {}

    numsites_keep = 0
    numsites_skip = 0
    chrseq = '0'

    with open(args.input, 'r') as f:
        reader = csv.reader(f, delimiter='\t')

        for row in reader:
            chrom = row[0]
            pos = int(row[1])
            ref = row[2]
            alt = row[3]
            sample = row[4]

            if chrom != chrseq:
                sequence = fasta_reader[chrom]
                chrseq = chrom

            if(len(alt) == 1 and len(ref)==1):
                mu_type = ref + alt
                category = getCategory(mu_type)
                if nbp > 0:
                    lseq = sequence[pos-(nbp+1):pos+nbp].seq
                else:
                    lseq = sequence[pos-1].seq
                    # eprint("lseq:", lseq)
                motif_a = getMotif(pos, lseq)
                subtype = str(category + "-" + motif_a)
                st = subtypes_dict[subtype]

                if sample not in samples_dict:
                    samples_dict[sample] = {}

                if subtype not in samples_dict[sample]:
                    samples_dict[sample][subtype] = 1
                else:
                    samples_dict[sample][subtype] += 1

        M = DataFrame(samples_dict).T.fillna(0).values
        samples = sorted(samples_dict)

    out = collections.namedtuple('Out', ['M', 'samples'])(M, samples)
    return out

###############################################################################
# get samples from input M matrix when using aggregation mode
###############################################################################
def getSamples(fh):
    samples = np.loadtxt(fh,
        dtype='S20',
        skiprows=1,
        delimiter='\t',
        usecols=(0,))
        
    util_log.debug(fh + " contains " + str(len(samples)) + " samples")

    return samples

###############################################################################
# aggregate M matrices from list of input files
###############################################################################
def aggregateM(inputM, subtypes_dict):
    colnames = ["ID"]
    M_colnames = colnames + list(sorted(subtypes_dict.keys()))
    colrange = range(1,len(M_colnames))

    if (inputM.lower().endswith('m_samples.txt') or 
            inputM.lower().endswith('m_regions.txt')):
        with open(inputM) as f:
            file_list = f.read().splitlines()

        # M output by sample
        if inputM.lower().endswith('m_samples.txt'):

            M_out = np.array([M_colnames])

            for mfile in file_list:
                samples = getSamples(mfile)

                M_it = np.loadtxt(mfile, skiprows=1, usecols=colrange)
                M_it = np.concatenate((np.array([samples]).T, M_it), axis=1)
                M_out = np.concatenate((M_out, M_it), axis=0)

            M = np.delete(M_out, 0, 0)
            M = np.delete(M, 0, 1)
            M = M.astype(np.float)

        # M output by region
        elif inputM.lower().endswith('m_regions.txt'):
            samples = getSamples(file_list[0])

            M_out = np.zeros((len(samples), len(M_colnames)-1))
            for mfile in file_list:
                M_it = np.loadtxt(mfile, skiprows=1, usecols=colrange)
                M_out = np.add(M_out, M_it)

            M = M_out.astype(np.float)

    else:
        samples = getSamples(inputM)
        M = np.loadtxt(inputM, skiprows=1, usecols=colrange)
        M = M.astype(np.float)

    out = collections.namedtuple('Out', ['M', 'samples'])(M, samples)
    return out

###############################################################################
# Class for fitting PCA or NMF models
###############################################################################
class DecompModel:
    def __init__(self, M_run, rank, seed, decomp):
        self.M_run = M_run
        self.rank = rank
        self.seed = seed
        self.decomp = decomp
        
        self.evar_dict = {}
        
        if self.decomp == "pca":
            # standarize input matrix
            X_std = StandardScaler().fit_transform(self.M_run)
            
            # run PCA
            pca = PCA(n_components = self.M_run.shape[1])
            W = pca.fit_transform(X_std)
            H = pca.components_.T * np.sqrt(pca.explained_variance_)
            
            if self.rank > 0:
                self.modrank = self.rank
                evar = np.cumsum(pca.explained_variance_ratio_)[self.rank-1]
                self.evar_dict[self.modrank] = evar

            elif self.rank == 0:
                util_log.debug("Finding optimal rank for " + decomp + " decomposition")
                evar_prev = 0
                i = 1
                for evar in np.cumsum(pca.explained_variance_ratio_):
                    self.modrank = i
                    # self.evar_list.append(evar)
                    self.evar_dict[self.modrank] = evar
                    if evar - evar_prev < 0.01:
                        self.modrank = i-1
                        evar = evar_prev
                        break
                    evar_prev = evar
                    util_log.debug("Explained variance for first " + 
                        str(i) + " " + decomp.upper() + " components: " + 
                        str(evar))
                    i += 1
                
            self.W = W[:,:self.modrank]
            self.H = H[:self.modrank,:]
        elif self.decomp == "nmf":
        
            if self.rank > 0:
                model = self.NMFmod(self.rank)
                self.modrank = self.rank
                
            elif self.rank == 0:
                util_log.debug("Finding optimal rank for " + decomp + " decomposition")
                self.evarprev = 0
                for i in range(1,6):
                    model = self.NMFmod(rank=i)
                    model_fit = model()
                    evar = model_fit.fit.evar()
                    self.modrank = i
            
                    if(i > 2 and evar - evarprev < 0.001):
                        model = self.NMFmod(rank=i-1)
                        self.modrank = i-1
                        break
                    
                    self.evar_dict[self.modrank] = evar
                    evarprev = evar
                    util_log.debug("Explained variance for first " + 
                        str(i) + " " + decomp.upper() + " components: " + 
                        str(evar))
            
            model_fit = model()
            self.evar_dict[self.modrank] = model_fit.fit.evar()
            self.W = model_fit.basis()
            self.H = model_fit.coef()
        
    # Specify NMF model
    # options can be added/modified per 
    # http://nimfa.biolab.si/nimfa.methods.factorization.nmf.html  
    def NMFmod(self, rank):
    
        prng = np.random.RandomState(self.seed)
        W_init = prng.rand(self.M_run.shape[0], rank)
        H_init = prng.rand(rank, self.M_run.shape[1])
        
        model = nimfa.Nmf(self.M_run,
            rank=rank,
            # seed=None,
            H=H_init,
            W=W_init,
            update="divergence",
            objective='div',
            n_run=1,
            max_iter=200)
        return model

###############################################################################
# write M matrix
###############################################################################
def writeM(M, M_path, subtypes_dict, samples):

    M_out = DataFrame(data=M,
                index=samples[0],
                columns=list(sorted(subtypes_dict.keys())))

    M_out.to_csv(M_path, index_label="ID", sep="\t")

###############################################################################
# write W matrix
###############################################################################
def writeW(W, W_path, samples):
    
    num_samples, num_sigs = W.shape
    W_out = DataFrame(data=W,
                index=samples[0],
                columns=["S" + str(i) for i in range(1,num_sigs+1)])
    
    W_out.to_csv(W_path, index_label="ID", sep="\t")

###############################################################################
# write H matrix
###############################################################################
def writeH(H, H_path, subtypes_dict):
    
    num_sigs, num_subtypes = H.shape
    H_out = DataFrame(data=H,
                index=["S" + str(i) for i in range(1,num_sigs+1)],
                columns=list(sorted(subtypes_dict.keys())))

    H_out.to_csv(H_path, index_label="Sig", sep="\t")

###############################################################################
# auto-generate R script
###############################################################################
def writeR(package, projectdir, matrixname):
    rscript_path = projectdir + "/" + "Helmsman_to_" + package + ".R"
    rscript = open(rscript_path, "w+")
    print("library(\"" + package + "\")", file=rscript)
    print("library(\"devtools\")", file=rscript)
    print("install_github(\"carjed/musigtools\")", file=rscript)
    print("library(\"musigtools\")", file=rscript)
    print("mu_counts <- read.table(\"" + 
        projectdir + "/" + 
        matrixname + ".txt\", header=T, stringsAsFactors=F)", file=rscript)
    print("msm <- format_counts(mu_counts, \"" + package + "\")", file=rscript)
    print("message(\"The mutation spectra matrix generated by Helmsman is " + 
        "now formatted for use with the " + package + " package, and loaded " + 
        "in a data frame named 'msm'. Please refer to the " + package + 
        " documentation for help with analyzing this matrix\")", 
            file=rscript)

