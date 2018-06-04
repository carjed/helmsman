#!/usr/bin/python

# system packages
from __future__ import print_function
import os
import sys
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
from scipy.stats import chisquare

# decomposition algorithms
import nimfa
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# outlier detection algorithms
from sklearn.neighbors import LocalOutlierFactor
from sklearn.covariance import EllipticEnvelope
from sklearn.ensemble import IsolationForest

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
    if re.match("^[ACGT]*$", mu_type):
        if (mu_type == "AC" or mu_type == "TG"):
            category = "A_C"
        if (mu_type == "AG" or mu_type == "TC"):
            category = "A_G"
        if (mu_type == "AT" or mu_type == "TA"):
            category = "A_T"
        if (mu_type == "CA" or mu_type == "GT"):
            category = "C_A"
        if (mu_type == "CG" or mu_type == "GC"):
            category = "C_G"
        if (mu_type == "CT" or mu_type == "GA"):
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

    if m1 < m2:
        motif_a = motif
    else:
        motif_a = altmotif

    return motif_a

###############################################################################
# define k-mer mutation subtypes
###############################################################################
def indexSubtypes(motiflength):
    categories = ["A_C", "A_G", "A_T", "C_G", "C_T", "C_A"]
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
        ext = [".A", ".C"]
        extr = list(np.repeat(ext,3))
        subtypes_list = [m+n for m,n in zip(categories,extr)]

    i = 0
    subtypes_dict = {}
    for subtype in sorted(subtypes_list):
        subtypes_dict[subtype] = i
        i += 1
        util_log.debug("subtype " + str(i) + " of " + 
            str(len(subtypes_dict.keys())) + " indexed: " + subtype)

    return subtypes_dict

###############################################################################
# Build dictionary with sample ID as key, group ID as value
###############################################################################
def indexGroups(groupfile):
    sg_dict = {}
    with open(groupfile) as sg_file:
        for line in sg_file:
           (key, val) = line.split()
           sg_dict[key] = val

    samples = sorted(list(set(sg_dict.values())))
    return samples

###############################################################################
# get samples from VCF file
###############################################################################
def getSamplesVCF(args, inputvcf):
    if args.samplefile:
        with open(args.samplefile) as f:
            keep_samples = f.read().splitlines()

        vcf_reader = VCF(inputvcf,
            mode='rb', gts012=True, lazy=True, samples=keep_samples)
        # vcf_reader.set_samples(keep_samples) # <- set_samples() subsets VCF
    else:
        vcf_reader = VCF(inputvcf,
            mode='rb', gts012=True, lazy=True)

    if args.groupfile:
        all_samples = vcf_reader.samples
        samples = indexGroups(args.groupfile)
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
        with open(args.samplefile) as f:
            keep_samples = f.read().splitlines()
        util_log.debug("VCF will be subset to " +
            str(len(keep_samples)) + "samples in " +
            args.samplefile)

        vcf_reader = VCF(inputvcf,
            mode='rb', gts012=True, lazy=True, samples=keep_samples)
        # vcf_reader.set_samples(keep_samples) # <- set_samples() subsets VCF
    else:
        vcf_reader = VCF(inputvcf,
            mode='rb', gts012=True, lazy=True)

    nbp = (args.length-1)//2

    # index samples
    if args.groupfile:
        all_samples = vcf_reader.samples
        samples = indexGroups(args.groupfile)
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

    batchit = 0
    sample_batch = []
    subtype_batch = []

    for record in vcf_reader:
        # debug--testing performance for triallelic sites
        # if(record.POS==91628): # triallelic site
        # if(record.POS==63549):
        #     eprint(acval)
        #     eprint(record.gt_types.tolist().index(1))

        # Filter by allele count, SNP status, and FILTER column
        # if len(record.ALT[0])==1:
        if record.is_snp and len(record.ALT)==1:
            # eprint("SNP check: PASS")
            acval = record.INFO['AC']
#             eprint(record.POS, acval)

            if ((acval<=args.maxac or args.maxac==0) and record.FILTER is None):
                # eprint(record.CHROM, record.POS, record.REF, record.ALT[0],
                    # acval, record.FILTER)

                # check and update chromosome sequence
                if record.CHROM != chrseq:
                    sequence = fasta_reader[record.CHROM]
                    chrseq = record.CHROM

                if nbp > 0:
                    lseq = sequence[record.POS-(nbp+1):record.POS+nbp].seq
                else:
                    lseq = sequence[record.POS-1].seq

                mu_type = record.REF + str(record.ALT[0])
                category = getCategory(mu_type)
                motif_a = getMotif(record.POS, lseq)
                subtype = str(category + "." + motif_a)

                if subtype in subtypes_dict:
                    st = subtypes_dict[subtype]

                    if args.groupfile:
                        sample = all_samples[record.gt_types.tolist().index(1)]

                        if sample in sg_dict:
                            sample_gp = sg_dict[sample]
                            ind = samples.index(sample_gp)
                            M[ind,st] += 1
                    else:
                        gt_new = record.gt_types
                        gt_new[gt_new == 3] = 0
                        M[:,st] = M[:,st]+gt_new

                    numsites_keep += 1

                else:
                    numsites_skip += 1

                if (numsites_keep%1000000 == 0):
                    util_log.debug(inputvcf + ": " + 
                        str(numsites_keep) + " sites counted")
                    # util_log.debug(str(numsites_skip) + " sites skipped")

            else:
                numsites_skip += 1

    util_log.info(inputvcf + ": " + 
        str(numsites_keep) + " sites counted")
    util_log.info(inputvcf + ": " + 
        str(numsites_skip) + " sites skipped")

    out = collections.namedtuple('Out', ['M', 'samples'])(M, samples)

    if par:
        return M
    else:
        return out

###############################################################################
# process tab-delimited text file, containing the following columns:
# CHR    POS    REF    ALT    SAMPLE_ID
###############################################################################
def processTxt(args, subtypes_dict):

    fasta_reader = Fasta(args.fastafile, read_ahead=1000000)

    nbp = (args.length-1)//2
    samples_dict = {}

    # M = np.zeros((len(samples), len(subtypes_dict)))
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

    if inputM.lower().endswith('nmf_m_spectra.txt'):
        samples = getSamples(inputM)
        M = np.loadtxt(inputM, skiprows=1, usecols=colrange)
        M = M.astype(np.float)
    else:
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
# Generate keep/drop lists
###############################################################################
class DetectOutliers:
    def __init__(self, M, samples, filtermode, threshold, projdir, seed):

        # outlier detection
        clf = LocalOutlierFactor(
            n_neighbors=20, 
            contamination=threshold)
        y_pred = clf.fit_predict(M)
        
        cee = EllipticEnvelope(
            contamination=threshold,
            random_state=seed)
        cee.fit(M)
        scores_pred = cee.decision_function(M)
        y_pred2 = cee.predict(M)
        
        cif = IsolationForest(
            contamination=threshold,
            random_state=seed)
        cif.fit(M)
        scores_pred = cif.decision_function(M)
        y_pred3 = cif.predict(M)
        
        outlier_methods = ["lof", "ee", "if"]
        ol_df = DataFrame(np.column_stack((y_pred, y_pred2, y_pred3)),
                   index=samples[0].tolist(),
                   columns=outlier_methods)
    
        keep_samples, drop_samples, drop_indices = ([] for i in range(3))
    
        omnibus_methods = ["any", "any2", "all"]
        if filtermode in omnibus_methods:
            dft = ol_df.sum(axis=1)
            dft = DataFrame(dft)
            if filtermode == "any":
                drop_samples = dft[dft[0] != 3].index.values.tolist()
                keep_samples = dft[dft[0] == 3].index.values.tolist()
            elif filtermode == "any2":
                drop_samples = dft[dft[0] <= -1].index.values.tolist()
                keep_samples = dft[dft[0] > -1].index.values.tolist()
            elif filtermode == "all":
                drop_samples = dft[dft[0] == -3].index.values.tolist()
                keep_samples = dft[dft[0] != -3].index.values.tolist()
            
        elif filtermode in outlier_methods:
            drop_samples = ol_df[ol_df[filtermode] == -1].index.values.tolist()
            keep_samples = ol_df[ol_df[filtermode] == 1].index.values.tolist()
            
        drop_bool = np.isin(samples[0], drop_samples)
        drop_indices = np.where(drop_bool)[0].tolist()
        
        self.keep = keep_samples
        self.drop = drop_samples
        self.drop_indices = drop_indices

###############################################################################
# write yaml config for diagnostic reports
###############################################################################
def writeReportConfig(paths, projdir, args):
    yaml_path = projdir + "/config.yaml"
    yaml = open(yaml_path, "w+")
    print("# Config file for doomsayer_diagnostics.r", file=yaml)
    
    for key in paths.keys():
        print(key + ": " + paths[key], file=yaml)

    print("staticplots: " + str(args.staticplots).lower(), file=yaml)

###############################################################################
# filter VCF input by kept samples
###############################################################################
def filterVCF(inputvcf, keep_samples):
    vcf = VCF(inputvcf, samples=keep_samples, mode='rb')

    print(vcf.raw_header.rstrip())
    for v in vcf:
        v.INFO['AC'] = str(v.num_het + v.num_hom_alt*2)

        if int(v.INFO['AC']) > 0:
            v.INFO['NS'] = str(v.num_called)
            v.INFO['AN'] = str(2*v.num_called)
            v.INFO['DP'] = str(np.sum(v.format('DP')))
            print(str(v).rstrip())

###############################################################################
# filter txt input by kept samples
###############################################################################
def filterTXT(inputtxt, keep_samples):
    with open(inputtxt, 'r') as f:
        reader = csv.reader(f, delimiter='\t')

        for row in reader:
            chrom = row[0]
            pos = row[1]
            ref = row[2]
            alt = row[3]
            sample = row[4]

            if sample in keep_samples:
                print("\t".join(row))
