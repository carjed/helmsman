#!/usr/bin/python

from __future__ import print_function
import os
import sys

sys.path.append(os.getcwd())

import shutil
import textwrap
import argparse
import warnings
import itertools
import timeit
import time
import random
import multiprocessing
import numpy as np
from joblib import Parallel, delayed
import subprocess
from distutils.dir_util import copy_tree
from util import *

###############################################################################
# Initialize pre-log, get version, and process args
###############################################################################
start = timeit.default_timer()

# get latest version from github tags
# via https://stackoverflow.com/questions/14989858
try:
    version = subprocess.check_output(["git", "describe"]).strip().decode('utf-8')
except:
    version = "[version not found]"

#-----------------------------------------------------------------------------
# Runtime control args
#-----------------------------------------------------------------------------
parser = argparse.ArgumentParser()

num_cores = multiprocessing.cpu_count()
parser.add_argument("-c", "--cpus",
                    help="number of CPUs. Must be integer value between 1 \
                        and "+str(num_cores),
                    nargs='?',
                    type=int,
                    choices=range(1,num_cores+1),
                    metavar='INT',
                    default=1)

parser.add_argument("-S", "--seed",
                    help="random seed for NMF decomposition",
                    nargs='?',
                    type=int,
                    metavar='INT',
                    default=int(start))

parser.add_argument("-v", "--verbose",
                    help="Enable verbose logging",
                    action="store_true")

parser.add_argument("-V", "--version",
                    action="version",
                    version='%(prog)s ' + version)

#-----------------------------------------------------------------------------
# Input args
#-----------------------------------------------------------------------------
mode_opts = ["vcf", "maf", "agg", "txt"]
parser.add_argument("-M", "--mode",
                    help="Mode for parsing input. Must be one of \
                        {"+", ".join(mode_opts)+ "}. \
                        Defaults to VCF mode.",
                    nargs='?',
                    type=str,
                    choices=mode_opts,
                    metavar='STR',
                    default="vcf")

parser.add_argument("-i", "--input",
                    help="In VCF mode (default) input file is a VCF \
                        or text file containing paths of multiple VCFs. \
                        Defaults to accept input from STDIN with \"--input -\". \
                        In aggregation mode, input file is a text file \
                        containing mutation subtype count matrices, \
                        or paths of multiple such matrices. \
                        In plain text mode, input file is tab-delimited text \
                        file containing 5 columns: CHR, POS, REF, ALT, ID",
                    required=True,
                    nargs='?',
                    type=str,
                    metavar='/path/to/input.vcf',
                    default=sys.stdin)

parser.add_argument("-w", "--rowwise",
                    help="Compile mutation spectra matrix from VCF files \
                        containing non-overlapping samples.",
                    action="store_true")

parser.add_argument("-f", "--fastafile",
                    help="reference fasta file",
                    nargs='?',
                    type=str,
                    metavar='/path/to/genome.fa',
                    default="chr20.fasta.gz")

parser.add_argument("-s", "--samplefile",
                    help="file with sample IDs to include (one per line)",
                    nargs='?',
                    metavar='/path/to/kept_samples.txt',
                    type=str)

parser.add_argument("-g", "--groupvar",
                    help="if --samplefile is provided with VCF input, or if \
                        input is MAF file, specify column name of the \
                        grouping variable to pool samples by. If left blank, \
                        matrix will be constructed per sample/tumor ID as usual",
                    nargs='?',
                    type=str,
                    metavar='STR')  

parser.add_argument("-u", "--impute",
                    help="if using VCF input mode, missing genotypes \
                        (i.e., \"./.\") will be imputed as the allele \
                        frequency of the samples with non-missing genotypes",
                    action="store_true")  

#-----------------------------------------------------------------------------
# Pre-filtering args
#-----------------------------------------------------------------------------
parser.add_argument("-C", "--minsnvs",
                    help="minimum # of SNVs per individual to be included \
                        in analysis. Default is 0.",
                    nargs='?',
                    type=int,
                    metavar='INT',
                    default=0)

parser.add_argument("-X", "--maxac",
                    help="maximum allele count for SNVs to keep in analysis. \
                        Defaults to 0 (all variants)",
                    nargs='?',
                    type=int,
                    metavar='INT',
                    default=0)

#-----------------------------------------------------------------------------
# Output args
#-----------------------------------------------------------------------------
parser.add_argument("-p", "--projectdir",
                    help="directory to store output files \
                        (do NOT include a trailing '/'). \
                        Defaults to " + os.getcwd() + "/helmsman_output",
                    nargs='?',
                    type=str,
                    metavar='/path/to/project_directory',
                    default="helmsman_output")

parser.add_argument("-m", "--matrixname",
                    help="filename prefix for M matrix [without extension]",
                    nargs='?',
                    type=str,
                    metavar='STR',
                    default="subtype_count_matrix")

package_opts = ["deconstructSigs", "maftools", "MutationalPatterns", 
                "SomaticSignatures", "signeR", "YAPSA"]
parser.add_argument("-k", "--package",
                    help="To use the mutation spectra matrix generated by \
                        Helmsman with a specific mutation signature analysis \
                        package, this option will print out the code necessary \
                        to load the Helmsman output into R and reformat for \
                        compatibility with one of the following packages: \
                        {"+", ".join(package_opts)+"}.",
                    nargs='?',
                    type=str,
                    choices=package_opts,
                    metavar='STR')
#-----------------------------------------------------------------------------
# Decomposition and outlier detection args
#-----------------------------------------------------------------------------
decomp_opts = ["nmf", "pca"]
parser.add_argument("-d", "--decomp", 
                    help="mode for matrix decomposition. Must be one of \
                        {"+", ".join(decomp_opts)+"}. \
                        Defaults to 'none'.",
                    nargs='?',
                    type=str,
                    choices=decomp_opts,
                    metavar='STR')

# rank_opts = range(2,11)
# ro_str = str(min(rank_opts)) + " and " + str(max(rank_opts))
parser.add_argument("-r", "--rank",
                    help="Rank for Matrix decomposition. \
                        If --decomp pca, will select first r components. \
                        Default [0] will force Helmsman to iterate through \
                        multiple ranks to find an optimal choice.",
                    nargs='?',
                    type=int,
                    # choices=rank_opts,
                    metavar='INT',
                    default=0)

motif_length_opts = [1,3,5,7]
mlo_str = ",".join(str(x) for x in motif_length_opts)

parser.add_argument("-l", "--length",
                    help="motif length. Allowed values are " + mlo_str,
                    nargs='?',
                    type=int,
                    choices=motif_length_opts,
                    metavar='INT',
                    default=3)

#-----------------------------------------------------------------------------
# initialize args and configure runtime logs
#-----------------------------------------------------------------------------
args = parser.parse_args()

# ignore warnings in sklearn 0.19.1 about covariance matrix when performing 
# outlier detection using elliptic envelope
# see https://github.com/scikit-learn/scikit-learn/issues/8811
# https://stackoverflow.com/questions/32612180
warnings.filterwarnings("ignore", category=RuntimeWarning)

if args.verbose:
    loglev = 'DEBUG' 
else:
    loglev = 'INFO'
    # ignore warning about covariance matrix not being full rank
    warnings.filterwarnings("ignore", category=UserWarning)    

util_log.setLevel(loglev)
log = getLogger("helmsman", level=loglev)

log.info("----------------------------------")
try:
    version = subprocess.check_output(["git", "describe"]).strip().decode('utf-8')
    log.info(sys.argv[0] + " " + str(version))
except:
    version = "[version not found]"
    log.warning(version)
log.info("----------------------------------")

if (args.mode == "maf" and not args.groupvar):
    args.groupvar = "Tumor_Sample_Barcode"

log.debug("Running with the following options:")
for arg in vars(args):
    log.debug(arg + ": " + str(getattr(args, arg)))

random.seed(args.seed)
log.info("random seed: " + str(args.seed))

#-----------------------------------------------------------------------------
# Initialize project directory
#-----------------------------------------------------------------------------
projdir = os.path.realpath(args.projectdir)

if not os.path.exists(args.projectdir):
    log.warn(projdir + " does not exist--creating now")
    os.makedirs(args.projectdir)
else:
    log.debug("All output files will be located in: " + projdir)

#-----------------------------------------------------------------------------
# index subtypes
#-----------------------------------------------------------------------------
subtypes_dict = indexSubtypes(args.length)

###############################################################################
# Build M matrix from inputs
###############################################################################
if args.mode == "vcf":
    if(args.input.lower().endswith(('.vcf', '.vcf.gz', '.bcf')) or
            args.input == "-"):
        par = False
        data = processVCF(args, args.input, subtypes_dict, par)
        M = data.M
        samples = np.array([data.samples], dtype=str)

    elif(args.input.lower().endswith(('.txt'))):
        par = True
        with open(args.input) as f:
            vcf_list = f.read().splitlines()

        results = Parallel(n_jobs=args.cpus) \
            (delayed(processVCF)(args, vcf, subtypes_dict, par) \
            for vcf in vcf_list)

        if args.rowwise:
            M = np.vstack(results)
            
            samples = np.array([])
            for vcf in vcf_list:
                samples = np.append(samples, getSamplesVCF(args, vcf))
            
        else:
            nrow, ncol = results[1].shape
            M = np.zeros((nrow, ncol))
    
            for M_sub in results:
                M = np.add(M, M_sub)
            samples = np.array([getSamplesVCF(args, vcf_list[1])])

elif args.mode == "maf":
    data = processMAF(args, subtypes_dict)
    M = data.M
    samples = np.array([data.samples], dtype=str)

elif args.mode == "txt":
    data = processTxt(args, subtypes_dict)
    M = data.M
    samples = np.array([data.samples], dtype=str)

elif args.mode == "agg":
    data = aggregateM(args.input, subtypes_dict)
    M = data.M
    samples = np.array([data.samples], dtype=str)

#-----------------------------------------------------------------------------
# Drop samples from M matrix with too few SNVs
#-----------------------------------------------------------------------------
if args.minsnvs > 0:

    lowsnv_samples = []
    highsnv_samples = []
    i = 0
    for row in M:
        if sum(M[i]) < args.minsnvs:
            lowsnv_samples.append(samples.flatten()[i])
        else:
            highsnv_samples.append(samples.flatten()[i])
        i += 1

    if len(lowsnv_samples) > 0:
        M = M[np.sum(M, axis=1)>=args.minsnvs,]
        samples = np.array([highsnv_samples])
        lowsnv_path = projdir + \
            "/helmsman_snvs_lt" + str(args.minsnvs) + ".txt"
        lowsnv_fh = open(lowsnv_path, "w")
        for sample in lowsnv_samples:
            lowsnv_fh.write("%s\n" % sample)
        lowsnv_fh.close()
        log.info(str(len(lowsnv_samples)) + " samples have fewer than " + 
            str(args.minsnvs) + " SNVs and will be dropped")

#-----------------------------------------------------------------------------
# Write M and M_f matrices
#-----------------------------------------------------------------------------
paths = {}
paths['M_path'] = projdir + "/" + args.matrixname + ".txt"

# M_f is the relative contribution of each subtype per sample
# adds 1e-4 to each count for error correction
# M_f = (M+1e-4)/(M.sum(axis=1))[:,None]
M_f = M/(M.sum(axis=1)+1e-8)[:,None]
paths['M_path_rates'] = projdir + "/" + args.matrixname + "_spectra.txt"

writeM(M, paths['M_path'], subtypes_dict, samples)
writeM(M_f, paths['M_path_rates'], subtypes_dict, samples)
log.debug("M matrix (spectra counts) saved to: " + paths['M_path'])
log.debug("M_f matrix (mutation spectra) saved to: " + paths['M_path_rates'])

###############################################################################
# Get matrix decomposition
###############################################################################

if args.decomp is not None:
    decomp_data = DecompModel(M_f, args.rank, args.seed, args.decomp)
    M_d = decomp_data.W
    
    # W matrix (contributions)
    paths['W_path'] = projdir + "/W_components.txt"
    writeW(decomp_data.W, paths['W_path'], samples)
    log.debug("W matrix saved to: " + paths['W_path'])
    
    # H matrix (loadings)
    paths['H_path'] = projdir + "/H_loadings.txt"
    writeH(decomp_data.H, paths['H_path'], subtypes_dict)
    log.debug("H matrix saved to: " + paths['H_path'])

###############################################################################
# auto-generate R script to pass data to MSA packages
###############################################################################
if args.package:
    writeR(args.package, args.projectdir, args.matrixname)
    log.info("To use this mutation spectra matrix with the " + args.package + 
        "R package, run the following command in R: \n" + 
        "\tsource(\"" + args.projectdir + 
        "/Helmsman_to_" + args.package + ".R\")")

###############################################################################
# Finish
###############################################################################
stop = timeit.default_timer()
tottime = round(stop - start, 2)
log.info("Total runtime: " + str(tottime) + " seconds")
