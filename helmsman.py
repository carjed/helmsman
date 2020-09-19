"""
Main script for running Helmsman
"""

from __future__ import print_function
import os
import sys
import argparse
import warnings
import timeit
import random
import multiprocessing
import subprocess
import numpy as np
# from joblib import Parallel, delayed
sys.path.append(os.getcwd())
import util


def main():
    #-----------------------------------------------------------------------------
    # Initialize pre-log, get version, and process args
    #-----------------------------------------------------------------------------
    start = timeit.default_timer()

    # get latest version from github tags
    # via https://stackoverflow.com/questions/14989858
    try:
        # v_dir = os.path.dirname(os.path.realpath(__file__)) + "/.git/refs/tags"
        # files = os.listdir(v_dir)
        # files = [os.path.join(v_dir, f) for f in files] # add path to each file
        # files.sort(key=lambda x: os.path.getmtime(x))
        # version = files[-1]
        # version = os.path.basename(version)
        version = "1.5.2"
    except AttributeError:
        version = "[version not found]"

    #-----------------------------------------------------------------------------
    # Runtime control args
    #-----------------------------------------------------------------------------
    parser = argparse.ArgumentParser()

    num_cores = multiprocessing.cpu_count()
    parser.add_argument(
        "-c",
        "--cpus",
        help="number of CPUs. Must be integer value between 1 \
                            and " + str(num_cores),
        nargs='?',
        type=int,
        choices=range(1, num_cores + 1),
        metavar='INT',
        default=1)

    parser.add_argument(
        "-S",
        "--seed",
        help="random seed for NMF decomposition",
        nargs='?',
        type=int,
        metavar='INT',
        default=int(start))

    parser.add_argument(
        "-v", "--verbose", help="Enable verbose logging", action="store_true")

    parser.add_argument(
        "-V", "--version", action="version", version='%(prog)s ' + version)

    #-----------------------------------------------------------------------------
    # Input args
    #-----------------------------------------------------------------------------
    mode_opts = ["vcf", "maf", "agg", "txt"]
    parser.add_argument(
        "-M",
        "--mode",
        help="Mode for parsing input. Must be one of \
                            {" + ", ".join(mode_opts) + "}. \
                            Defaults to VCF mode.",
        nargs='?',
        type=str,
        choices=mode_opts,
        metavar='STR',
        default="vcf")

    parser.add_argument(
        "-i",
        "--input",
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

    parser.add_argument(
        "-w",
        "--rowwise",
        help="Compile mutation spectra matrix from VCF files \
                            containing non-overlapping samples.",
        action="store_true")

    parser.add_argument(
        "-f",
        "--fastafile",
        help="reference fasta file",
        nargs='?',
        type=str,
        metavar='/path/to/genome.fa',
        default="chr20.fasta.gz")

    parser.add_argument(
        "-s",
        "--samplefile",
        help="file with sample IDs to include (one per line)",
        nargs='?',
        metavar='/path/to/kept_samples.txt',
        type=str)

    parser.add_argument(
        "-g",
        "--groupvar",
        help="if --samplefile is provided with VCF input, or if \
                            input is MAF file, specify column name of the \
                            grouping variable to pool samples by. If left blank, \
                            matrix will be constructed per sample/tumor ID as usual",
        nargs='?',
        type=str,
        metavar='STR')

    parser.add_argument(
        "-H",
        "--haploid",
        help=
        "By default, Helmsman assumes diploid genotypes. For data containing \
                            haploid genotypes (e.g., male X chromosome), you must use \
                            this flag to ensure alleles are properly counted",
        action="store_true")

    parser.add_argument(
        "-u",
        "--impute",
        help="if using VCF input mode, missing genotypes \
                            (i.e., \"./.\") will be imputed as the allele \
                            frequency of the samples with non-missing genotypes",
        action="store_true")

    #-----------------------------------------------------------------------------
    # Pre-filtering args
    #-----------------------------------------------------------------------------
    parser.add_argument(
        "-C",
        "--minsnvs",
        help="minimum # of SNVs per individual to be included \
                            in analysis. Default is 0.",
        nargs='?',
        type=int,
        metavar='INT',
        default=0)

    parser.add_argument(
        "-X",
        "--maxac",
        help="maximum allele count for SNVs to keep in analysis. \
                            Defaults to 0 (all variants)",
        nargs='?',
        type=int,
        metavar='INT',
        default=0)

    #-----------------------------------------------------------------------------
    # Output args
    #-----------------------------------------------------------------------------
    parser.add_argument(
        "-p",
        "--projectdir",
        help="directory to store output files \
                            (do NOT include a trailing '/'). \
                            Defaults to " + os.getcwd() + "/helmsman_output",
        nargs='?',
        type=str,
        metavar='/path/to/project_directory',
        default="helmsman_output")

    parser.add_argument(
        "-m",
        "--matrixname",
        help="filename prefix for M matrix [without extension]",
        nargs='?',
        type=str,
        metavar='STR',
        default="subtype_count_matrix")

    package_opts = [
        "deconstructSigs", "maftools", "MutationalPatterns",
        "SomaticSignatures", "signeR", "YAPSA"
    ]
    parser.add_argument(
        "-k",
        "--package",
        help="To use the mutation spectra matrix generated by \
                            Helmsman with a specific mutation signature analysis \
                            package, this option will print out the code necessary \
                            to load the Helmsman output into R and reformat for \
                            compatibility with one of the following packages: \
                            {" + ", ".join(package_opts) + "}.",
        nargs='?',
        type=str,
        choices=package_opts,
        metavar='STR')
    #-----------------------------------------------------------------------------
    # Decomposition and outlier detection args
    #-----------------------------------------------------------------------------
    decomp_opts = ["nmf", "pca"]
    parser.add_argument(
        "-d",
        "--decomp",
        help="mode for matrix decomposition. Must be one of \
                            {" + ", ".join(decomp_opts) + "}. \
                            Defaults to 'none'.",
        nargs='?',
        type=str,
        choices=decomp_opts,
        metavar='STR')

    # rank_opts = range(2,11)
    # ro_str = str(min(rank_opts)) + " and " + str(max(rank_opts))
    parser.add_argument(
        "-r",
        "--rank",
        help="Rank for Matrix decomposition. \
                            If --decomp pca, will select first r components. \
                            Default [0] will force Helmsman to iterate through \
                            multiple ranks to find an optimal choice.",
        nargs='?',
        type=int,
        # choices=rank_opts,
        metavar='INT',
        default=0)

    motif_length_opts = [1, 3, 5, 7]
    mlo_str = ",".join(str(x) for x in motif_length_opts)

    parser.add_argument(
        "-l",
        "--length",
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

    util.util_log.setLevel(loglev)
    log = util.get_logger("helmsman", level=loglev)

    log.info("----------------------------------")
    try:
        # version = subprocess.check_output(["git",
        #                                   "describe"]).strip().decode('utf-8')
        # v_dir = os.path.dirname(os.path.realpath(__file__)) + "/.git/refs/tags"
        # files = os.listdir(v_dir)
        # files = [os.path.join(v_dir, f) for f in files] # add path to each file
        # files.sort(key=lambda x: os.path.getmtime(x))
        # version = files[-1]
        # version = os.path.basename(version)
        version = "1.4.2"
        log.info("%s %s", sys.argv[0], version)
    except AttributeError:
        version = "[version not found]"
        log.warning(version)
    log.info("----------------------------------")

    if (args.mode == "maf" and not args.groupvar):
        args.groupvar = "Tumor_Sample_Barcode"

    log.debug("Running with the following options:")
    for arg in vars(args):
        log.debug("%s : %s", arg, getattr(args, arg))

    random.seed(args.seed)
    log.info("random seed: %s", str(args.seed))

    #-----------------------------------------------------------------------------
    # Initialize project directory
    #-----------------------------------------------------------------------------
    projdir = os.path.realpath(args.projectdir)

    if not os.path.exists(args.projectdir):
        log.warning("%s does not exist--creating now", projdir)
        os.makedirs(args.projectdir)
    else:
        log.debug("All output files will be located in: %s", projdir)

    #-----------------------------------------------------------------------------
    # index subtypes
    #-----------------------------------------------------------------------------
    subtypes_dict = util.indexSubtypes(args.length)

    #-----------------------------------------------------------------------------
    # Build M matrix from inputs
    #-----------------------------------------------------------------------------
    data_in = util.processInput(args.mode, args, subtypes_dict)
    data = data_in.data
    count_matrix = data.M
    samples = np.array([data.samples], dtype=str)

    #-----------------------------------------------------------------------------
    # Drop samples from M matrix with too few SNVs
    #-----------------------------------------------------------------------------
    if args.minsnvs > 0:

        lowsnv_samples = []
        highsnv_samples = []
        i = 0
        for i in range(0, count_matrix.shape[0]):
            if sum(count_matrix[i]) < args.minsnvs:
                lowsnv_samples.append(samples.flatten()[i])
            else:
                highsnv_samples.append(samples.flatten()[i])
            i += 1

        if lowsnv_samples:
            count_matrix = count_matrix[np.sum(count_matrix, axis=1) >= args.
                                        minsnvs, ]
            samples = np.array([highsnv_samples])
            lowsnv_path = projdir + \
                "/helmsman_snvs_lt" + str(args.minsnvs) + ".txt"
            lowsnv_fh = open(lowsnv_path, "w")
            for sample in lowsnv_samples:
                lowsnv_fh.write("%s\n" % sample)
            lowsnv_fh.close()
            log.info("%s samples have fewer than %s SNVs and will be dropped",
                     len(lowsnv_samples), args.minsnvs)

    #-----------------------------------------------------------------------------
    # Get matrix decomposition and write output to files
    #-----------------------------------------------------------------------------
    paths = {
        'M_path': projdir + "/" + args.matrixname + ".txt",
        'M_path_rates': projdir + "/" + args.matrixname + "_spectra.txt",
        'W_path': projdir + "/W_components.txt",
        'H_path': projdir + "/H_loadings.txt"
    }

    dat_out = util.writeOutput(paths, samples, subtypes_dict)

    try:
        dat_out.writeM(count_matrix)
        log.debug("Spectra count matrix saved to: %s", paths['M_path'])
        log.debug("Spectra frequency matrix saved to: %s",
                  paths['M_path_rates'])
    except IOError:
        log.warning("could not write W matrix")

    if args.decomp is not None:
        decomp_data = util.DecompModel(count_matrix, args.rank, args.seed,
                                       args.decomp)

        try:
            dat_out.writeW(decomp_data)
            log.debug("W matrix saved to: %s", paths['W_path'])
        except IOError:
            log.warning("could not write W matrix")
        try:
            dat_out.writeH(decomp_data)
            log.debug("H matrix saved to: %s", paths['H_path'])
        except IOError:
            log.warning("could not write H matrix")

    #-----------------------------------------------------------------------------
    # auto-generate R script to pass data to MSA packages
    #-----------------------------------------------------------------------------
    if args.package:
        util.writeR(args.package, args.projectdir, args.matrixname)
        log_message = "To use this mutation spectra matrix" + \
            "with the {} R package, ".format(args.package) + \
            "run the following command in R: \n" + \
            "\n\t\tsource(\"{}/Helmsman_to_{}.R\")\n".format(args.projectdir, args.package)
        log.info(log_message)

    stop = timeit.default_timer()
    tottime = round(stop - start, 2)
    log.info("Total runtime: %s seconds", tottime)


main()
