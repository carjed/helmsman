
![](assets/helmsman_logo2.png)

[![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](http://www.jedidiahcarlson.com/docs/helmsman) [![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg?style=flat)](https://opensource.org/licenses/MIT) [![Build Status](https://travis-ci.org/carjed/helmsman.svg?branch=master)](https://travis-ci.org/carjed/helmsman)

> O be swift—<br />
we have always known you wanted us.

_(from 'The Helmsman' by Hilda Doolittle [1886 - 1961])_

------------------------------------

# Introduction

**Helmsman** is a utility for rapidly and efficiently generating mutation spectra matrices from massive next-generation sequencing datasets. See [Alexandrov et al., _Cell Reports_, 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3588146/) for a detailed explanation of why mutation signature analysis methods are useful, and how they work.

Currently, the majority of mutation signature analysis methods are implemented as R packages. To generate the mutation spectra matrix, these packages must read the entire dataset, containing information about every individual mutation, into memory. This is generally not a problem for datasets containing only a few thousand SNVs or a few dozen samples, but for much larger datasets, it is extremely easy to exceed the physical memory capacity of your machine (not to mention that even for small files, these functions can be quite slow to run--sometimes over 30 minutes for a 3Mb file!).

Other mutation signature analysis tools do not even provide these convenience functions, and leave it to the user to coerce their data into program-specific formats. Not only is this inconvenient, but it impedes users' ability to port their data between tools and take advantage of the unique features provided by different packages.

Helmsman aims to alleviate these performance barriers and lack of standardization. Helmsman was initially written to evaluate patterns of variation in massive whole-genome datasets, containing tens or hundreds of millions of SNVs observed in tens of thousands of individuals. Generating the mutation spectra matrix for such data is impossible with R-based implementations, but Helmsman 

Helmsman includes several other convenient features, including:

- ability to pool samples together when generating the mutation spectra matrix
- aggregate data from multiple VCF files
- bare-bones (and really fast!) non-negative matrix factorization, to extract signatures without even relying on external packages
- able to run in parallel

------------------------------------

# Setup

## Using Conda (recommended)

The easiest way to start using Helmsman is to create a Conda environment, based on the dependencies specified in the `env.yml` file:

```{sh}
git clone https://github.com/carjed/helmsman.git
cd helmsman

conda env create -n helmsman -f env.yml
source activate helmsman
```

## Using pip

If you do not have Conda on your system, the prerequisites for Doomsayer can also be installed using `pip`:

```{sh}
git clone https://github.com/carjed/doomsayer.git
cd doomsayer

pip install -r pip_reqs.txt
```

# Quick Start

Suppose we have a Variant Call format (VCF) file named `input.vcf`, containing the genotypes of N individuals at each somatic mutation identified. You will also need the corresponding reference genome. With the following command, With the following command, Helmsman will parse the VCF file and write a file under `/path/to/output/` containing the Nx96 mutation spectra matrix:

```{sh}
python helmsman.py --input /path/to/input.vcf --fastafile /path/to/reference_genome.fasta --projectdir /path/to/output/
```

<!-- ### Citation
If you use Helmsman in your research, please cite our [paper](#) (citation pending). -->

# Usage

```
usage: helmsman.py [-h] [-c [INT]] [-S [INT]] [-v] [-V] [-M [STR]] -i
                   [/path/to/input.vcf] [-w] [-f [/path/to/genome.fa]]
                   [-g [/path/to/sample_batches.txt]]
                   [-s [/path/to/kept_samples.txt]] [-C [INT]] [-X [INT]]
                   [-p [/path/to/project_directory]] [-m [STR]] [-d [STR]]
                   [-r [INT]] [-t [FLOAT]] [-l [INT]]

optional arguments:
  -h, --help            show this help message and exit
  -c [INT], --cpus [INT]
                        number of CPUs. Must be integer value between 1 and 10
  -S [INT], --seed [INT]
                        random seed for NMF and outlier detection
  -v, --verbose         Enable verbose logging
  -V, --version         show program's version number and exit
  -M [STR], --mode [STR]
                        Mode for parsing input. Must be one of {vcf, agg,
                        txt}. Defaults to VCF mode.
  -i [/path/to/input.vcf], --input [/path/to/input.vcf]
                        In VCF mode (default) input file is a VCF or text file
                        containing paths of multiple VCFs. Defaults to accept
                        input from STDIN with "--input -". In aggregation
                        mode, input file is a text file containing mutation
                        subtype count matrices, or paths of multiple such
                        matrices. In plain text mode, input file is tab-
                        delimited text file containing 5 columns: CHR, POS,
                        REF, ALT, ID
  -w, --rowwise         Compile mutation spectra matrix from VCF files
                        containing non-overlapping samples.
  -f [/path/to/genome.fa], --fastafile [/path/to/genome.fa]
                        reference fasta file
  -g [/path/to/sample_batches.txt], --groupfile [/path/to/sample_batches.txt]
                        two-column tab-delimited file containing sample IDs
                        (column 1) and group membership (column 2) for pooled
                        analysis
  -s [/path/to/kept_samples.txt], --samplefile [/path/to/kept_samples.txt]
                        file with sample IDs to include (one per line)
  -C [INT], --minsnvs [INT]
                        minimum # of SNVs per individual to be included in
                        analysis. Default is 0.
  -X [INT], --maxac [INT]
                        maximum allele count for SNVs to keep in analysis.
                        Defaults to 0 (all variants)
  -p [/path/to/project_directory], --projectdir [/path/to/project_directory]
                        directory to store output files (do NOT include a
                        trailing '/'). Defaults to /mnt/norbert/home/jedidiah/
                        projects/helmsman/doomsayer_output
  -m [STR], --matrixname [STR]
                        filename prefix for M matrix [without extension]
  -d [STR], --decomp [STR]
                        mode for matrix decomposition. Must be one of {nmf,
                        pca}. Defaults to pca.
  -r [INT], --rank [INT]
                        Rank for Matrix decomposition. If --decomp pca, will
                        select first r components. Default [0] will force
                        Doomsayer to iterate through multiple ranks to find an
                        optimal choice.
  -t [FLOAT], --threshold [FLOAT]
                        threshold for fraction of potential outliers
  -l [INT], --length [INT]
                        motif length. Allowed values are 1,3,5,7
```