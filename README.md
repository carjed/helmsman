
![](assets/helmsman_logo2.png)

[![DOI](https://zenodo.org/badge/136064814.svg)](https://zenodo.org/badge/latestdoi/136064814) [![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](http://www.jedidiahcarlson.com/docs/helmsman) [![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg?style=flat)](https://opensource.org/licenses/MIT) [![Build Status](https://travis-ci.org/carjed/helmsman.svg?branch=master)](https://travis-ci.org/carjed/helmsman)

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

[logo](https://www.kisspng.com/png-computer-icons-icon-design-2740670/)
[logo font](http://www.1001fonts.com/sailor-scrawl-font.html)