<img src="https://raw.githubusercontent.com/carjed/helmsman/master/assets/Helmsman_white_bg.png" width="400" height="400">

[![DOI](https://zenodo.org/badge/136064814.svg)](https://zenodo.org/badge/latestdoi/136064814) [![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](http://www.jedidiahcarlson.com/docs/helmsman) [![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg?style=flat)](https://opensource.org/licenses/MIT) [![Binder](https://img.shields.io/badge/launch-binder-d06681.svg?style=flat)](https://mybinder.org/v2/gh/carjed/helmsman/master) [![Build Status](https://travis-ci.org/carjed/helmsman.svg?branch=master)](https://travis-ci.org/carjed/helmsman)

> O be swift—<br />
we have always known you wanted us.

_(from 'The Helmsman' by Hilda Doolittle [1886 - 1961])_

------------------------------------

# Introduction

Helmsman is a utility for rapidly and efficiently generating mutation spectra matrices from massive next-generation sequencing datasets, for use in a wide range of mutation signature analysis tools. See [Alexandrov et al., _Cell Reports_, 2013](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3588146/) for a detailed explanation of how mutation signature analysis methods work, and why they are useful in studying cancer genomes.

Currently, the majority of mutation signature analysis methods are implemented as R packages. To generate the mutation spectra matrix, these packages must read the entire dataset (containing information about every individual mutation), into memory. This is generally not a problem for datasets containing only a few thousand SNVs or a few dozen samples, but for much larger datasets, it is extremely easy to exceed the physical memory capacity of your machine. Depending on the dimensions of the data, even very small files can take a very long time to process in these R packages--sometimes over 30 minutes for a 3Mb file!

Other mutation signature analysis tools do not even provide these convenience functions, and leave it to the user to coerce their data into program-specific formats. Not only is this inconvenient, but it impedes users' ability to port their data between tools and take advantage of the unique features provided by different packages.

Helmsman aims to alleviate these performance barriers and lack of standardization. Helmsman was initially written to evaluate patterns of variation in massive whole-genome datasets, containing tens or hundreds of millions of SNVs observed in tens of thousands of individuals. Generating the mutation spectra matrix for such data is virtually impossible with R-based implementations, but Helmsman provides a very fast and scalable solution for analyzing these datasets.

If you plan to use the output of Helmsman in other mutation signature analysis packages, Helmsman can automatically generate a small R script with all the code necessary to read the mutation spectra matrix and format it for compatibility with existing tools, using functions from the [musigtools](https://github.com/carjed/musigtools) package.

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

If you do not have Conda on your system, the prerequisites for Helmsman can also be installed with `pip` inside of a python3 virtual environment (assuming you have [`virtualenv`](https://virtualenv.pypa.io/en/latest/) installed):

```{sh}
git clone https://github.com/carjed/helmsman.git
cd helmsman

virtualenv -p python3 helmsman_env
source helmsman_env/bin/activate

pip install -r pip_reqs.txt
```

It is also possible to forgo the virtual environment setup and use `pip` to install the necessary dependencies in your global `site-packages` directory, but this is not recommended as doing so may cause dependency conflicts between Helmsman and other programs/packages.

## Docker

For more flexible deployment options, Helmsman is available as a Docker container. The following command will pull and run the preconfigured image from the [Docker Hub](https://hub.docker.com/):

```
docker run -d --name helmsman \
  -v /path/to/local/data:/data \ # map directory containing input data
  -p 8888:8888 \ # expose jupyter notebook on port 8888
  start-notebook.sh --NotebookApp.token='' \ # start with token disabled
  carjed/helmsman
```

You may also clone this repository and build the dockerfile locally, using the following commands:

```{sh}
git clone https://github.com/carjed/helmsman.git
cd helmsman

docker build -t latest --force-rm .

docker run -d --name helmsman \
  -p 8888:8888 \
  start-notebook.sh --NotebookApp.token='' \
  helmsman
```

# Quick Start

Suppose we have a Variant Call format (VCF) file named `input.vcf`, containing the genotypes of N individuals at each somatic mutation identified. You will also need the corresponding reference genome. With the following command, With the following command, Helmsman will parse the VCF file and write a file under `/path/to/output/` containing the Nx96 mutation spectra matrix:

```{sh}
python helmsman.py --input /path/to/input.vcf --fastafile /path/to/reference_genome.fasta --projectdir /path/to/output/
```

# Citation
If you use Helmsman in your research, please cite our [paper](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-5264-y) published in BMC Genomics:

> Carlson J, Li JZ, Z&ouml;llner S. Helmsman: fast and efficient mutation signature analysis for massive sequencing datasets. *BMC Genomics.* 2018;19:845. [`doi:10.1186/s12864-018-5264-y`](http://dx.doi.org/10.1186/s12864-018-5264-y)

-------------

The Helmsman mascot was designed by Robert James Russell—view more of his work at http://www.robertjamesrussell.com/art/ and follow him on Twitter at [@robhollywood](https://twitter.com/robhollywood)!
