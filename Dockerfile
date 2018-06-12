###############################################################################
# Modified from 
# https://github.com/jupyter/docker-stacks/blob/master/base-notebook/Dockerfile
###############################################################################
FROM jupyter/scipy-notebook:c7fb6660d096

LABEL maintainer="Jedidiah Carlson <jed.e.carlson@gmail.com>"

###############################################################################
# install apt dependencies
###############################################################################
USER root
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    fonts-dejavu \
    tzdata \
    gfortran \
    gcc && apt-get clean && \
    rm -rf /var/lib/apt/lists/*

###############################################################################
# copy gh repo to home dir
###############################################################################
USER ${NB_USER}
COPY . ${HOME}

USER root
RUN chown -R ${NB_UID} ${HOME}

USER ${NB_USER}

###############################################################################
# install dependencies from config file
###############################################################################
RUN conda env update -n root -f env.yml && \
conda clean -tipsy
