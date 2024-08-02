##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

FROM ubuntu:24.04 AS os-base
##FROM ubuntu:20.10
##FROM ubuntu:20.04

MAINTAINER BigOmics "support@bigomics.ch"

# =====================================================================
# Creates base docker image for Omics Playground. Generally not for
# direct deployment but as base image for further build. Install
# necessary Debian packages and R packages 
# =====================================================================

ENV DEBIAN_FRONTEND noninteractive

RUN apt update && apt install -y \
    locales git
    
# Set the locale to UTF8
RUN sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen && locale-gen
ENV LANG=en_US.UTF-8 LANGUAGE=en_US:en LC_ALL=en_US.UTF-8 

WORKDIR /
COPY dev/install_ubuntu.sh install_ubuntu.sh
RUN sh install_ubuntu.sh

#------------------------------------------------------------
# Install R packages that are required
#------------------------------------------------------------
FROM os-base AS r-base1
WORKDIR /

# Install basic stuff and R
#RUN apt install -y r-base r-base-dev 
RUN R -e "install.packages('remotes')"

# install Playbase with dependencies using RStudio binary packages
RUN echo Downloading Playbase
COPY dev/rspm.R rspm.R
RUN R -e "source('rspm.R'); remotes::install_github('bigomics/playbase',dependencies=TRUE)"

## add any missing packages manually------------------------
FROM r-base1 AS r-base2
#RUN R -e "BiocManager::install(c('rhdf5', 'monocle', 'mygene'))"
#RUN R -e "install.packages(c('RSpectra'))"

## playdata: perhaps not needed because playbase already states
## playdata as dependency in DESCRIPTION file.
#RUN git clone -b main --single-branch --depth 1 \
#    https://github.com/bigomics/playdata.git \
#    && R CMD INSTALL playdata && rm -fr playdata

#------------------------------------------------------------
# Clean up when done.
#------------------------------------------------------------
#RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
#RUN chmod ugo+rwX /tmp && chmod o+t /tmp

