##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

## Start from base image, update git code and add data
## folder. Create docker ready to be deployed.

#------------------------------------------------------------
# Start from lastest base image
#------------------------------------------------------------

ARG DOCKER_TAG=bigomics/omicsplayground:latest

FROM ${DOCKER_TAG}  

#------------------------------------------------------------
# Set env variables 
#------------------------------------------------------------
WORKDIR /
ARG BRANCH=master
ARG BRANCH_PLAYBASE=main
ARG GITHUB_SHA=master

ENV LANG=en_US.UTF-8 LANGUAGE=en_US:en LC_ALL=en_US.UTF-8 
RUN sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen && locale-gen

#------------------------------------------------------------
# Install de-novo bigomics stuff (as R packages)
#------------------------------------------------------------

RUN R -e "remotes::install_github('bigomics/PCSF',dependencies=FALSE)"
RUN R -e "remotes::install_github('bigomics/fgsea')"


# copy bigomics/playbase code from github checkout
COPY . /playbase

RUN R CMD INSTALL /playbase 

# add entrypoint to execute when the docker container starts up

ENTRYPOINT ["/playbase/.github/workflows/build-test-push/entrypoint.sh"]