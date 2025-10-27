##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

## Start from base image, update git code and add data
## folder. Create docker ready to be deployed.

#------------------------------------------------------------
# Start from lastest base image
#------------------------------------------------------------

ARG DOCKER_TAG=bigomics/omicsplayground:devel

FROM ${DOCKER_TAG}   

#------------------------------------------------------------
# Set env variables 
#------------------------------------------------------------
WORKDIR /

ENV LANG=en_US.UTF-8 LANGUAGE=en_US:en LC_ALL=en_US.UTF-8 
RUN sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen && locale-gen
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

#------------------------------------------------------------
# Install de-novo bigomics stuff (as R packages)
#------------------------------------------------------------
RUN Rscript -e "remotes::install_github('lifs-tools/rgoslin')"
RUN Rscript -e "remotes::install_github('ncats/RAMP-DB')"

# copy bigomics/playbase code from github checkout
COPY . /playbase

RUN set -e && \
    R CMD INSTALL /playbase

# add entrypoint to execute when the docker container starts up

# echo list files in current directory

WORKDIR /playbase
RUN chmod +x /playbase/.github/workflows/unittest/entrypoint.sh

ENTRYPOINT ["/playbase/.github/workflows/unittest/entrypoint.sh"]
