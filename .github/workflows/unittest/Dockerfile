##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics Sagl. All rights reserved.
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
# RUN Rscript -e "remotes::install_github('lifs-tools/rgoslin')"
# RUN Rscript -e "remotes::install_github('ncats/RAMP-DB')"

# The base image may not have playbase.ingest yet, and R CMD INSTALL does
# not resolve Remotes: install it first.
RUN Rscript -e "remotes::install_github('bigomics/playbase.ingest', dependencies = FALSE, upgrade = 'never'); if (!requireNamespace('playbase.ingest', quietly = TRUE)) stop('playbase.ingest failed to install')"

# copy bigomics/playbase code from github checkout
COPY . /playbase

# private playbase.preprocess is not in the base image. Named explicitly:
# dev/write_description.R drops Remotes packages from Imports, so
# install_deps() would skip it. PAT comes in as a build secret, not a layer.
# install_github only warns on failure, hence the explicit check.
RUN --mount=type=secret,id=gh_pat \
    GITHUB_PAT=$(cat /run/secrets/gh_pat) \
    Rscript -e "remotes::install_github('bigomics/playbase.preprocess', upgrade = 'never'); if (!requireNamespace('playbase.preprocess', quietly = TRUE)) stop('playbase.preprocess failed to install')"

RUN set -e && \
    R CMD INSTALL /playbase

# add entrypoint to execute when the docker container starts up

# echo list files in current directory

WORKDIR /playbase
RUN chmod +x /playbase/.github/workflows/unittest/entrypoint.sh

ENTRYPOINT ["/playbase/.github/workflows/unittest/entrypoint.sh"]
