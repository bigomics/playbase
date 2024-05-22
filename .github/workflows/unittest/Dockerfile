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

ENV LANG=en_US.UTF-8 LANGUAGE=en_US:en LC_ALL=en_US.UTF-8 
RUN sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen && locale-gen
SHELL ["/bin/bash", "-o", "pipefail", "-c"]

#------------------------------------------------------------
# Install de-novo bigomics stuff (as R packages)
#------------------------------------------------------------

# copy bigomics/playbase code from github checkout
COPY . /playbase

# This can be removed once new base is built
RUN R -e "remotes::install_version('matrixStats',version='1.1.0',force=TRUE)"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/Matrix.utils/Matrix.utils_0.9.8.tar.gz', type = 'source', repos = NULL)"
RUN R -e "remotes::install_github('cran/Matrix', ref = '23d9d53')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/Matrix.utils/Matrix.utils_0.9.8.tar.gz', repos = NULL, type = 'source')"
RUN R -e "BiocManager::install('orthogene', dependencies = NA, force = FALSE, ask = FALSE, update = FALSE)"
RUN R -e "install.packages('plsRcox', dependencies = FALSE)"
RUN R -e "remotes::install_version('SeuratObject', '4.1.4', repos = c('https://satijalab.r-universe.dev', getOption('repos')))"
RUN R -e "remotes::install_version('Seurat', '4.4.0', repos = c('https://satijalab.r-universe.dev', getOption('repos')))"
RUN R -e "BiocManager::install('org.Rn.eg.db', dependencies = NA, force = FALSE, ask = FALSE, update = FALSE)"
RUN R -e "BiocManager::install('AnnotationHub',dependencies = NA, force = FALSE, ask = FALSE, update = FALSE)"
RUN R -e "devtools::install_version('dbplyr', version = '2.3.4')"

RUN set -e && \
    R CMD INSTALL /playbase

# add entrypoint to execute when the docker container starts up

# echo list files in current directory

WORKDIR /playbase
RUN chmod +x /playbase/.github/workflows/unittest/entrypoint.sh

ENTRYPOINT ["/playbase/.github/workflows/unittest/entrypoint.sh"]