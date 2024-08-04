##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

FROM playbase-os

#------------------------------------------------------------
# Install R packages that are required
#------------------------------------------------------------
WORKDIR /

# Install basic stuff and R
RUN apt install -y r-base r-base-dev 
COPY dev/rspm.R /root/.Rprofile

# if a github PAT is available put it here
COPY github_pat /root/.Renviron

RUN R -e "install.packages('remotes')"
RUN R -e "remotes::install_version('BiocManager', version='1.30.23')"
RUN R -e "BiocManager::install(version='3.18')"

COPY dev/install_conda.R /tmp/install_conda.R
RUN Rscript /tmp/install_conda.R

## add problematics packages manually
RUN R -e "BiocManager::install('msa')"
RUN R -e "BiocManager::install('orthogene')"
RUN R -e "BiocManager::install('MSnbase')" 

# install Playbase with dependencies using RStudio binary packages
RUN echo Downloading Playbase
RUN R -e "remotes::install_github('bigomics/playbase',dependencies=TRUE)"

## add any missed packages manually
RUN R -e "BiocManager::install(c('topGO'))"
RUN R -e "remotes::install_github('bigomics/PCSF')"

#------------------------------------------------------------
# Clean up when done.
#------------------------------------------------------------
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN chmod ugo+rwX /tmp && chmod o+t /tmp
RUN rm -f /.Renviron /root/.Renviron 
