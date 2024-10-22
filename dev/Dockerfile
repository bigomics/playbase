##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

FROM playbase-rbase AS playbase

#------------------------------------------------------------
# Install playbase 
#------------------------------------------------------------

RUN echo ********* INSTALLING PLAYBASE ***************
WORKDIR /playbase

COPY DESCRIPTION ./DESCRIPTION
COPY dev/*.R Makevars* dev/
COPY dev/Rprofil[e] /root/.Rprofile
# if a github PAT is available put it in Renviron (but remove after!!!)
COPY dev/Renviro[n] /root/.Renviron


RUN Rscript dev/install_playbase.R 'github'

# update dependencies (should mostly be done..)
RUN Rscript dev/install_dependencies.R 

# reinstall some failed ones
#RUN R -e "remotes::install_github('bigomics/PCSF')"
RUN R -e "BiocManager::install('TileDBArray')"

#------------------------------------------------------------
# Clean up when done.
#------------------------------------------------------------
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN chmod ugo+rwX /tmp && chmod o+t /tmp
RUN rm -f /.Renviron /root/.Renviron /playbase/.Renviron 
RUN rm -fr /tmp/Rtmp*