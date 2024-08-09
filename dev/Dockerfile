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

# pre-install some needed packages
RUN R -e "BiocManager::install('AnnotationDbi')"
#RUN R -e "BiocManager::install('RSQLite')"
RUN R -e "BiocManager::install('topGO')"
RUN R -e "BiocManager::install('org.Hs.eg.db')"

# preinstall dependencies
RUN Rscript dev/install_dependencies.R 

#RUN Rscript dev/install_playbase.R 'github'
RUN R -e "remotes::install_github('bigomics/playbase')"

# update all dependencies (should mostly be done..)
RUN Rscript dev/install_dependencies.R 

# reinstall some failed ones
#RUN R -e "remotes::install_github('bigomics/PCSF')"

#------------------------------------------------------------
# Clean up when done.
#------------------------------------------------------------
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN chmod ugo+rwX /tmp && chmod o+t /tmp
RUN rm -f /.Renviron /root/.Renviron /playbase/.Renviron 
RUN rm -fr /tmp/Rtmp*