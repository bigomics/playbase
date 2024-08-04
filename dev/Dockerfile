##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

FROM playbase-os AS playbase

#------------------------------------------------------------
# Install R packages that are required
#------------------------------------------------------------
# Install basic stuff and R
RUN apt install -y r-base r-base-dev 

RUN mkdir -p /playbase/dev
COPY dev/ /playbase/dev/

WORKDIR /playbase
RUN cp dev/rspm.R /root/.Rprofile

# if a github PAT is available put it here
COPY github_pat /root/.Renviron

RUN R -e "install.packages('remotes')"
RUN R -e "remotes::install_version('BiocManager', version='1.30.23')"
RUN R -e "BiocManager::install(version='3.18')"

RUN Rscript dev/install_conda.R && rm -fr ~/.local/share/r-miniconda/pkgs*
RUN Rscript dev/install_playbase.R && rm -fr /tmp/Rtmp*

#------------------------------------------------------------
# Clean up when done.
#------------------------------------------------------------
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN chmod ugo+rwX /tmp && chmod o+t /tmp
RUN rm -f /.Renviron /root/.Renviron 
