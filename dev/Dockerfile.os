##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2020 BigOmics Analytics Sagl. All rights reserved.
##

FROM ubuntu:24.04 AS playbase-os

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

# disable c++ version warning (msa package)
RUN mv /usr/include/c++/13/bits/c++0x_warning.h \
    /usr/include/c++/13/bits/c++0x_warning.h.DISABLED && \
    touch /usr/include/c++/13/bits/c++0x_warning.h 

#------------------------------------------------------------
# Clean up when done.
#------------------------------------------------------------
RUN apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN chmod ugo+rwX /tmp && chmod o+t /tmp
