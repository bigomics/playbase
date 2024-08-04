## ---------------------------------------------------------------------
## Install Playbase and dependencies from GitHub as packages using R
## package manager
## ---------------------------------------------------------------------

if(basename(getwd())=="dev") setwd('..')

# This file is supposed to run from the root Playground folder
if (basename(getwd()) != "playbase") {
    stop("This file is supposed to run from the root Playbase folder")
}

source("dev/rspm.R")
require <- function(pkg) (pkg %in% installed.packages()[,'Package'])

if(!require("remotes")) install.packages('remotes')
if(!require("BiocManager")) {
  remotes::install_version('BiocManager', version='1.30.23')
  if(!BiocManager::version()=="3.18") BiocManager::install(version='3.18')
}

## add problematic packages first
if(!require("msa")) BiocManager::install('msa')
if(!require("orthogene")) BiocManager::install('orthogene')
if(!require("MSnbase")) BiocManager::install('MSnbase')

# install Playbase with dependencies
remotes::install_github('bigomics/playbase',dependencies=TRUE)

## add any missed packages manually
if(!require("topGO")) BiocManager::install(c('topGO'))
if(!require("PCSF"))  remotes::install_github('bigomics/PCSF')

# remove big non-used packages
try(remove.packages('BSgenome.Hsapiens.UCSC.hg38'))
try(remove.packages('EnsDb.Hsapiens.v86'))
