## ---------------------------------------------------------------------
## Install Playbase and dependencies from GitHub as packages using R
## package manager
## ---------------------------------------------------------------------
#!/usr/bin/env Rscript

# This file is supposed to run from the root Playground folder
if(basename(getwd())=="dev") setwd('..')
if (basename(getwd()) != "playbase") {
    stop("This file is supposed to run from the root Playbase folder")
}

# test if there is at least one argument: if not, return an error
args = commandArgs(trailingOnly=TRUE)

source = "github"
if (length(args)==0) {
  source <- "github"
} else if (length(args)>=1) {
  source <- as.character(args[1])
}
if(!(source %in% c("local","github","rcmd"))) {
  stop("source argument must be 'local', 'github' or 'rcmd'")
}
message("installing playbase from ", source)

source("dev/rspm.R")
source("dev/functions.R")

require <- function(pkg) (pkg %in% installed.packages()[,'Package'])

if(!require("remotes")) install.packages('remotes')
if(!require("BiocManager")) {
  remotes::install_version('BiocManager', version='1.30.23')
  if(!BiocManager::version()=="3.18") BiocManager::install(version='3.18')
}

## add problematic packages first
if(any(c("RccpEigen","ggforce") %in% missing.imports)) {
  if(!dir.exists("~/.R")) dir.create("~/.R")
  if(!file.exists("~/.R/Makevars")) file.copy("dev/Makevars","~/.R/Makevars.save")
  file.copy("dev/Makevars.no-error","~/.R/Makevars")
  BiocManager::install(c("RccpEigen","ggforce"),ask=FALSE)
  file.remove("~/.R/Makevars")
  if(!file.exists("~/.R/Makevars.save")) file.copy("dev/Makevars.save","~/.R/Makevars")  
}
if(!require("msa")) BiocManager::install('msa')
if(!require("orthogene")) BiocManager::install('orthogene')
if(!require("MSnbase")) BiocManager::install('MSnbase')


## Install dependencies first
install_dependencies( use.remotes=TRUE )

## Install playbase without dependencies
if(source == 'github') {
  remotes::install_github('bigomics/playbase',dependencies=FALSE)
} else if(source == 'local') {
  remotes::install_local(paste='.', dependencies=FALSE)
} else if(source == 'rcmd') {
  system("R CMD INSTALL .")
}

## add any missed packages manually
if(!require("topGO")) BiocManager::install(c('topGO'))
if(!require("org.Hs.eg.db")) BiocManager::install(c('org.Hs.eg.db'))
if(!require("PCSF"))  remotes::install_github('bigomics/PCSF')

# remove big non-used packages
remove.pkg <- function(p) if(require(p)) try(remove.packages(p))
remove.pkg('BSgenome.Hsapiens.UCSC.hg38')
remove.pkg('EnsDb.Hsapiens.v86')
