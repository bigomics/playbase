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

source("dev/rspm.R")
source("dev/functions.R")

require <- function(pkg) (pkg %in% installed.packages()[,'Package'])

if(!require("remotes")) install.packages('remotes')
if(!require("BiocManager")) {
  remotes::install_version('BiocManager', version='1.30.23')
  if(!BiocManager::version()=="3.18") BiocManager::install(version='3.18')
}

P <- installed.packages()
imports <- gsub("[\n ]","",strsplit(P["playbase","Imports"], split=",")[[1]])
imports <- trimws(sub("@.*","",imports))
missing <- setdiff( imports, rownames(P) )

if(length(missing)) {
  message("Missing packages (",length(missing),"):",paste(missing,collapse=" "))
  for(p in missing) {
    if(!require(p)) BiocManager::install(p, ask=FALSE, dependencies=TRUE)
  }
} else {
  message("No missing packages")
}

## update all packages
message("Updating packages...")
#update.packages(ask=FALSE, lib.loc="/usr/local/lib/R/site-library")
BiocManager::install(ask=FALSE, lib.loc="/usr/local/lib/R/site-library")
