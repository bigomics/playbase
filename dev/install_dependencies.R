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
use.remotes <- ifelse(dir.exists('R'), FALSE, TRUE)
install_dependencies(use.remotes) 
