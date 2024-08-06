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
  source <- ifelse(dir.exists('R'), "local", "github")
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

missing.imports <- NULL

if(dir.exists('R')) {
  pkg <- scan_packages(path='R')
  missing.imports <- pkg$imports.missing
}

if(is.null(missing.imports) && file.exists('DESCRIPTION')) {
  desc <- scan_description('.') 
  missing.imports <- desc$imports
}

## suppress ultra-verbose packages that link to RccpEigen
if(!is.null(missing.imports) && length(missing.imports)) {
  missing.imports <- install_silent(
    missing.imports,
    linkto = c("RcppEigen","BH")
  )
}

## add problematic packages first
if(!require("msa")) BiocManager::install('msa')
if(!require("orthogene")) BiocManager::install('orthogene')
if(!require("MSnbase")) BiocManager::install('MSnbase')

## Install dependencies first
message('>>> installing dependencies...')
if(source == "github") {
  install_dependencies( use.remotes = TRUE )
} else {
  install_dependencies( use.remotes = FALSE )
}

## Install playbase 
if(source == 'github') {
  message('>>> installing playbase from github')
  remotes::install_github('bigomics/playbase',dependencies=FALSE)
} else if(source == 'local') {
  message('>>> installing playbase from local R folder (install_local)')
  remotes::install_local(paste='.', dependencies=FALSE)
} else if(source == 'rcmd') {
  message('>>> installing playbase from local R folder (R CMD INSTALL)')  
  system("R CMD INSTALL .")
}

## try again any missed packages manually
desc <- scan_description("/usr/local/lib/R/site-library/playbase")
installed.pkgs <- installed.packages()[,"Package"]
missing.imports <- setdiff( desc$imports, installed.pkgs )
missing.remotes <- desc$remotes[setdiff( names(desc$remotes), installed.pkgs )]

if(length(missing.imports)) {
  message(">>> Retrying failed imports: ", paste(missing.imports,collapse=" "))
  BiocManager::install( missing.imports )
}
if(length(missing.remotes)) {
  message(">>> Retrying failed remotes: ", paste(names(missing.remotes),collapse=" "))  
  for(r in missing.remotes) {
    remotes::install_url(r)
  }
}

## final few that often fail
if(!require("topGO")) BiocManager::install(c('topGO'))
if(!require("org.Hs.eg.db")) BiocManager::install(c('org.Hs.eg.db'))
if(!require("PCSF"))  install_silent('github::bigomics/PCSF')

# remove big non-used(?) packages
remove.pkg <- function(p) if(require(p)) try(remove.packages(p))
remove.pkg('BSgenome.Hsapiens.UCSC.hg38')  ## 787Mb
remove.pkg('EnsDb.Hsapiens.v86')           ## 333Mb
