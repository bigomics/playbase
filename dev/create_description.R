##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

# This file is supposed to run from the root Playground folder
if (basename(getwd()) != "playbase") {
    stop("This file is supposed to run from the root Playbase folder")
}

## ---------------------------------------------------------------------
## Automatically scan all used packages and install
## ---------------------------------------------------------------------
## We use renv to detect dependencies. Renv is looking for library and
## require statements in the r/R source files.
## install.packages("renv")
renv.out <- renv::dependencies(path = "R", root = getwd(), errors = "ignored")
pkg.used <- unique(renv.out$Package)
pkg.used <- setdiff(pkg.used, "playbase")

## Define remote locations or versions
github_url <- function(repo) paste0("github::",repo)
github_url <- function(repo) {
  paste0("url::https://github.com/",repo,"/archive/HEAD.zip")
}

remotes.url <- c(
  "PCSF"      = github_url("bigomics/PCSF"),
  "playdata"  = github_url("bigomics/playdata"),
  "infercnv"  = github_url("broadinstitute/infercnv"),
  "EPIC"      = github_url("GfellerLab/EPIC"),
  "SuperCell" = github_url("GfellerLab/SuperCell"),
  "NNLM"      = github_url("linxihui/NNLM"),
  "iheatmapr" = github_url("ropensci/iheatmapr"),
  "Azimuth"   = github_url("satijalab/azimuth"),
  "org.Pf.plasmo.db" = "url::https://bioconductor.org/packages/3.14/data/annotation/src/contrib/org.Pf.plasmo.db_3.14.0.tar.gz",
  "KEGG.db"   = "url::https://bioconductor.org/packages/3.11/data/annotation/src/contrib/KEGG.db_3.2.4.tar.gz"
)

#pkg.remotes <- pkg.remotes[names(pkg.remotes) %in% pkg.used]
pkg.remotes <- remotes.url[names(remotes.url) %in% pkg.used]
pkg.imports <- setdiff(pkg.used, names(pkg.remotes))

if(file.exists("DESCRIPTION")) {
  stop("ERROR: existing DESCRIPTION file. Please remove first.")
}
desc.header <- readLines("dev/DESCRIPTION.header")
desc.file <- "DESCRIPTION"

write(desc.header, file=desc.file)

write("Imports:", file=desc.file, append=TRUE)
write(paste0("    ",sort(pkg.imports),","), file=desc.file, append=TRUE)

write("Remotes:", file=desc.file, append=TRUE)
write(paste0("    ",sort(pkg.remotes),","), file=desc.file, append=TRUE)


