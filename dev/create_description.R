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
cat("RENV:: building dependencies...\n")
renv.out <- renv::dependencies(path = "R", root = getwd(), errors = "ignored")
pkg.used <- unique(renv.out$Package)
pkg.used <- setdiff(pkg.used, "playbase")
cat("RENV:: done!\n")

## Define remote locations or versions
pkg.remotes <- c(
  "playdata" = "bigomics/playdata",
  "NNLM" = "linxihui/NNLM",
  "KEGG.db" = "ncullen93/KEGG.db",
  "iheatmapr" = "ropensci/iheatmapr",
  "SuperCell" = "GfellerLab/SuperCell",
  "Azimuth" = "satijalab/azimuth",
  "org.Pf.plasmo.db" = "url::https://www.bioconductor.org/packages/3.14/data/annotation/src/contrib/org.Pf.plasmo.db_3.14.0.tar.gz",
  "infercnv" = "broadinstitute/infercnv@infercnv-v1.3.3",
  "PCSF" = "bigomics/PCSF",
  "EPIC" = "GfellerLab/EPIC"
)

remotes.url <- c(
  "PCSF" = "url::https://github.com/bigomics/PCSF/archive/HEAD.zip",
  "playdata" = "url::https://github.com/bigomics/playdata/archive/HEAD.zip",
  "infercnv" = "url::https://github.com/broadinstitute/infercnv/archive/HEAD.zip",
  "EPIC" = "url::https://github.com/GfellerLab/EPIC/archive/HEAD.zip",
  "SuperCell" = "url::https://github.com/GfellerLab/SuperCell/archive/HEAD.zip",
  "NNLM" = "url::https://github.com/linxihui/NNLM/archive/HEAD.zip",
  "iheatmapr" = "url::https://github.com/ropensci/iheatmapr/archive/HEAD.zip",
  "Azimuth" = "url::https://github.com/satijalab/azimuth/archive/HEAD.zip",
  "org.Pf.plasmo.db" = "url::https://bioconductor.org/packages/3.14/data/annotation/src/contrib/org.Pf.plasmo.db_3.14.0.tar.gz",
  "KEGG.db" = "url::https://bioconductor.org/packages/3.11/data/annotation/src/contrib/KEGG.db_3.2.4.tar.gz"
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


