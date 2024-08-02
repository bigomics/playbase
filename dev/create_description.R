##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


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
  "xbioc" = "renozao/xbioc",
  "SuperCell" = "GfellerLab/SuperCell",
  "Azimuth" = "satijalab/azimuth",
  "org.Pf.plasmo.db" = "url::https://www.bioconductor.org/packages//2.12/data/annotation/src/contrib/org.Pf.plasmo.db_2.9.0.tar.gz",
  "infercnv" = "broadinstitute/infercnv@infercnv-v1.3.3",
  "PCSF" = "bigomics/PCSF",
  "EPIC" = "GfellerLab/EPIC"
)

pkg.remotes <- pkg.remotes[names(pkg.remotes) %in% pkg.used]
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


