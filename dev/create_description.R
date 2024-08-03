##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

# This file is supposed to run from the root Playground folder
if (basename(getwd()) != "playbase") {
    stop("This file is supposed to run from the root Playbase folder")
}

if(!require("renv")) install.packages("renv")
if(!require("remotes")) install.packages("remotes")
if(!require("devtools")) install.packages("devtools")

## ---------------------------------------------------------------------
## Automatically scan all used packages and install
## ---------------------------------------------------------------------
## We use renv to detect dependencies. Renv is looking for library and
## require statements in the r/R source files.
renv.out <- renv::dependencies(path = "R", root = getwd(), errors = "ignored")
pkg.used <- unique(renv.out$Package)
pkg.used <- setdiff(pkg.used, "playbase")

## Define remote locations or versions
github_url <- function(repo) paste0("github::",repo)
github_url <- function(repo) {
  if(grepl("@",repo)) {
    branch <- sub(".*@","",repo)
    repo <- sub("@.*","",repo)    
    paste0("url::https://github.com/",repo,"/archive/refs/heads/",branch,".zip")
  } else {
    paste0("url::https://github.com/",repo,"/archive/HEAD.zip")
  }
}

add_github <- function(repo) {
  pkg.name <- gsub(".*[/]|@.*","",repo)
  remotes.url[pkg.name] <<- github_url(repo)
}
remotes.url <- c(
  "KEGG.db" = "url::https://bioconductor.org/packages/3.11/data/annotation/src/contrib/KEGG.db_3.2.4.tar.gz",
  "org.Pf.plasmo.db" = "url::https://bioconductor.org/packages/3.14/data/annotation/src/contrib/org.Pf.plasmo.db_3.14.0.tar.gz",
  "Azimuth" = "url::https://github.com/satijalab/azimuth/archive/HEAD.zip"
)
add_github("bigomics/PCSF")
add_github("bigomics/playdata")
add_github("bigomics/playbase")
add_github("bigomics/bigdash")
add_github("bigomics/bigLoaders")
add_github("bigomics/fgsea")
add_github("bigomics/wizardR")
add_github("GfellerLab/EPIC")
add_github("broadinstitute/infercnv")
add_github("GfellerLab/SuperCell")
add_github("linxihui/NNLM")
add_github("Coolgenome/iTALK")
add_github("wt2015-github/FastGGM")
add_github("satijalab/azimuth")
add_github("JohnCoene/waiter")
add_github("JohnCoene/firebase@omics")
add_github("JohnCoene/bsutils")
add_github("ropensci/iheatmapr")
##add_github("rstudio/bslib@v0.6.1")
add_github("rstudio/htmltools")
add_github("bigomics/biomaRt")
add_github("Bioconductor/BiocFileCache")
add_github("cysouw/qlcMatrix")
add_github("cole-trapnell-lab/leidenbase")
add_github('cole-trapnell-lab/monocle3')
add_github('bartongroup/Proteus')
add_github('cran/riverplot')
add_github('Ironholds/rgeolocate')


pkg.remotes <- remotes.url[names(remotes.url) %in% pkg.used]
pkg.imports <- setdiff(pkg.used, names(pkg.remotes))

if(file.exists("DESCRIPTION")) {
  hash <- substring(tempfile(),nchar(tempfile())-3, nchar(tempfile()))
  file.copy("DESCRIPTION",paste0("DESCRIPTION.save.",hash))
  message("WARNING: overwriting existing DESCRIPTION file.")
}

desc.header <- readLines("dev/DESCRIPTION.header")
desc.file <- "DESCRIPTION"

write(desc.header, file=desc.file)

write("Imports:", file=desc.file, append=TRUE)
write(paste0("    ",sort(pkg.imports),","), file=desc.file, append=TRUE)

write("Remotes:", file=desc.file, append=TRUE)
write(paste0("    ",sort(pkg.remotes),","), file=desc.file, append=TRUE)


