##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

require <- function(pkg) (pkg %in% installed.packages()[,'Package'])

scan_packages <- function(path='R') {
  
  ## ---------------------------------------------------------------------
  ## Automatically scan all used packages and install
  ## ---------------------------------------------------------------------
  ## We use renv to detect dependencies. Renv is looking for library and
  ## require statements in the r/R source files.
  renv.out <- renv::dependencies(path = path, root = getwd(), errors = "ignored")
  pkg.used <- sort(unique(renv.out$Package))

  ## scan roxygen2 @import and @importFrom directives (renv misses these)
  r.files <- list.files(path, pattern = "\\.[rR]$", full.names = TRUE, recursive = TRUE)
  roxygen.pkgs <- c()
  for(f in r.files) {
    lines <- readLines(f, warn = FALSE)
    importfrom.lines <- grep("^#'\\s*@importFrom\\s", lines, value = TRUE)
    pkgs <- gsub("^#'\\s*@importFrom\\s+(\\S+).*", "\\1", importfrom.lines)
    roxygen.pkgs <- c(roxygen.pkgs, pkgs)
    import.lines <- grep("^#'\\s*@import\\s", lines, value = TRUE)
    import.lines <- grep("@importFrom", import.lines, value = TRUE, invert = TRUE)
    pkgs <- gsub("^#'\\s*@import\\s+", "", import.lines)
    pkgs <- unlist(strsplit(pkgs, "\\s+"))
    roxygen.pkgs <- c(roxygen.pkgs, pkgs)
  }
  pkg.used <- unique(c(pkg.used, roxygen.pkgs))

  ## manually add some missing packages
  pkg.used <- c(pkg.used, c("sf"))
  pkg.used <- setdiff(pkg.used, "playbase")

  ## Define remote locations or versions
  github_url <- function(repo) paste0("github::",repo)
  github_url <- function(repo) {
    if(grepl("@",repo)) {
      branch <- sub(".*@","",repo)
      repo <- sub("@.*","",repo)    
      paste0("url::https://github.com/",repo,"/archive/refs/tags/",branch,".zip")
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
    "Azimuth" = "url::https://github.com/satijalab/azimuth/archive/refs/heads/master.zip"
  )

  force.remotes.url <- c(
    "rjson" = "url::https://cran.r-project.org/src/contrib/Archive/rjson/rjson_0.2.21.tar.gz",
    "rms" = "url::https://cran.r-project.org/src/contrib/Archive/rms/rms_6.8-0.tar.gz"    
  )
  
  ## commented out entries are now in standard CRAN/cBio repo
  add_github("bigomics/PCSF")
  add_github("bigomics/playdata")
  add_github("bigomics/playbase")
  add_github("bigomics/bigdash")
  add_github("bigomics/bigLoaders")
  add_github("bigomics/wizardR")
  add_github("GfellerLab/EPIC")
  add_github("broadinstitute/infercnv@infercnv-v1.3.3")
  add_github("GfellerLab/SuperCell")
  add_github("linxihui/NNLM")
  add_github("Coolgenome/iTALK")
  add_github("wt2015-github/FastGGM")
  add_github("JohnCoene/firebase@omics")
  add_github("JohnCoene/bsutils")
  add_github("Bioconductor/BiocFileCache")
  add_github('cole-trapnell-lab/monocle3')
  add_github('bartongroup/Proteus')
  add_github('cran/riverplot')
  add_github('Ironholds/rgeolocate')

  pkg.remotes <- remotes.url[names(remotes.url) %in% pkg.used]
  pkg.imports <- setdiff(pkg.used, names(pkg.remotes))

  pkg.installed <- installed.packages()[,'Package']
  pkg.missing <- setdiff( c(pkg.imports,names(pkg.remotes)), pkg.installed)
  missing.imports <- setdiff(pkg.imports, pkg.installed)
  missing.remotes <- pkg.remotes[!(names(pkg.remotes) %in% pkg.installed)]
  missing.remotes <- c(missing.remotes, force.remotes.url)
  
  list(
    used = pkg.used,
    installed = pkg.installed,
    missing = pkg.missing,
    imports = pkg.imports,
    remotes = pkg.remotes,    
    missing.imports = missing.imports,
    missing.remotes = missing.remotes
  )

}

##org="Hs";version="3.21.0"
install_newest_orgdb <- function(org="Hs", version="3.19.1") {
  db <- sub("Hs",org,'org.Hs.eg.db')
  if(require(db) && packageVersion(db) >= version) return(NULL)  
  CMD <- paste0("cd /tmp && wget https://bioconductor.org/packages/release/data/annotation/src/contrib/org.Hs.eg.db_",version,".tar.gz && tar xvfz org.Hs.eg.db_",version,".tar.gz && cd org.Hs.eg.db && sed -i 's/, AnnotationDbi.*//' DESCRIPTION && R CMD INSTALL .")
  if(org!="Hs") CMD <- gsub("Hs",org,CMD)
  system(CMD)
}

install_dependencies <- function(use.remotes=FALSE) {
  
  require <- function(pkg) (pkg %in% installed.packages()[,'Package'])
  remove.pkg <- function(p) if(require(p, character.only=TRUE)) try(remove.packages(p))

  if(!require("remotes")) install.packages('remotes')
  if(!require("BiocManager")) {
    remotes::install_version('BiocManager', version='1.30.23')
    if(!BiocManager::version()=="3.18") BiocManager::install(version='3.18')
  }

  if(use.remotes) {
    # install dependencies using remotes
    remotes::install_deps('.', dependencies = c("Imports","Remotes"))
  } else {
    ## Here we handle missing dependencies ourselves. Better control of
    ## skipping packages that are already installed.
    pkg <- scan_packages('R')
    if( length(pkg$missing.imports) || length(pkg$missing.remotes) ) {
      message("Installing missing dependencies...")
      for(p in pkg$missing.imports) {
        if(!require(p)) BiocManager::install(p, ask=FALSE, dependencies=TRUE)
      }
      for(p in names(pkg$missing.remotes)) {
        url <- sub("url::","",pkg$missing.remotes[p])
        if(!require(p)) remotes::install_url(url, ask=FALSE, dependencies=TRUE)
      }
    } else {
      message("All dependencies aleady installed")
    }
  }

  ## overwrite
  ## since 31.5.2025 (IK)
  install_newest_orgdb("Hs","3.21.0")  
  install_newest_orgdb("Mm","3.21.0")
  install_newest_orgdb("Rn","3.21.0")   
}

install_silent <- function(pkg.list, linkto=NULL, force=FALSE) {
  ## suppress ultra-verbose packages that issue loads of warnings
  ## during compilation.
  if(!is.null(linkto)) {
    P <- available.packages()
    pkg.noisy <- c("PCSF","SeuratObject","RSpectra")
    for( k in linkto) {
      pkg1 <- rownames(P)[grep(k,as.character(P[,"LinkingTo"]))]
      pkg.noisy <- c(pkg.noisy, k, pkg1)
    }
    pkg <- intersect(pkg.list, pkg.noisy)
  } else {
    pkg <- pkg.list
  }
  
  if(length(pkg) > 0) {
    message("> Pre-installing ultra-verbose packages: ", paste(pkg,collapse=" "))
    if(!dir.exists("~/.R")) dir.create("~/.R")
    if(file.exists("~/.R/Makevars")) file.copy("~/.R/Makevars","~/.R/Makevars.save")
    file.copy("dev/Makevars.no-error","~/.R/Makevars")
    for(p in pkg) {
      if(grepl("url::",p)) {
        p <- sub("url::","",p)        
        remotes::install_url(p, ask=FALSE, force=force)
      } else if(grepl("github::",p)) {
        p <- sub("github::","",p)
        remotes::install_github(p, ask=FALSE, force=force)
      } else {
        BiocManager::install(p, ask=FALSE, force=force)
      }
    }
    file.remove("~/.R/Makevars")
    if(file.exists("~/.R/Makevars.save")) file.copy("~/.R/Makevars.save","~/.R/Makevars")
    pkg.list <- setdiff(pkg.list, pkg)
  } else {
    print("> No ultra-verbose packages!")
  }
  pkg.list
}

scan_description <- function(path=NULL) {

  if(is.null(path)) {
    search_paths = c('.', '/usr/lib/R/library/playbase',
                     '/usr/local/lib/R/site-library/playbase')
    sel <- which(file.exists(file.path(search_paths,"DESCRIPTION")))
    if(length(sel) > 0) {
      path <- search_paths[sel[1]]
    }
  }
  if(is.null(path)) stop("could not find playbase DESCRIPTION file")
  desc <- readLines(file.path(path,'DESCRIPTION'))
  sel1 <- (grep("Imports:",desc)) : (grep("Remotes:",desc)-1)
  imports <- desc[sel1]
  imports <- setdiff(unlist(strsplit( imports, split='[ ]')), c("","Imports:"))
  imports <- trimws(gsub("[,]","",imports))
  imports
  
  sel2 <- (grep("Remotes:",desc)+1) : length(desc)
  remotes <- trimws(gsub("[,]","",desc[sel2]))
  remotes <- grep("^url::", remotes, value=TRUE)
  remotes <- sub("^url::","",remotes)
  remotes.names <- gsub(".*contrib/|.*github.com/","",remotes)
  remotes.names <- gsub("_[0-9].*|/archive.*","",remotes.names)
  remotes.names <- gsub(".*/","",remotes.names)  
  names(remotes) <- remotes.names
  
  list(
    imports = imports,
    remotes = remotes
  )  
}
