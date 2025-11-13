##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


## ================================================================================
## PGXINFO methods (preferred API)
## ================================================================================

#' Add new pgx object to PGX-table
#'
#' @param pgxinfo The existing pgxinfo data frame containing dataset metadata
#' @param pgx The pgx object containing updated metadata to add
#' @param remove.old Logical indicating whether to remove existing entries for the same dataset. Default is TRUE.
#'
#' @return Updated pgxinfo data frame with additional rows from pgx
#'
#' @description
#' Updates the pgxinfo dataset metadata table with information from a new pgx object.
#'
#' @details
#' This function takes an existing pgxinfo data frame and a pgx object as input.
#' It extracts the metadata stored in pgx$info and appends it as a new row to the pgxinfo table.
#'
#' If remove.old is TRUE, it will first remove any existing rows for the same dataset
#' before appending the new row. This avoids duplicating information for the same dataset.
#'
#' The updated pgxinfo data frame containing all dataset metadata is returned.
#'
#' @export
pgxinfo.add <- function(pgxinfo, pgx, remove.old = TRUE) {
  cond <- grep("title|source|group|batch|sample|patient|donor|repl|clone|cluster|lib.size|^[.]",
    colnames(pgx$samples),
    invert = TRUE, value = TRUE
  )
  cond <- grep("title|source|batch|sample|patient|donor|repl|clone|cluster|lib.size|^[.]",
    colnames(pgx$samples),
    invert = TRUE, value = TRUE
  )

  organism <- pgx.getOrganism(pgx)

  this.date <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  date <- ifelse(is.null(pgx$date), this.date, as.character(pgx$date))
  dataset.name <- pgx$name

  creator <- ifelse(!is.null(pgx$creator), pgx$creator, "")

  this.info <- c(
    dataset = dataset.name,
    ## author = "", ## add author? maintainer? owner??
    creator = creator,
    datatype = ifelse(is.null(pgx$datatype), "", pgx$datatype),
    description = ifelse(is.null(pgx$description), "", pgx$description),
    organism = organism,
    nsamples = nrow(pgx$samples),
    nfeatures = nrow(pgx$X),
    nsets = nrow(pgx$gsetX),
    conditions = paste(cond, collapse = " "),
    date = as.character(date),
    path = NULL
  )

  ## force to be character...
  if (!is.null(pgxinfo) && NCOL(pgxinfo) > 0 && nrow(pgxinfo) > 0) {
    if ("date" %in% colnames(pgxinfo)) {
      pgxinfo$date <- as.character(pgxinfo$date)
    }
    if (!"creator" %in% colnames(pgxinfo)) {
      pgxinfo$creator <- ""
    }
    which.factor <- which(sapply(pgxinfo, is.factor))
    which.factor
    for (i in which.factor) {
      pgxinfo[, i] <- as.character(pgxinfo[, i])
    }

    ## remove existing entries??
    if (remove.old && nrow(pgxinfo) > 0) {
      d1 <- sub("[.]pgx$", "", pgxinfo$dataset)
      d2 <- sub("[.]pgx$", "", this.info["dataset"])
      if (!is.null(d2) && !is.na(d2) && d2 %in% d1 && d2 != "") {
        sel <- which(d1 != d2)
        pgxinfo <- pgxinfo[sel, , drop = FALSE]
      }
    }

    ## merge with same columns
    info.cols <- colnames(pgxinfo)
    info.cols <- unique(c(info.cols, names(this.info)))
    this.info <- this.info[match(info.cols, names(this.info))]
    names(this.info) <- info.cols
    pgxinfo1 <- pgxinfo
    for (f in setdiff(info.cols, colnames(pgxinfo1))) {
      pgxinfo1[[f]] <- NA
    }
    match(info.cols, colnames(pgxinfo1))
    pgxinfo1 <- pgxinfo1[, match(info.cols, colnames(pgxinfo1)), drop = FALSE]
    colnames(pgxinfo1) <- info.cols
    pgxinfo <- rbind(pgxinfo1, this.info)
  } else {
    pgxinfo <- data.frame(rbind(this.info))
  }
  rownames(pgxinfo) <- NULL

  return(pgxinfo)
}


#' @title Read dataset information file (aka PGX info)
#'
#' @param pgx.dir Character string specifying the path to the PGX directory containing .pgx files.
#' @param file Character string specifying the filename for the dataset info file. Default is "datasets-info.csv".
#' @param match Logical indicating whether to match .pgx filenames. Default is TRUE.
#'
#' @return Data frame containing dataset information for PGX files.
#'
#' @description Reads the dataset information CSV file for a directory of PGX files.
#'
#' @details This function reads the dataset info CSV located in a PGX directory. It contains metadata like datatype, organism, sample sizes, etc.
#'
#' The \code{pgx.dir} argument specifies the path to the PGX directory containing .pgx files.
#'
#' The \code{file} argument specifies the filename for the dataset info file. By default this is "datasets-info.csv".
#'
#' If \code{match=TRUE}, the function filters the info to only datasets matching .pgx files in the directory.
#'
#' @export
pgxinfo.read <- function(pgx.dir, file = "datasets-info.csv", match = TRUE, use.cache = FALSE, clean.description = TRUE) {
  ##  pgx.dir="~/Playground/pgx";file = "datasets-info.csv"
  ##  pgx.dir="~/Downloads";file = "datasets-info.csv"
  if (match) {
    pgx.files <- dir(pgx.dir, pattern = "[.]pgx$")
    if (length(pgx.files) == 0) {
      return(NULL)
    } ## no files!
  }

  pgxinfo <- NULL
  pgxinfo.file <- file.path(pgx.dir, file)

  ## do not use fread.csv or fread here!! see issue #441
  pgxinfo <- NULL
  if (use.cache) {
    pgxinfo <- cached.csv(pgxinfo.file, stringsAsFactors = FALSE, row.names = NULL, sep = ",")
  } else if (file.exists(pgxinfo.file)) {
    pgxinfo <- utils::read.csv(pgxinfo.file, stringsAsFactors = FALSE, row.names = NULL, sep = ",")
  }
  if (is.null(pgxinfo)) {
    return(NULL)
  }

  if ("ngenes" %in% tolower(colnames(pgxinfo)))
    colnames(pgxinfo)[tolower(colnames(pgxinfo)) == "ngenes"] <- "nfeatures"

  pgxinfo$X <- NULL ## delete first column
  pgxinfo <- pgxinfo[which(!is.na(pgxinfo$dataset)), ] ## remove NA
  if (match && nrow(pgxinfo)) {
    pgx.files1 <- sub("[.]pgx$", "", pgx.files)
    pgxinfo.datasets <- sub("[.]pgx$", "", pgxinfo$dataset)
    sel <- pgxinfo.datasets %in% pgx.files1
    pgxinfo <- pgxinfo[sel, , drop = FALSE]
  }
  info.colnames <- c(
    "dataset", "datatype", "description", "nsamples",
    "nfeatures", "nsets", "conditions", "organism", "date",
    "creator"
  )
  if (is.null(pgxinfo)) {
    aa <- rep(NA, length(info.colnames))
    names(aa) <- info.colnames
    pgxinfo <- data.frame(rbind(aa))[0, ]
  }
  ## add missing columns fields
  missing.cols <- setdiff(info.colnames, colnames(pgxinfo))
  for (s in missing.cols) pgxinfo[[s]] <- rep(NA, nrow(pgxinfo))
  ii <- match(info.colnames, colnames(pgxinfo))
  pgxinfo <- pgxinfo[, ii]

  if (clean.description) {
    ## clean up description: remove special characters from description
    ## (other columns too??)
    pgxinfo$description <- gsub("[\"\']", " ", pgxinfo$description) ## remove quotes (important!!)
    pgxinfo$description <- gsub("[\n]", ". ", pgxinfo$description) ## replace newline
    pgxinfo$description <- trimws(gsub("[ ]+", " ", pgxinfo$description)) ## remove ws
  }

  return(pgxinfo)
}

#' @title Write dataset information file (aka PGX info)
#'
#' @param pgxinfo Character string specifying the pgxinfo dataframe
#' @param pgx.dir Character string specifying the target folder
#' @param info.file Character string specifying the filename for the dataset info file. Default is "datasets-info.csv".
#'
#' @description Write the dataset information CSV file with some safety processing
#'
#' @details This function writes the dataset info CSV located in a PGX directory. It contains metadata like datatype, organism, sample sizes, etc.
#'
#' The \code{pgxinfo} argument specifies the PGXinfo dataframe
#' The \code{pgx.dir} argument specifies the target folder
#' The \code{info.file} argument specifies the filename for the dataset info file. By default this is "datasets-info.csv".
#'
#' @export
pgxinfo.write <- function(pgxinfo, pgx.dir, info.file = "datasets-info.csv") {
  info.file <- basename(info.file)
  rownames(pgxinfo) <- NULL
  pgxinfo <- data.frame(pgxinfo, check.names = FALSE)

  ## remove special characters from description (other columns too??)
  pgxinfo$description <- gsub("[\"\']", " ", pgxinfo$description) ## remove quotes (important!!)
  pgxinfo$description <- gsub("[\n]", ". ", pgxinfo$description) ## replace newline
  pgxinfo$description <- trimws(gsub("[ ]+", " ", pgxinfo$description)) ## remove ws

  file1 <- file.path(pgx.dir, info.file)
  utils::write.csv(pgxinfo, file = file1)
}


#' @title Check if PGX metadata needs update
#'
#' @param pgx.dir Path to the PGX directory containing .pgx files
#' @param check.sigdb Logical indicating whether to check the sigdb file. Default is TRUE.
#' @param verbose Logical indicating whether to print verbose messages. Default is TRUE.
#'
#' @return Logical indicating if any metadata files need update.
#'
#' @description Checks if the main PGX metadata files need to be updated against the .pgx files.
#' All pgx files in the datafolder should be included.
#'
#' @details This function checks if the main metadata files in a PGX directory need to be updated:
#'
#' - datasets-info.csv
#' - datasets-allFC.csv
#' - datasets-sigdb.h5
#'
#' It compares the .pgx files in the directory against the metadata files.
#' If any .pgx files are missing from the metadata, it will return TRUE indicating the metadata needs to be updated.
#'
#' Set verbose=TRUE to print debug messages. Set check.sigdb=FALSE to skip checking the sigdb file.
#'
#' @export
pgxinfo.needUpdate <- function(
    pgx.dir,
    check.sigdb = TRUE,
    verbose = TRUE) {
  if (!dir.exists(pgx.dir)) {
    stop(paste("[pgxinfo.needUpdate] FATAL ERROR : folder", pgx.dir, "does not exist"))
  }

  pgx.dir <- pgx.dir[1] ## only one folder!!!
  pgx.files <- dir(pgx.dir, pattern = "[.]pgx$")
  if (length(pgx.files) == 0) {
    return(FALSE)
  } ## no files!

  ## all public datasets
  pgx.files <- sub("[.]pgx$", "", pgx.files) ## strip pgx
  dbg("[pgxinfo.needUpdate] number of PGX files: ", length(pgx.files))

  ## ----------------------------------------------------------------------
  ## If an allFC file exists
  ## ----------------------------------------------------------------------
  allfc.file <- "datasets-allFC.csv"
  info.file <- "datasets-info.csv"
  sigdb.file <- "datasets-sigdb.h5"

  allfc.file1 <- file.path(pgx.dir, allfc.file)
  has.fc <- file.exists(allfc.file1)

  info.file1 <- file.path(pgx.dir, info.file)
  has.info <- file.exists(info.file1)

  sigdb.file1 <- file.path(pgx.dir, sigdb.file)
  has.sigdb <- file.exists(sigdb.file1)

  if (verbose) {
    dbg("[pgxinfo.needUpdate] has datasets-allFC.csv : ", has.fc)
    dbg("[pgxinfo.needUpdate] has datasets-info.csv  : ", has.info)
    if (check.sigdb) dbg("[pgxinfo.needUpdate] has datasets-sigdb.h5  : ", has.sigdb)
  }

  ## ----------------------------------------------------------------------
  ## If an allFC exists, check if it is done for all PGX files
  ## ----------------------------------------------------------------------

  fc.complete <- info.complete <- h5.complete <- FALSE
  fc.missing <- info.missing <- h5.missing <- pgx.files

  if (has.fc) {
    if (verbose) message("[pgxinfo.needUpdate] checking which pgx already done in allFC...")
    allFC <- data.table::fread(allfc.file1, check.names = FALSE, nrows = 1) ## HEADER!!!
    fc.files <- gsub("^\\[|\\].*", "", colnames(allFC)[-1])
    fc.files <- sub("[.]pgx$", "", fc.files) ## strip pgx
    fc.complete <- all(pgx.files %in% fc.files)
    fc.missing <- setdiff(pgx.files, fc.files)
  }

  if (has.info) {
    if (verbose) message("[pgxinfo.needUpdate] checking which pgx already in PGX info...")
    ## do not use fread! quoting bug
    pgxinfo <- utils::read.csv(info.file1, stringsAsFactors = FALSE, row.names = NULL, sep = ",")
    pgxinfo$X <- NULL ## delete first column
    pgxinfo <- pgxinfo[which(!is.na(pgxinfo$dataset)), ] ## remove NA
    info.files <- unique(sub(".pgx$", "", pgxinfo$dataset))
    info.complete <- all(pgx.files %in% info.files)
    info.missing <- setdiff(pgx.files, info.files)
    info.complete
  }

  if (has.sigdb && check.sigdb) {
    if (verbose) message("[pgxinfo.needUpdate] checking which pgx already in sigdb...")
    H <- rhdf5::h5ls(sigdb.file1)
    h5.ok1 <- all(c("matrix", "colnames", "rownames", "data") %in% H$name)
    h5.ok2 <- all(c("/clustering", "/data", "/enrichment", "/signature") %in% H$group)
    h5.ok <- (h5.ok1 && h5.ok2)
    h5.complete <- FALSE
    if (h5.ok) {
      cn <- rhdf5::h5read(sigdb.file1, "data/colnames")
      h5.files <- gsub("^\\[|\\].*", "", cn)
      h5.files <- sub("[.]pgx$", "", h5.files) ## strip pgx
      h5.complete <- all(pgx.files %in% h5.files)
      h5.missing <- setdiff(pgx.files, h5.files)
    }
  }

  ## Return checks
  has.files <- (has.fc && has.info)
  is.complete <- (fc.complete && info.complete)
  if (check.sigdb) {
    has.files <- has.files && has.sigdb
    is.complete <- is.complete && h5.complete
  }
  return(!has.files || !is.complete)
}


#' @title Update PGX dataset folder metadata
#'
#' @param pgx.dir Path to the PGX dataset folder containing .pgx files
#' @param force Logical indicating whether to force update even if not needed. Default is FALSE.
#' @param delete.old Logical indicating whether to delete old metadata files. Default is FALSE.
#' @param new.pgx Character vector of new .pgx files added to the folder. Default is NULL.
#' @param update.sigdb Logical indicating whether to update sigdb file. Default is TRUE.
#' @param verbose Logical indicating whether to print status messages. Default is TRUE.
#'
#' @return NULL. Metadata files in pgx.dir are updated if needed.
#'
#' @description Updates the main metadata files in a PGX dataset folder if needed.
#'
#' @details This function checks if the main metadata files in a PGX folder need to be updated
#' based on new or changed .pgx files. It compares the .pgx files to the existing metadata.
#'
#' The metadata files updated are:
#' - datasets-info.csv
#' - datasets-allFC.csv
#' - datasets-sigdb.h5
#'
#' If force=TRUE, it will update regardless of whether changes are detected.
#' If delete.old=TRUE, existing metadata files will be deleted before writing new files.
#'
#' @export
pgxinfo.updateDatasetFolder <- function(pgx.dir,
                                        force = FALSE,
                                        delete.old = FALSE,
                                        new.pgx = NULL,
                                        update.sigdb = TRUE,
                                        verbose = TRUE) {
  if (!dir.exists(pgx.dir)) {
    stop(paste("[pgxinfo.updateDatasetFolder] FATAL ERROR : folder", pgx.dir, "does not exist"))
  }

  pgx.dir <- pgx.dir[1] ## only one folder!!!
  pgx.files <- dir(pgx.dir, pattern = "[.]pgx$")
  pgx.files <- sub("[.]pgx$", "", pgx.files) ## strip pgx

  allfc.file <- file.path(pgx.dir, "datasets-allFC.csv")
  info.file <- file.path(pgx.dir, "datasets-info.csv")
  sigdb.file <- file.path(pgx.dir, "datasets-sigdb.h5")
  tsne.file <- file.path(pgx.dir, "datasets-tsne.csv")

  if (length(pgx.files) == 0) {
    file.remove(info.file)
    file.remove(allfc.file)
    # should return FALSE?
    return(NULL)
  }

  ## ----------------------------------------------------------------------
  ## subroutines
  ## ----------------------------------------------------------------------
  updateTSNE <- function(tsne.file, sigdb.file, allFC) {
    h5 <- rhdf5::h5ls(sigdb.file)
    has.tsne2d <- ("tsne2d" %in% h5$name)
    if (has.tsne2d) {
      cn <- rhdf5::h5read(sigdb.file, "data/colnames")
      tsne <- rhdf5::h5read(sigdb.file, "clustering/tsne2d")
      if (length(cn) != nrow(tsne)) {
        dbg("[pgxinfo.updateDatasetFolder] ***ERROR**** length mismatch!")
      }
      rownames(tsne) <- cn
    } else {
      tsne <- pgx.clusterBigMatrix(
        as.matrix(allFC),
        methods = "tsne",
        dims = 2,
        perplexity = 50,
        center.features = TRUE,
        scale.features = FALSE,
        reduce.sd = 1000,
        reduce.pca = 50,
        find.clusters = FALSE
      )$tsne2d
    }
    colnames(tsne) <- paste0("tsne.", 1:ncol(tsne))
    utils::write.csv(tsne, file = tsne.file)
  }

  ## ----------------------------------------------------------------------
  ## If an allFC file exists
  ## ----------------------------------------------------------------------

  has.fc <- file.exists(allfc.file)
  has.info <- file.exists(info.file)
  has.sigdb <- file.exists(sigdb.file)
  has.tsne <- file.exists(tsne.file)

  if (force) {
    has.fc <- FALSE
    has.info <- FALSE
    has.sigdb <- FALSE
    has.tsne <- FALSE
  }

  ## ----------------------------------------------------------------------
  ## Early check if signature/sig100.up/down is formated correctly
  ## From Axel we see some weird cases where it gets broken
  ## ----------------------------------------------------------------------
  if (has.sigdb) {
    h5 <- rhdf5::h5ls(sigdb.file)
    has.sig100.up <- ("sig100.up" %in% h5$name)
    has.sig100.dn <- ("sig100.dn" %in% h5$name)
    if (has.sig100.up && has.sig100.dn) {
      sig100.up <- rhdf5::h5read(sigdb.file, "signature/sig100.up")
      sig100.dn <- rhdf5::h5read(sigdb.file, "signature/sig100.dn")
      ## Check if sig100.up/down are properly formatted (should be character gene names, not numbers)
      if (is.numeric(sig100.up) || is.numeric(sig100.dn)) {
        dbg("[pgxinfo.updateDatasetFolder] ***WARNING*** sig100.up/dn are numeric, should be character gene names")
        file.remove(sigdb.file)
        has.sigdb <- FALSE ## force recreation
      } else if (is.character(sig100.up) && is.character(sig100.dn)) {
        ## Check if they are numbers formatted as characters
        up.numeric <- suppressWarnings(!is.na(as.numeric(sig100.up[1, 1])))
        dn.numeric <- suppressWarnings(!is.na(as.numeric(sig100.dn[1, 1])))
        if (up.numeric || dn.numeric) {
          dbg("[pgxinfo.updateDatasetFolder] ***WARNING*** sig100.up/dn contain numeric values as characters, should be gene names")
          file.remove(sigdb.file)
          has.sigdb <- FALSE ## force recreation
        }
      }
    } else {
      has.sigdb <- FALSE ## force recreation
    }
  }

  ## ----------------------------------------------------------------------
  ## Check if all files are missing or have entries to be deleted
  ## ----------------------------------------------------------------------

  pgxinfo <- NULL
  fc.missing <- pgx.files
  info.missing <- pgx.files
  h5.missing <- pgx.files
  fc.delete <- c()
  info.delete <- c()
  h5.delete <- c()

  allFC <- NULL
  fc.files <- NULL
  if (has.fc) {
    ## read just header
    allFC <- read.csv(allfc.file, check.names = FALSE, nrows = 1)
    colnames(allFC) <- gsub("[\"\'\n\t]", "", colnames(allFC))

    fc.files <- gsub("\\].*", "", colnames(allFC)[-1])
    fc.files <- gsub("^.*\\[", "", fc.files)
    fc.missing <- setdiff(pgx.files, fc.files)
    fc.delete <- setdiff(fc.files, pgx.files)
    allFC <- NULL
  }

  pgxinfo.files <- NULL
  if (has.info) {
    ## do not use fread! quoting bug
    pgxinfo <- utils::read.csv(info.file, stringsAsFactors = FALSE, row.names = NULL, sep = ",")
    pgxinfo$X <- NULL ## delete first column
    pgxinfo <- pgxinfo[which(!is.na(pgxinfo$dataset)), ] ## remove NA
    pgxinfo.files <- unique(sub(".pgx$", "", pgxinfo$dataset))
    info.missing <- setdiff(pgx.files, pgxinfo.files)
    info.delete <- setdiff(pgxinfo.files, pgx.files)
  }

  valid.h5 <- function(h5.file) {
    if (!file.exists(h5.file)) {
      return(FALSE)
    }
    H <- rhdf5::h5ls(h5.file)
    ok1 <- all(c("matrix", "colnames", "rownames", "data") %in% H[, "name"])
    ok2 <- all(c("/clustering", "/data", "/enrichment", "/signature") %in% H[, "group"])
    ok1 && ok2
  }
  valid.h5(sigdb.file)

  h5.files <- NULL
  if (has.sigdb && valid.h5(sigdb.file)) {
    cn <- rhdf5::h5read(sigdb.file, "data/colnames")
    h5.files <- gsub("\\].*", "", cn)
    h5.files <- gsub("^.*\\[", "", h5.files)
    h5.files <- sub("[.]pgx$", "", h5.files) ## strip pgx
    h5.files
    h5.missing <- setdiff(pgx.files, h5.files)
    h5.delete <- setdiff(h5.files, pgx.files)
  }

  tsne.files <- NULL
  tsne.missing <- pgx.files
  if (file.exists(tsne.file)) {
    tsne <- read.csv(tsne.file, row.names = 1)
    tsne.files <- gsub("^\\[|\\].*", "", rownames(tsne))
    tsne.missing <- setdiff(pgx.files, tsne.files)
  }

  all(pgx.files %in% pgxinfo.files)
  all(pgx.files %in% fc.files)
  all(pgx.files %in% h5.files)
  all(pgx.files %in% tsne.files)

  ## ----------------------------------------------------------------------
  ## Check if it is done for all PGX files
  ## ----------------------------------------------------------------------

  ## files to be done either for allFC or missing in INFO
  pgx.missing <- unique(c(fc.missing, info.missing))
  pgx.delete <- unique(c(fc.delete, info.delete))

  pgxinfo.changed <- FALSE
  pgxfc.changed <- FALSE

  ## these pgx need forced update
  if (!is.null(new.pgx)) {
    new.pgx <- sub(".pgx$", "", new.pgx)
    new.pgx <- intersect(new.pgx, pgx.files) ## only existing pgx
    pgx.delete <- union(pgx.delete, new.pgx)
    pgx.missing <- union(pgx.missing, new.pgx)
    fc.missing <- union(fc.missing, new.pgx)
    info.missing <- union(info.missing, new.pgx)
  }

  ## nothing to do???
  if (length(pgx.missing) == 0 && length(pgx.delete) == 0 && length(h5.missing) == 0 &&
    length(tsne.missing) == 0) {
    if (verbose) message("[pgxinfo.updateDatasetFolder] All files complete. No update required.")
    return()
  }

  ## pgxinfo and allFC OK, but sigdb not upto date
  if (update.sigdb && length(pgx.missing) == 0 && length(pgx.delete) == 0 &&
    (length(h5.missing) || length(tsne.missing))) {
    allFC <- fread.csv(allfc.file, row.names = 1, check.names = FALSE)
    allFC <- as.matrix(allFC)
    if (length(h5.missing)) {
      if (verbose) message("[updateDatasetFolder] missing sigdb. Creating new sigdb file.")
      if (file.exists(sigdb.file)) unlink(sigdb.file)
      pgx.createSignatureDatabaseH5.fromMatrix(sigdb.file, X = allFC)
    }
    if (length(tsne.missing)) {
      if (verbose) message("[updateDatasetFolder] missing tsne. Creating new tsne file.")
      updateTSNE(tsne.file, sigdb.file, allFC)
    }
    return()
  }

  ## ----------------------------------------------------------------------
  ## Remove to-be-deleted entries
  ## ----------------------------------------------------------------------
  allFC <- NULL

  ## Reread allFC file. Before we only read the header.
  if (!force && file.exists(allfc.file)) {
    allFC <- fread.csv(allfc.file, row.names = 1, check.names = FALSE)
    allFC <- as.matrix(allFC)
    colnames(allFC) <- gsub("[\"\']", "", colnames(allFC)) ## safety
    if (ncol(allFC) == 0) allFC <- NULL
  }

  ## remove from allFC matrix
  if (!is.null(allFC) && length(pgx.delete)) {
    allfc.pgx <- gsub("^.*\\[|\\].*", "", colnames(allFC))
    del <- which(allfc.pgx %in% pgx.delete)
    if (length(del)) {
      allFC <- allFC[, -del, drop = FALSE]
      pgxfc.changed <- TRUE
    }
  }

  ## remove pgx entry from  info file
  if (!is.null(pgxinfo) && length(pgx.delete)) {
    del <- which(pgxinfo$dataset %in% pgx.delete)
    if (length(del)) {
      pgxinfo <- pgxinfo[-del, , drop = FALSE]
      pgxinfo.changed <- TRUE
    }
  }

  ## remove from H5 sigdb file
  if (has.sigdb && length(h5.delete)) {
    cn <- rhdf5::h5read(sigdb.file, "data/colnames")
    h5.files <- gsub("^\\[|\\].*", "", cn)
    h5.files <- sub("[.]pgx$", "", h5.files) ## strip pgx
    del <- which(h5.files %in% h5.delete)
    if (length(del)) {
      ##      cn[del] <- paste("[DELETED]", sub(".*\\] ", "", cn[del]))
      cn[del] <- paste("[DELETED]", cn[del])
      rhdf5::h5delete(sigdb.file, "data/colnames")
      rhdf5::h5write(cn, sigdb.file, "data/colnames")
      if (delete.old) {
        sigdb.removeDataset(sigdb.file, "DELETED")
        dbg("[pgxinfo.updateDatasetFolder] deleting H5 entries")
      }
    }
  }

  ## ----------------------------------------------------------------------
  ## For all new PGX files, load the PGX file and get the meta FC
  ## matrix.
  ## ----------------------------------------------------------------------
  info.cols <- NULL
  missing.FC <- list()

  if (length(pgx.missing)) {
    dbg("[updateDatasetFolder] missing pgx = ", pgx.missing)
    pgxfile <- pgx.missing[1]
    pgx <- NULL
    for (pgxfile in pgx.missing) {
      cat(".")

      ## load the PGX file
      pgxfile1 <- file.path(pgx.dir, pgxfile)
      pgxfile1 <- paste0(sub("[.]pgx$", "", pgxfile1), ".pgx")
      pgx <- try(local(get(load(pgxfile1, verbose = 0)))) ## override any name

      if (inherits(pgx, "try-error")) {
        message(paste("[updateDatasetFolder] ERROR in loading PGX file:", pgxfile1, ". skipping\n"))
        next()
      }

      if (!pgx.checkObject(pgx)) {
        message(paste("[updateDatasetFolder] INVALID PGX object", pgxfile, ". Skipping"))
        next()
      }

      ## ---------------------------------------------
      ## extract the meta FC matrix
      ## ---------------------------------------------
      if (pgxfile %in% fc.missing) {
        meta <- pgx.getMetaFoldChangeMatrix(pgx, what = "meta")
        ## rownames(meta$fc) <- toupper(rownames(meta$fc)) ## human genes
        F <- collapse_by_humansymbol(meta$fc, pgx$genes)
        if (nrow(F) == 1 && rownames(F) == "---") {
          F <- meta$fc
        }
        missing.FC[[pgxfile]] <- F
        pgxfc.changed <- TRUE
      }

      ## ---------------------------------------------
      ## compile the info for update
      ## ---------------------------------------------
      if (pgxfile %in% info.missing) {
        pgx$name <- sub(".pgx$", "", pgxfile) ## force filename as name
        pgxinfo <- pgxinfo.add(pgxinfo, pgx)
        pgxinfo.changed <- TRUE
      }
    }
    remove(pgx)
  }

  ## ----------------------------------------------------------------------
  ## Update the INFO meta file
  ## ----------------------------------------------------------------------
  if (pgxinfo.changed) {
    dbg("[pgxinfo.updateDatasetFolder] updating pgxinfo file")
    pgxinfo.write(pgxinfo, pgx.dir = pgx.dir, info.file = basename(info.file))
    Sys.chmod(info.file, "0666")
  }

  ## ----------------------------------------------------------------------
  ## Update the ALL.FC meta file
  ## ----------------------------------------------------------------------
  pgxfc.changed
  if (pgxfc.changed && length(missing.FC) > 0) {
    ## find most common genes
    all.gg <- toupper(as.character(unlist(sapply(missing.FC, rownames))))
    gg.tbl <- table(all.gg)

    ## Conform the multiple metaFC matrices
    gg <- names(gg.tbl)

    missing.FC <- lapply(missing.FC, function(x) {
      x <- x[match(gg, toupper(rownames(x))), , drop = FALSE]
      rownames(x) <- gg
      return(x)
    })

    ## append file name in front of contrast names
    id <- paste0("[", sub("[.]pgx$", "", names(missing.FC)), "]")
    id
    for (i in 1:length(missing.FC)) {
      ## clean up string
      fcnames <- colnames(missing.FC[[i]])
      fcnames <- gsub("[\"\']", " ", fcnames) ## remove quotes (important!!)
      fcnames <- gsub("[\n]", ". ", fcnames) ## replace newline
      fcnames <- trimws(gsub("[ ]+", " ", fcnames)) ## remove ws
      colnames(missing.FC[[i]]) <- paste0(id[i], " ", fcnames)
    }
    allFC.new <- do.call(cbind, missing.FC)
    allFC.new <- as.matrix(allFC.new)

    if (is.null(allFC)) {
      allFC <- allFC.new
    } else {
      ## Add any new FC profiles to the existing allFC
      gg <- sort(unique(c(rownames(allFC), rownames(allFC.new))))
      j1 <- match(gg, rownames(allFC))
      j2 <- match(gg, rownames(allFC.new))
      allFC <- allFC[j1, , drop = FALSE]
      allFC.new <- allFC.new[j2, , drop = FALSE]
      allFC <- cbind(allFC, allFC.new)
      rownames(allFC) <- gg
    }

    ## restrict to 20000 genes
    ## allfc.sd <- apply(allFC, 1, stats::sd, na.rm = TRUE)
    allfc.sd <- matrixStats::rowSds(allFC, na.rm = TRUE)
    allfc.nna <- rowMeans(!is.na(allFC))
    jj <- Matrix::head(order(-allfc.sd * allfc.nna), 20000)
    allFC <- allFC[jj, , drop = FALSE]
    pgxfc.changed <- TRUE
  }

  ## save modified allFC
  if (pgxfc.changed) {
    ## check for duplicates
    dbg("[pgxinfo.updateDatasetFolder] Updating allFC file")
    allFC <- allFC[, !duplicated(colnames(allFC)), drop = FALSE]
    allFC <- allFC[, order(colnames(allFC)), drop = FALSE]
    allFC <- round(allFC, digits = 4)
    AA <- data.frame(rownames = rownames(allFC), allFC, check.names = FALSE)
    data.table::fwrite(AA, file = allfc.file)
    Sys.chmod(allfc.file, "0666")
    remove(AA)
  }

  ## update user sigdb or create if not exists
  update.sigdb2 <- (!file.exists(sigdb.file) || pgxfc.changed || length(h5.missing))
  update.sigdb2
  if (update.sigdb2) {
    dbg("[pgxinfo.updateDatasetFolder] Updating sigdb file")

    if (!file.exists(sigdb.file) || force) {
      ## NEED RETHINK!!!! HERE???
      if (file.exists(sigdb.file)) unlink(sigdb.file)
      dbg("[pgxinfo.updateDatasetFolder] creating signature DB")
      pgx.createSignatureDatabaseH5.fromMatrix(sigdb.file, X = allFC)
      dbg("[pgxinfo.updateDatasetFolder] done. created sigdb.file = ", sigdb.file)
    } else if (pgxfc.changed) {
      cn <- rhdf5::h5read(sigdb.file, "data/colnames")
      fc.add <- any(!colnames(allFC) %in% cn)
      fc.del <- any(!cn %in% colnames(allFC))
      if (fc.add) {
        pgx.createSignatureDatabaseH5.fromMatrix(sigdb.file, X = allFC)
        cn <- rhdf5::h5read(sigdb.file, "data/colnames") ## new size
      }

      if (fc.del) {
        del <- which(!cn %in% colnames(allFC))
        ##  cn[del] <- paste("[DELETED]", sub(".*\\] ", "", cn[del]))
        cn[del] <- paste("[DELETED]", cn[del])
        rhdf5::h5delete(sigdb.file, "data/colnames")
        rhdf5::h5write(cn, sigdb.file, "data/colnames")
        Sys.sleep(0.1) ## time to write...
        if (delete.old) {
          dbg("[pgxinfo.updateDatasetFolder] deleting sigDB entry")
          sigdb.removeDataset(sigdb.file, "DELETED")
          Sys.sleep(0.1) ## allow write time
        } else {
          dbg("[pgxinfo.updateDatasetFolder] deletions tagged but not removed")
        }
      }
    } else {
      ## no change
    }
  }

  ## -------------------------------------------------------
  ## check tsne file. update tsne file from H5
  ## -------------------------------------------------------
  tsne.ok <- FALSE
  if (file.exists(tsne.file)) {
    dbg("[pgxinfo.updateDatasetFolder] checking TSNE file...")
    tsne <- read.csv(tsne.file, row.names = 1)
    dim(tsne)
    cn <- rhdf5::h5read(sigdb.file, "data/colnames") ## new size
    cn <- grep("DELETED", cn, invert = TRUE)
    tsne.ok <- ifelse(all(cn %in% rownames(tsne)), TRUE, FALSE)
  }

  ## if needed update tsne file from H5
  if (!tsne.ok || !file.exists(tsne.file) || pgxfc.changed) {
    dbg("[pgxinfo.updateDatasetFolder] updating TSNE...")
    updateTSNE(tsne.file, sigdb.file, allFC)
  }

  dbg("[pgxinfo.updateDatasetFolder] done!")
  return()
}

#' Delete pgx entries in pgx info objects
#'
#' @param pgx.dir The folder containing pgxinfo metadata
#' @param pgxname The name of the pgx object for which metadata to delete
#' @param purge.h5 Logical indicating whether to remove entry in big H5 file
#'
#' @return NULL
#'
#' @description
#' Removes entries the pgxinfo metadata files from a pgxname
#'
#' @details
#' This function takes a pgxname as input and removes all entries in the pgx metadata files.
#'
#'
#' @export
pgxinfo.delete <- function(pgx.dir, pgxname, purge.h5 = FALSE) {
  allfc.file <- file.path(pgx.dir, "datasets-allFC.csv")
  info.file <- file.path(pgx.dir, "datasets-info.csv")
  sigdb.file <- file.path(pgx.dir, "datasets-sigdb.h5")
  tsne.file <- file.path(pgx.dir, "datasets-tsne.csv")

  allFC <- NULL

  ## Reread allFC file. Before we only read the header.
  if (file.exists(allfc.file)) {
    allFC <- read.csv(allfc.file, row.names = 1, nrow = 1, check.names = FALSE)
    allFC <- as.matrix(allFC)
    allfc.pgx <- gsub("^\\[|\\].*", "", colnames(allFC))
    del <- which(allfc.pgx == pgxname)
    if (length(del)) {
      allFC <- fread.csv(allfc.file, row.names = 1, check.names = FALSE)
      allFC <- allFC[, -del, drop = FALSE]
      allFC <- round(allFC, digits = 4)
      AA <- data.frame(rownames = rownames(allFC), allFC, check.names = FALSE)
      data.table::fwrite(AA, file = allfc.file)
      Sys.chmod(allfc.file, "0666")
      remove(AA)
    }
  }

  ## remove from PGX info file
  if (file.exists(info.file)) {
    pgxinfo <- utils::read.csv(info.file, stringsAsFactors = FALSE, row.names = NULL, sep = ",")
    pgxinfo$X <- NULL ## delete first column
    pgxinfo <- pgxinfo[which(!is.na(pgxinfo$dataset)), ] ## remove NA
    del <- which(pgxinfo$dataset == pgxname)
    if (length(del)) {
      pgxinfo <- pgxinfo[-del, , drop = FALSE]
      rownames(pgxinfo) <- NULL
      pgxinfo <- data.frame(pgxinfo, check.names = FALSE)
      utils::write.csv(pgxinfo, file = info.file)
    }
  }

  ## remove from H5 sigdb file
  if (file.exists(sigdb.file)) {
    cn <- rhdf5::h5read(sigdb.file, "data/colnames")
    h5.files <- gsub("^\\[|\\].*", "", cn)
    h5.files <- sub("[.]pgx$", "", h5.files) ## strip pgx
    del <- which(h5.files == pgxname)
    if (length(del)) {
      ##      cn[del] <- paste("[DELETED]", sub(".*\\] ", "", cn[del]))
      cn[del] <- paste("[DELETED]", cn[del])
      rhdf5::h5delete(sigdb.file, "data/colnames")
      rhdf5::h5write(cn, sigdb.file, "data/colnames")
    }
    if (purge.h5) {
      sigdb.removeDataset(sigdb.file, "DELETED")
    }
  }

  ## remove from PGX info file
  if (file.exists(tsne.file)) {
    pos <- utils::read.csv(tsne.file, stringsAsFactors = FALSE, row.names = 1)
    pos.files <- gsub("^\\[|\\].*", "", rownames(pos))
    pos.files <- sub("[.]pgx$", "", pos.files) ## strip pgx
    del <- which(pos.files %in% c(pgxname, "DELETED"))
    if (length(del)) {
      pos <- pos[-del, , drop = FALSE]
      utils::write.csv(pos, file = tsne.file)
    }
  }
}
