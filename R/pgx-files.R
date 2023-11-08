##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' @title Load an R Object from a File
#'
#' @description This function loads an R object from a file and makes it available in the local environment.
#'
#' @param file A character string specifying the path to the file containing the R object to be loaded.
#' @param verbose An integer specifying the level of verbosity for printing messages during the loading process (default is 0).
#'
#' @details The function uses the `load` function from base R to load the R object from the specified file.
#' The loaded object is then made available in the local environment using the `local` function.
#'
#' @return The function does not return a value, but makes the loaded R object available in the local environment.
#'
#' @export
pgx.load <- function(file, verbose = 0) {
  local(get(load(file, verbose = verbose)))
}


#' @describeIn pgx.save deprecated version to save pgx/ngs file.
#' @export
ngs.save <- function(ngs, file, update.date = TRUE, light = TRUE, system = FALSE) {
  message("WARNING: ngs.save() is deprecated. please use pgx.save()")
  pgx.save(ngs, file = file, update.date = update.date, light = light, system = system)
}


#' Save PGX object to file
#'
#' @title Save PGX Object
#'
#' @param pgx PGX object to save
#' @param file File path to save PGX object to
#' @param update.date Logical indicating whether to update date field. Default TRUE.
#' @param light Logical indicating whether to save a light version without some large fields. Default TRUE.
#' @param system Logical indicating whether to keep system-level objects like omicsnet. Default FALSE.
#'
#' @return NULL. The PGX object is saved to the specified file path.
#'
#' @description Saves a PGX object to an R binary file.
#'
#' @details This function saves a PGX object to a .rda or .RData file.
#' By default it saves a light version without some large objects like gene set matrices.
#' It also optionally updates the date field and removes system-level objects like the omicsnet.
#' The PGX object can later be loaded back into R using \code{pgx.load()}.
#'
#' @export
pgx.save <- function(pgx, file, update.date = TRUE, light = TRUE, system = FALSE) {
  if (update.date || is.null(pgx$date)) pgx$date <- Sys.Date()

  if (light) {
    ## ------- make a light version
    pgx$gx.meta$outputs <- NULL
    pgx$gset.meta$outputs <- NULL
    pgx$model.parameters$efit <- NULL
    pgx$gmt.all <- NULL
    pgx$collections <- NULL
    pgx$gset.meta$matrices <- NULL
  }
  if (system == FALSE) {
    ## remove system (is big...)
    pgx$omicsnet <- NULL
    pgx$omicsnet.reduced <- NULL
  }
  sort(sapply(pgx, utils::object.size)) / 1e9
  sum(sapply(pgx, utils::object.size)) / 1e9

  cat(">>> saving PGX file to", file, "\n")
  file <- iconv(file, from = "", to = "ASCII//TRANSLIT")
  save(pgx, file = file)
}


#' @title Check if an Object Exists in an HDF5 File
#'
#' @description This function checks if a specified object exists in an HDF5 file.
#'
#' @param h5.file A character string specifying the path to the HDF5 file.
#' @param obj A character string specifying the name of the object to check for existence in the HDF5 file.
#'
#' @details The function uses the `h5ls` function from the rhdf5 package to list all objects in the specified HDF5 file.
#' It then constructs the full path for each object by concatenating the group and name columns of the resulting data frame, separated by a forward slash.
#' The function checks if the specified object name is present in this list of full paths, ignoring any leading forward slashes.
#'
#' @return A logical value indicating whether the specified object exists in the HDF5 file (TRUE) or not (FALSE).
#'
#' @export
h5exists <- function(h5.file, obj) {
  xobjs <- apply(rhdf5::h5ls(h5.file)[, 1:2], 1, paste, collapse = "/")
  obj %in% gsub("^/|^//", "", xobjs)
}


#' @title Save matrix to HDF5
#'
#' @param X The matrix to save
#' @param h5.file Path to the HDF5 file
#' @param chunk Chunk size for chunked storage. Default NULL for no chunking.
#'
#' @return NULL. The matrix is saved to the HDF5 file.
#'
#' @description Saves a matrix to an HDF5 file for efficient storage and retrieval.
#'
#' @details This function saves a matrix \code{X} to an HDF5 file at \code{h5.file}.
#' It first deletes any existing file at that path, then creates a new HDF5 file.
#'
#' The matrix is saved under the "data/matrix" group. Chunked storage can be used
#' by specifying a \code{chunk} size. This allows efficient access of subsets of
#' the matrix.
#'
#' The matrix is saved using lossless compression level 7. The HDF5 file remains
#' open after writing, and should be closed using \code{rhdf5::h5close()} after use.
#'
#' @export
pgx.saveMatrixH5 <- function(X, h5.file, chunk = NULL) {
  if (file.exists(h5.file)) unlink(h5.file)
  rhdf5::h5createFile(h5.file)
  rhdf5::h5createGroup(h5.file, "data")
  X <- as.matrix(X)

  if (is.null(chunk)) {
    rhdf5::h5write(X, h5.file, "data/matrix")
  } else {
    rhdf5::h5createDataset(
      h5.file, "data/matrix",
      c(nrow(X), ncol(X)),
      chunk = chunk,
      level = 7
    )
    rhdf5::h5write(
      X,
      file = h5.file,
      name = "data/matrix",
      index = list(1:nrow(X), 1:ncol(X))
    )
  }
  rhdf5::h5write(rownames(X), h5.file, "data/rownames")
  rhdf5::h5write(colnames(X), h5.file, "data/colnames")

  rhdf5::h5closeAll()
}


#' @title Read PGX Options
#'
#' @param file The path to the PGX options file. Default is "./OPTIONS".
#'
#' @return A named list containing the PGX options.
#'
#' @description Reads PGX analysis options from a file.
#'
#' @details This function reads a simple text file containing PGX analysis options,
#' one option per line in the format:
#'
#' \code{option=value}
#'
#' Options include parameters like:
#'
#' \code{fdr=0.05} - FDR threshold
#' \code{logfc=1} - Log fold-change threshold
#'
#' The options file allows saving a set of parameters for easily re-running
#' an analysis with the same settings.
#'
#' @export
pgx.readOptions <- function(file = "./OPTIONS") {
  if (!file.exists(file)) {
    return(NULL)
  }
  P <- utils::read.table(file, sep = "=", row.names = 1)
  opt.names <- trimws(rownames(P))
  opt <- list(P[, 1])
  opt <- sapply(opt, trimws)
  opt <- as.list(opt)
  names(opt) <- opt.names
  opt <- sapply(opt, strsplit, split = "[;]")
  ## convert character to R types
  opt <- lapply(opt, utils::type.convert, as.is = TRUE)
  opt
}


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

  organism <- NULL
  if ("organism" %in% names(pgx)) {
    organism <- pgx$organism
  } else {
    organism <- pgx.getOrganism(pgx$counts)
  }

  this.date <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  date <- ifelse(is.null(pgx$date), this.date, as.character(pgx$date))
  dataset.name <- pgx$name

  creator <- ifelse("creator" %in% names(pgx), pgx$creator, "")

  this.info <- c(
    dataset = dataset.name,
    ## author = "", ## add author? maintainer? owner??
    creator = creator,
    datatype = ifelse(is.null(pgx$datatype), "", pgx$datatype),
    description = ifelse(is.null(pgx$description), "", pgx$description),
    organism = organism,
    nsamples = nrow(pgx$samples),
    ngenes = nrow(pgx$X),
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
pgxinfo.read <- function(pgx.dir, file = "datasets-info.csv", match = TRUE) {
  ##  pgx.dir="~/Playground/pgx";file = "datasets-info.csv"
  ##  pgx.dir="~/Downloads";file = "datasets-info.csv"
  pgx.files <- dir(pgx.dir, pattern = "[.]pgx$")
  if (length(pgx.files) == 0) {
    return(NULL)
  } ## no files!

  pgxinfo <- NULL
  pgxinfo.file <- file.path(pgx.dir, file)
  if (file.exists(pgxinfo.file)) {
    ## do not use fread.csv or fread here!! see issue #441
    pgxinfo <- utils::read.csv(pgxinfo.file, stringsAsFactors = FALSE, row.names = NULL, sep = ",")
    pgxinfo$X <- NULL ## delete first column
    pgxinfo <- pgxinfo[which(!is.na(pgxinfo$dataset)), ] ## remove NA
    if (match && nrow(pgxinfo)) {
      pgx.files1 <- sub("[.]pgx$", "", pgx.files)
      pgxinfo.datasets <- sub("[.]pgx$", "", pgxinfo$dataset)
      sel <- pgxinfo.datasets %in% pgx.files1
      pgxinfo <- pgxinfo[sel, , drop = FALSE]
    }
  }
  info.colnames <- c(
    "dataset", "datatype", "description", "nsamples",
    "ngenes", "nsets", "conditions", "organism", "date",
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
  return(pgxinfo)
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
    allfc.file1 <- file.path(pgx.dir, allfc.file)
    info.file1 <- file.path(pgx.dir, info.file)
    file.remove(info.file1)
    file.remove(allfc.file1)
    # should return FALSE?
    return(NULL)
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
  if (has.fc) {
    allFC <- data.table::fread(allfc.file, check.names = FALSE, nrows = 1) ## HEADER!!!
    fc.files <- gsub("^\\[|\\].*", "", colnames(allFC)[-1])
    fc.missing <- setdiff(pgx.files, fc.files)
    fc.delete <- setdiff(fc.files, pgx.files)
    allFC <- NULL
  }

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

  if (has.sigdb && valid.h5(sigdb.file)) {
    cn <- rhdf5::h5read(sigdb.file, "data/colnames")
    h5.files <- gsub("^\\[|\\].*", "", cn)
    h5.files <- sub("[.]pgx$", "", h5.files) ## strip pgx
    h5.missing <- setdiff(pgx.files, h5.files)
    h5.delete <- setdiff(h5.files, pgx.files)
  }

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

  dbg("[pgxinfo.updateDatasetFolder] pgx.missing = ", pgx.missing)
  dbg("[pgxinfo.updateDatasetFolder] pgx.delete = ", pgx.delete)
  dbg("[pgxinfo.updateDatasetFolder] h5.missing = ", h5.missing)
  dbg("[pgxinfo.updateDatasetFolder] h5.delete = ", h5.delete)

  ## nothing to do???
  if (length(pgx.missing) == 0 && length(pgx.delete) == 0 && length(h5.missing) == 0) {
    if (verbose) message("[pgxinfo.updateDatasetFolder] All files complete. No update required.")
    return()
  }

  ## pgxinfo and allFC OK, but sigdb not upto date
  if (update.sigdb && length(pgx.missing) == 0 && length(pgx.delete) == 0 &&
    length(h5.missing) > 0) {
    allFC <- fread.csv(allfc.file, row.names = 1, check.names = FALSE)
    allFC <- as.matrix(allFC)
    if (file.exists(sigdb.file)) unlink(sigdb.file)
    if (verbose) message("[updateDatasetFolder] missing sigdb. Creating new sigdb file.")
    pgx.createSignatureDatabaseH5.fromMatrix(sigdb.file, X = allFC)
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
  }

  ## remove from allFC matrix
  if (!is.null(allFC) && length(pgx.delete)) {
    allfc.pgx <- gsub("^\\[|\\].*", "", colnames(allFC))
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
        rownames(meta$fc) <- toupper(rownames(meta$fc))
        missing.FC[[pgxfile]] <- meta$fc
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
    rownames(pgxinfo) <- NULL
    pgxinfo <- data.frame(pgxinfo, check.names = FALSE)
    utils::write.csv(pgxinfo, file = info.file)
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
    id <- paste0("[", sub("[.]pgx", "", names(missing.FC)), "]")
    id
    for (i in 1:length(missing.FC)) {
      colnames(missing.FC[[i]]) <- paste0(id[i], " ", colnames(missing.FC[[i]]))
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
    allfc.sd <- apply(allFC, 1, stats::sd, na.rm = TRUE)
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
  update.sigdb2 <- (!file.exists(sigdb.file)) || pgxfc.changed
  update.sigdb2
  if (update.sigdb2) {
    dbg("[pgxinfo.updateDatasetFolder] Updating sigdb file")

    if (!file.exists(sigdb.file) || force) {
      ## NEED RETHINK!!!! HERE???
      if (file.exists(sigdb.file)) unlink(sigdb.file)
      dbg("[pgxinfo.updateDatasetFolder] creating signature DB")
      pgx.createSignatureDatabaseH5.fromMatrix(sigdb.file, X = allFC)
    } else if (pgxfc.changed) {
      cn <- rhdf5::h5read(sigdb.file, "data/colnames")
      fc.add <- any(!colnames(allFC) %in% cn)
      fc.del <- any(!cn %in% colnames(allFC))
      if (fc.add) {
        pgx.createSignatureDatabaseH5.fromMatrix(sigdb.file, X = allFC)
        cn <- rhdf5::h5read(sigdb.file, "data/colnames")  ## new size
      }      
      
      if (fc.del) {
        del <- which(!cn %in% colnames(allFC))
        ##        cn[del] <- paste("[DELETED]", sub(".*\\] ", "", cn[del]))
        cn[del] <- paste("[DELETED]", cn[del])
        rhdf5::h5delete(sigdb.file, "data/colnames")
        rhdf5::h5write(cn, sigdb.file, "data/colnames")
        Sys.sleep(0.1)  ## time to write...
        if (delete.old) {
          dbg("[pgxinfo.updateDatasetFolder] deleting sigDB entry")
          sigdb.removeDataset(sigdb.file, "DELETED")
          Sys.sleep(0.1)  ## allow write time
        } else {
          dbg("[pgxinfo.updateDatasetFolder] deletions tagged but not removed")
        }
      }
    } else {
      ## no change
    }

    ## update tsne file from H5
    tsne.file <- file.path(pgx.dir, "datasets-tsne.csv")
    if (!file.exists(tsne.file) || pgxfc.changed) {
      cn <- rhdf5::h5read(sigdb.file, "data/colnames")
      tsne <- rhdf5::h5read(sigdb.file, "clustering/tsne2d")
      length(cn)
      dim(tsne)
      if(length(cn) != nrow(tsne)) {
        dbg("[pgxinfo.updateDatasetFolder] ***ERROR**** length mismatch!")
      }
      rownames(tsne) <- cn
      colnames(tsne) <- paste0("tsne.", 1:ncol(tsne))
      utils::write.csv(tsne, file = tsne.file)
    }
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
