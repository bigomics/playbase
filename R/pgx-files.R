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
  pgx <- local(get(load(file, verbose = verbose)))
  pgx$filename <- file
  pgx
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
  if (is(pgx, "reactivevalues")) pgx <- shiny::reactiveValuesToList(pgx)
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
  pgx$filename <- NULL

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
  if (is.null(rownames(X))) rownames(X) <- 1:nrow(X)
  if (is.null(colnames(X))) colnames(X) <- paste0("V", 1:ncol(X))
  rhdf5::h5write(rownames(X), h5.file, "data/rownames")
  rhdf5::h5write(colnames(X), h5.file, "data/colnames")

  rhdf5::h5closeAll()
}


#' @title Read (sub-)matrix from HDF5 file
#'
#' @export
h5.readMatrix <- function(h5.file, rows = NULL, cols = NULL,
                          matrixid='data/matrix', rowid='data/rownames',
                          colid='data/colnames' ) {
  cn <- rhdf5::h5read(h5.file, colid)
  rn <- rhdf5::h5read(h5.file, rowid)
  rowidx <- 1:length(rn)
  colidx <- 1:length(cn)
  if (!is.null(rows)) {
    if(all(is.integer(rows))) rows <- rn[rows]
    rowidx <- match(intersect(rows, rn), rn)
  }
  if (!is.null(cols)) {
    if(all(is.integer(cols))) cols <- cn[cols]
    colidx <- match(intersect(cols, cn), cn)
  }
  nr <- length(rowidx)
  nc <- length(colidx)
  message("[sigdb.getConnectivityMatrix] reading large H5 file: ", nr, "x", nc, "")
  X <- rhdf5::h5read(h5.file, matrixid, index = list(rowidx, colidx))
  rownames(X) <- rn[rowidx]
  colnames(X) <- cn[colidx]
  return(X)
}

#' @export
h5.rownames <- function(h5.file, id = NULL) {
  if(is.null(id)) {
    tags <- sub('^/','',apply(rhdf5::h5ls(sigdb)[,1:2],1,paste,collapse='/'))
    ids <- c('data/rownames','meta/genes')
    id <- intersect(ids,tags)
  }
  id  
  if(is.null(id)) {
    warning("could not resolve rownames")
    return(NULL)
  }
  rhdf5::h5read(h5.file, id)
}

#' @export
h5.colnames <- function(h5.file, id=NULL) {
  if(is.null(id)) {
    tags <- sub('^/','',apply(rhdf5::h5ls(sigdb)[,1:2],1,paste,collapse='/'))
    ids <- c('data/colnames','meta/sampleid','meta/sample')
    id <- intersect(ids,tags)
  }
  id
  if(is.null(id)) {
    warning("could not resolve colnames")
    return(NULL)
  }
  rhdf5::h5read(h5.file, id)
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
pgx.readOptions <- function(file = "./OPTIONS", default = NULL) {
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
  if (!is.null(default)) {
    ## add default options not in options
    ## which.add <- c("a","b","c")
    which.add <- setdiff(names(default), names(opt))
    if (length(which.add)) {
      cat("[pgx.readOptions] adding missing defaults:", paste(which.add, collapse = " "), "\n")
      opt <- c(opt, default[which.add])
    }
    which.extra <- setdiff(names(opt), names(default))
    if (length(which.extra)) {
      cat("[pgx.readOptions] warning: extra options:", paste(which.extra, collapse = " "), "\n")
    }
  }
  opt <- opt[order(names(opt))]
  opt
}


#' @title Cache a CSV file locally or on S3
#'
#' @param file Character string specifying the path to the CSV file to cache
#' @param bucket Character string specifying the S3 bucket name if caching on S3. Default is NULL for local caching.
#' @param FUN Function to use for reading the CSV file. Default is read.csv.
#' @param ... Other arguments passed to FUN.
#'
#' @return A data frame containing the contents of the cached CSV file.
#'
#' @details This function checks if a cached copy of the CSV file exists locally in /tmp or on S3.
#' If not, it reads the file using FUN(), caches it locally or on S3, and returns the contents.
#' On subsequent calls it returns the cached copy instead of re-reading the file.
#'
#' @examples
#' \dontrun{
#' # Cache locally
#' df <- cached.csv("data.csv")
#'
#' # Cache on S3
#' df <- cached.csv("s3://mybucket/data.csv", bucket = "mybucket")
#' }
#' @export
cached.csv <- function(file, FUN = read.csv, force = FALSE, ...) {
  file2 <- sub(getwd(), "", file)
  file2 <- paste0("cache-", gsub("[-._/]", "", file2), ".rds")
  file2
  ## cache.file <- file.path("/tmp", file2)
  cache.dir <- ifelse(dir.exists("cache"), "./cache", tempdir())
  cache.file <- file.path(cache.dir, file2)
  if (force || !file.exists(cache.file)) {
    if (!file.exists(file)) {
      return(NULL)
    }
    message("[cached.csv] reading from file ", file)
    csv <- FUN(file, ...)
    message("[cached.csv] saving cache at ", cache.file)
    saveRDS(csv, file = cache.file)
  } else {
    message("[cached.csv] reading from cache ", cache.file)
    csv <- readRDS(cache.file)
  }

  return(csv)
}

#' @export
cached.csv.s3 <- function(file, bucket, FUN = read.csv, force = FALSE, ...) {
  file2 <- sub(getwd(), "", file)
  file2 <- paste0("cache-", gsub("[-._/]", "", file2), ".rds")
  file2
  cache.dir <- ifelse(dir.exists("cache"), "./cache", tempdir())
  cache.file <- file.path(cache.dir, file2)
  if (force || !file.exists(cache.file)) {
    if (!file.exists(file)) {
      return(NULL)
    }
    message("[cached.csv.s3] reading from S3 bucket=", bucket, "   file=", file)
    obj <- aws.s3::get_object(object = file, bucket = bucket)
    csv <- data.table::fread(rawToChar(obj))
    message("[cached.csv] saving cache at ", cache.file)
    saveRDS(csv, file = cache.file)
  } else {
    message("[cached.csv.s3] reading from cache ", cache.file)
    csv <- readRDS(cache.file)
  }
  csv
}

#' @export
sampleMatrixFromNames <- function(names) {
  samples <- do.call(rbind, strsplit(names, split = "[_]"))

  ## remove unique columns (probably sample ID)
  id.cols <- which(apply(samples, 2, function(x) !any(duplicated(x))))
  group.cols <- which(apply(samples, 2, function(x) any(duplicated(x))))

  ## give rownames and columnames
  rownames(samples) <- names
  colnames(samples) <- paste0("V", 1:ncol(samples))
  colnames(samples)[group.cols] <- paste0("group", 1:length(group.cols))
  if (length(group.cols) == 1) colnames(samples)[group.cols[1]] <- "group"
  colnames(samples)[id.cols] <- paste0("id", 1:length(id.cols))
  if (length(id.cols) == 1) colnames(samples)[id.cols[1]] <- "id"

  ## samples <- type.convert(data.frame(samples), as.is=TRUE)
  samples
}
