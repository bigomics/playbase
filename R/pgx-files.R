##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' @export
pgx.load <- function(file, verbose = 0) {
  local(get(load(file, verbose = verbose)))
}

#' @export
ngs.save <- function(ngs, file, update.date = TRUE, light = TRUE, system = FALSE) {
  message("warning: ngs.save() is deprecated. please use pgx.save()")
  pgx.save(ngs, file = file, update.date = update.date, light = light, system = system)
}

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
  sort(sapply(pgx, object.size)) / 1e9
  sum(sapply(pgx, object.size)) / 1e9

  cat(">>> saving PGX file to", file, "\n")
  file <- iconv(file, from = "", to = "ASCII//TRANSLIT")
  save(pgx, file = file)
}


#' @export
h5exists <- function(h5.file, obj) {
  xobjs <- apply(rhdf5::h5ls(h5.file)[, 1:2], 1, paste, collapse = "/")
  obj %in% gsub("^/|^//", "", xobjs)
}


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

#' @export
pgx.readMatrixH5 <- function(h5.file, select = NULL, rows = NULL) {
  if (is.null(select) && is.null(rows)) {
    X <- rhdf5::h5read(h5.file, "data/matrix")
    rn <- rhdf5::h5read(h5.file, "data/rownames")
    cn <- rhdf5::h5read(h5.file, "data/colnames")
  }
  if (!is.null(select) || !is.null(rows)) {
    X <- rhdf5::h5read(h5.file, "data/matrix", index = list(rows, select))
    rn <- rhdf5::h5read(h5.file, "data/rownames", index = list(rows))
    cn <- rhdf5::h5read(h5.file, "data/colnames", index = list(select))
  }
  rownames(X) <- rn
  colnames(X) <- cn
  X[which(X < -999999)] <- NA
  as.matrix(X)
}

#' @export
pgx.readOptions <- function(file = "./OPTIONS") {
  if (!file.exists(file)) {
    return(NULL)
  }
  opt <- read.table(file, sep = "=", row.names = 1)
  opt <- gsub("^[ ]*|[ ]*$", "", apply(opt, 1, c)) ## strip leading/post spaces
  opt <- sapply(opt, list)
  opt <- sapply(opt, strsplit, split = "[;]")
  is.bool <- sapply(opt, function(x) all(tolower(x) %in% c("true", "false")))
  is.bool
  opt[is.bool] <- sapply(opt[is.bool], function(x) tolower(x) %in% c("true"))
  names(opt) <- trimws(names(opt))
  opt
}

#' @export
pgx.readDatasetProfiles <- function(pgx.dir, file = "datasets-allFC.csv", verbose = TRUE) {
  if (!dir.exists(pgx.dir)) {
    stop(paste("[readDatasetProfiles1] FATAL ERROR : folder", pgx.dir, "does not exist"))
  }
  fn <- file.path(pgx.dir, file)
  fn
  if (!file.exists(fn)) {
    stop("FATAL : could not find profiles matrix. please create first with initDatasetFolder().\n")

    return()
  } else {
    if (verbose) message("[readDatasetProfiles1] Found existing dataset profiles matrix")
  }

  allFC <- fread.csv(file = file.path(pgx.dir, file), row.names = 1, check.names = FALSE)
  allFC <- as.matrix(allFC)
  if (verbose) message("[readDatasetProfiles1] dataset profiles matrix : dim=", dim(allFC))
  return(allFC)
}

#' @export
pgx.scanInfoFile.DEPRECATED <- function(
    pgx.dir,
    file = "datasets-info.csv",
    allfc.file = "datasets-allFC.csv",
    info.file = "datasets-info.csv",
    force = FALSE,
    new.pgx = NULL,
    verbose = TRUE) {
  INITDATASETFOLDER <- TRUE

  pgx.files <- dir(pgx.dir, pattern = "[.]pgx$")
  if (length(pgx.files) == 0) {
    out <- list(INITDATASETFOLDER = FALSE)
    return(out) ## no files!
  }

  ## only run pgx.initDatasetFolder if pgx are changed

  if (!dir.exists(pgx.dir)) {
    stop(paste("[initDatasetFolder] FATAL ERROR : folder", pgx.dir, "does not exist"))
  }

  ## all public datasets
  pgx.dir <- pgx.dir[1] ## only one folder!!!
  pgx.files <- dir(pgx.dir, pattern = "[.]pgx$")
  pgx.files <- sub("[.]pgx$", "", pgx.files) ## strip pgx

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

  allfc.file1 <- file.path(pgx.dir, allfc.file)
  has.fc <- file.exists(allfc.file1)

  info.file1 <- file.path(pgx.dir, info.file)
  has.info <- file.exists(info.file1)

  ## ----------------------------------------------------------------------
  ## If an allFC file exits, check if it is done for all PGX files
  ## ----------------------------------------------------------------------

  pgxinfo <- NULL
  pgx.missing0 <- pgx.files
  pgx.missing1 <- pgx.files
  pgx.delete0 <- c()
  pgx.delete1 <- c()

  allFC <- NULL
  if (!force && has.fc) {
    if (verbose) message("[initDatasetFolder] checking which pgx files already done in allFC...")
    allFC <- data.table::fread(allfc.file1, check.names = FALSE, nrows = 1) ## HEADER!!!
    pgx.done <- gsub("^\\[|\\].*", "", colnames(allFC)[-1])
    pgx.missing0 <- setdiff(pgx.missing0, pgx.done)
    pgx.delete0 <- setdiff(pgx.done, pgx.files)
    allFC <- NULL
  }

  if (!force && has.info) {
    if (verbose) message("[initDatasetFolder] checking which pgx files already in PGX info...")
    ## do not use fread! quoting bug
    pgxinfo <- read.csv(info.file1, stringsAsFactors = FALSE, row.names = 1, sep = ",")
    pgxinfo.files <- unique(sub(".pgx$", "", pgxinfo$dataset))
    pgx.missing1 <- setdiff(pgx.missing1, pgxinfo.files)
    pgx.delete1 <- setdiff(pgxinfo.files, pgx.files)
  }

  ## ----------------------------------------------------------------------
  ## Check if it is done for all PGX files
  ## ----------------------------------------------------------------------

  ## files to be done either for allFC or missing in INFO
  pgx.missing <- unique(c(pgx.missing0, pgx.missing1))
  pgx.delete <- unique(c(pgx.delete0, pgx.delete1))

  pgxinfo.changed <- FALSE
  pgxfc.changed <- FALSE

  ## these pgx need forced update
  if (!is.null(new.pgx)) {
    new.pgx <- sub(".pgx$", "", new.pgx)
    new.pgx <- intersect(new.pgx, pgx.files) ## only existing pgx
    pgx.delete <- union(pgx.delete, new.pgx)
    pgx.missing <- union(pgx.missing, new.pgx)
    pgx.missing0 <- union(pgx.missing0, new.pgx)
    pgx.missing1 <- union(pgx.missing1, new.pgx)
    sel1 <- which(pgxinfo$dataset %in% new.pgx)
    if (length(sel1)) {
      pgxinfo <- pgxinfo[-sel1, ]
      pgxinfo.changed <- TRUE
    }
  }
  if (length(pgx.missing) == 0 && length(pgx.delete) == 0) {
    if (verbose) message("[initDatasetFolder] no update required. use FORCE=1 for forced update.")
    INITDATASETFOLDER <- FALSE
  }

  if (verbose) message("[initDatasetFolder] folder has ", length(pgx.missing), " new PGX files")
  if (verbose) message("[initDatasetFolder] info-file has ", length(pgx.delete), " old items")

  return(
    list(
      INITDATASETFOLDER = INITDATASETFOLDER,
      pgxinfo = pgxinfo,
      pgx.files = pgx.files,
      pgxinfo.changed = pgxinfo.changed,
      pgxfc.changed = pgxfc.changed,
      info.file1 = info.file1,
      pgx.missing = pgx.missing,
      pgx.missing0 = pgx.missing0,
      pgx.missing1 = pgx.missing1
    )
  )
}

#' @export
pgx.initDatasetFolder.DEPRECATED <- function(pgx.dir,
                                             allfc.file = "datasets-allFC.csv",
                                             info.file = "datasets-info.csv",
                                             info.file1 = NULL,
                                             force = FALSE,
                                             delete.old = FALSE,
                                             pgxinfo = NULL,
                                             pgx.files = NULL,
                                             pgxinfo.changed = NULL,
                                             pgxfc.changed = NULL,
                                             pgx.missing = NULL,
                                             pgx.missing1 = NULL,
                                             pgx.missing0 = NULL,
                                             new.pgx = NULL,
                                             update.sigdb = TRUE,
                                             verbose = TRUE) {
  ## ----------------------------------------------------------------------
  ## Reread allFC file. Before we only read the header.
  ## ----------------------------------------------------------------------

  allFC <- NULL
  allfc.file1 <- file.path(pgx.dir, allfc.file)
  if (!force && file.exists(allfc.file1) && length(pgx.missing) > 0) {
    allFC <- fread.csv(allfc.file1, row.names = 1, check.names = FALSE)
  }

  ## these pgx need forced update, so remove
  if (!is.null(allFC) && !is.null(new.pgx)) {
    new.pgx <- sub(".pgx$", "", new.pgx)
    allfc.pgx <- gsub("^\\[|\\].*", "", colnames(allFC))
    sel2 <- which(allfc.pgx %in% new.pgx)
    if (length(sel2)) {
      allFC <- allFC[, -sel2, drop = FALSE]
      pgxfc.changed <- TRUE
    }
  }

  ## ----------------------------------------------------------------------
  ## For all new PGX files, load the PGX file and get the meta FC
  ## matrix.
  ## ----------------------------------------------------------------------
  info.cols <- NULL
  missing.FC <- list()
  message("[initDatasetFolder] missing pgx = ", paste(pgx.missing, collapse = " "))
  pgxfile <- pgx.missing[1]

  ngs <- NULL
  for (pgxfile in pgx.missing) {
    cat(".")
    pgxfile1 <- file.path(pgx.dir, pgxfile)
    pgxfile1 <- paste0(sub("[.]pgx$", "", pgxfile1), ".pgx")

    pgx <- try(local(get(load(pgxfile1, verbose = 0)))) ## override any name

    if ("try-error" %in% class(pgx)) {
      message(paste("[initDatasetFolder] ERROR in loading PGX file:", pgxfile1, ". skipping\n"))
      next()
    }

    if (!pgx.checkObject(pgx)) {
      message(paste("[initDatasetFolder] INVALID PGX object", pgxfile, ". Skipping"))
      next()
    }

    ## ---------------------------------------------
    ## extract the meta FC matrix
    ## ---------------------------------------------

    if (pgxfile %in% pgx.missing0) {
      meta <- pgx.getMetaFoldChangeMatrix(pgx, what = "meta")
      rownames(meta$fc) <- toupper(rownames(meta$fc))
      missing.FC[[pgxfile]] <- meta$fc
      pgxfc.changed <- TRUE
    }

    ## ---------------------------------------------
    ## compile the info for update
    ## ---------------------------------------------
    if (pgxfile %in% pgx.missing1) {
      pgx$name <- sub(".pgx$", "", pgxfile) ## force filename as name
      pgxinfo <- pgx.updateInfoPGX(pgxinfo, pgx)
      pgxinfo.changed <- TRUE
    }
  }

  pgx <- NULL

  ## ----------------------------------------------------------------------
  ## Update the INFO meta file
  ## ----------------------------------------------------------------------
  rownames(pgxinfo) <- NULL
  pgxinfo <- data.frame(pgxinfo, check.names = FALSE)

  ## remove unneccessary entries if forced.
  sel.delete <- which(!sub(".pgx$", "", pgxinfo$dataset) %in% pgx.files)
  if (length(sel.delete) && delete.old) {
    pgxinfo <- pgxinfo[-sel.delete, , drop = FALSE]
    pgxinfo.changed <- TRUE
  }

  if (pgxinfo.changed) {
    if (verbose) message("[initDatasetFolder] writing updated PGX.INFO file to ", info.file1, "...")
    write.csv(pgxinfo, file = info.file1)
    Sys.chmod(info.file1, "0666")
  }

  ## ----------------------------------------------------------------------
  ## Update the ALL.FC meta file
  ## ----------------------------------------------------------------------

  ## remove unneccessary entries if forced.
  if (!is.null(allFC) && delete.old) {
    fc.done <- gsub("^\\[|\\].*", "", colnames(allFC))
    sel.deleteFC <- which(!fc.done %in% pgx.files)
    if (length(sel.deleteFC)) {
      allFC <- allFC[, -sel.deleteFC, drop = FALSE]
      pgxfc.changed <- TRUE
    }
  }

  if (length(missing.FC) == 0 && !pgxfc.changed) {
    ## no change in info
    dbg("[initDatasetFolder] allFC complete. no change needed.")
    return(NULL)
  }

  if (length(missing.FC) > 0) {
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
    allfc.sd <- apply(allFC, 1, sd, na.rm = TRUE)
    allfc.nna <- rowMeans(!is.na(allFC))
    jj <- Matrix::head(order(-allfc.sd * allfc.nna), 20000)
    allFC <- allFC[jj, , drop = FALSE]
    pgxfc.changed <- TRUE
  }

  ## save modified allFC
  if (pgxfc.changed) {
    ## check for duplicates
    if (verbose) message("[initDatasetFolder] allFC changed. updating file: ", allfc.file1)
    allFC <- allFC[, !duplicated(colnames(allFC)), drop = FALSE]
    allFC <- allFC[, order(colnames(allFC)), drop = FALSE]
    allFC <- round(allFC, digits = 4)
    AA <- data.frame(rownames = rownames(allFC), allFC, check.names = FALSE)

    data.table::fwrite(AA, file = allfc.file1)
    Sys.chmod(allfc.file1, "0666")
    remove(AA)
  }

  ## update user sigdb or create if not exists
  sigdb <- file.path(pgx.dir, "datasets-sigdb.h5")
  if (update.sigdb && (!file.exists(sigdb) || pgxfc.changed)) {
    ## NEED RETHINK!!!! HERE???
    if (file.exists(sigdb)) unlink(sigdb)
    if (verbose) message("[initDatasetFolder] creating signature DB to", sigdb, "...")
    pgx.createSignatureDatabaseH5.fromMatrix(sigdb, X = allFC)
    if (verbose) message("[initDatasetFolder] add enrichment signature to", sigdb, "...")
    pgx.addEnrichmentSignaturesH5(sigdb, X = allFC, methods = "rankcor")
  }

  # ## do not use fread.csv or fread here!! see issue #441

  pgxinfo <- playbase::pgxinfo.read(pgx.dir)
  return(pgxinfo)
}


#' Update PGX-table with new pgx object.
#'
#' @export
pgx.updateInfoPGX <- function(pgxinfo, pgx, remove.old = TRUE) {
  cond <- grep("title|source|group|batch|sample|patient|donor|repl|clone|cluster|lib.size|^[.]",
    colnames(pgx$samples),
    invert = TRUE, value = TRUE
  )
  cond <- grep("title|source|batch|sample|patient|donor|repl|clone|cluster|lib.size|^[.]",
    colnames(pgx$samples),
    invert = TRUE, value = TRUE
  )
  cond

  is.mouse <- (mean(grepl("[a-z]", pgx$genes$gene_name)) > 0.8)
  organism <- c("human", "mouse")[1 + is.mouse]
  if ("organism" %in% names(pgx)) organism <- pgx$organism

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

  pgxinfo
}

## ================================================================================
## PGXINFO methods (preferred API)
## ================================================================================

#' Delete pgx entry in datasets-info table in files (WIP)
#'
#' @export
pgxinfo.deletePgx <- function(pgx.dir, pgxname,
                              delete.fc = FALSE) {
  info.file <- file.path(pgx.dir, "datasets-info.csv")
  pgxname <- sub("[.]pgx$", "", pgxname)
  pgxinfo <- read.csv(info.file, row.names = 1)
  pgxinfo

  info_datasets <- sub("[.]pgx$", "", pgxinfo$dataset)
  sel <- which(info_datasets == pgxname)
  sel
  if (!is.null(sel) && length(sel) > 0) {
    pgxinfo <- pgxinfo[-sel, , drop = FALSE]
    write.csv(pgxinfo, file = info.file)
  }

  ## Should we also delete the entry in allFC and sigdb? This will
  ## take some overhead, so if not needed better skip.
  allfc.file <- file.path(pgx.dir, "datasets-allFC.csv")
  tsne.file <- file.path(pgx.dir, "datasets-tsne.csv")
  h5.file <- file.path(pgx.dir, "datasets-sigdb.h5")

  ## delete columns from allFC file
  if (delete.fc && file.exists(allfc.file)) {
    allFC <- data.table::fread(allfc.file, check.names = FALSE) ## HEADER!!!
    allFC <- as.matrix(allFC, rownames = 1)
    del <- grep(paste0("\\[", pgxname, "\\]"), colnames(allFC))
    del
    if (length(del) > 0) {
      allFC <- allFC[, -del, drop = FALSE]
      allFC <- round(allFC, digits = 4)
      allFC1 <- data.frame(gene = rownames(allFC), allFC, check.names = FALSE)
      data.table::fwrite(allFC1, allfc.file) ## HEADER!!!
    }
  }

  ## delete rows from t-SNE file
  if (delete.fc && file.exists(tsne.file)) {
    tsne <- data.table::fread(tsne.file, check.names = FALSE) ## HEADER!!!
    tsne <- as.matrix(tsne, rownames = 1)
    del <- grep(paste0("\\[", pgxname, "\\]"), rownames(tsne))
    if (length(del)) {
      tsne <- tsne[-del, , drop = FALSE]
      write.csv(tsne, file = tsne.file)
    }
  }

  ## delete dataset from H5 file
  if (delete.fc && file.exists(h5.file)) {
    sigdb.removeDataset(h5.file, pgxname)
  }
}


#' Update PGX-info file with new pgx object (WIP)
#'
#' @export
pgxinfo.addPgx <- function(pgx.dir, pgx, file = "datasets-info.csv",
                           update.fc = TRUE) {
  info.file <- file.path(pgx.dir, file)
  pgxinfo <- read.csv(info.file, row.names = 1)
  pgxinfo <- pgx.updateInfoPGX(pgxinfo, pgx, remove.old = TRUE)
  write.csv(pgxinfo, file = info.file)

  allfc.file <- file.path(pgx.dir, "datasets-allFC.csv")
  tsne.file <- file.path(pgx.dir, "datasets-tsne.csv")
  h5.file <- file.path(pgx.dir, "datasets-sigdb.h5")
  pgxname <- sub("[.]pgx$", "", pgx$name)
  pgxname

  ## add FC columns to allFC file
  if (update.fc && !file.exists(allfc.file)) {
    ## complete update of all files
    pgx.initDatasetFolder(pgx.dir)
    return()
  }

  ## add FC columns to allFC file
  allfC <- NULL
  if (update.fc && file.exists(allfc.file)) {
    allFC <- data.table::fread(allfc.file, check.names = FALSE) ## HEADER!!!
    allFC <- as.matrix(allFC, rownames = 1)
    del <- grep(paste0("\\[", pgxname, "\\]"), colnames(allFC))
    if (length(del) > 0) {
      allFC <- allFC[, -del, drop = FALSE]
    }
    F1 <- playbase::pgx.getMetaMatrix(pgx)$fc
    F1 <- F1[match(rownames(allFC), rownames(F1)), ]
    colnames(F1) <- paste0("[", pgxname, "] ", colnames(F1))
    allFC <- cbind(allFC, F1)
    allFC <- round(allFC, digits = 4)
    allFC1 <- data.frame(gene = rownames(allFC), allFC, check.names = FALSE)
    data.table::fwrite(allFC1, allfc.file) ## HEADER!!!

    ## NEED RETHINK!!!! This could be done perhaps more efficient
    ## when updating with one extra dataset.
    pgx.createSignatureDatabaseH5.fromMatrix(h5.file,
      X = allFC,
      update.only = FALSE
    )

    ## update t-SNE file (in the future we do not need this file)
    ## Just get the tSNE from sigdb
    cn <- rhdf5::h5read(h5.file, "data/colnames")
    tsne <- rhdf5::h5read(h5.file, "clustering/tsne2d")
    rownames(tsne) <- cn
    colnames(tsne) <- paste0("tsne.", 1:ncol(tsne))
    write.csv(tsne, file = tsne.file)
  }
}

#' Read dataset information file (aka PGX info)
#'
#' @export
pgxinfo.read <- function(pgx.dir, file = "datasets-info.csv", match = TRUE) {
  pgx.files <- dir(pgx.dir, pattern = "[.]pgx$")
  if (length(pgx.files) == 0) {
    return(NULL)
  } ## no files!

  pgxinfo <- NULL
  pgxinfo.file <- file.path(pgx.dir, file)
  if (file.exists(pgxinfo.file)) {
    ## do not use fread.csv or fread here!! see issue #441
    pgxinfo <- read.csv(pgxinfo.file, stringsAsFactors = FALSE, row.names = 1, sep = ",")
    if (match) {
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


#' Checks if any of metadata files (info, allFC, sigdb) need
#' update. All pgx files in the datafolder should be included.
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

  fc.complete <- info.complete <- h5.complete <- TRUE
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
    pgxinfo <- read.csv(info.file1, stringsAsFactors = FALSE, row.names = 1, sep = ",")
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
    if (h5.ok) {
      cn <- rhdf5::h5read(sigdb.file1, "data/colnames")
      h5.files <- gsub("^\\[|\\].*", "", cn)
      h5.files <- sub("[.]pgx$", "", h5.files) ## strip pgx
      h5.complete <- all(pgx.files %in% h5.files)
      h5.missing <- setdiff(pgx.files, h5.files)
    }
  }

  if (verbose) {
    dbg("[pgxinfo.needUpdate] fc.complete = ", fc.complete)
    dbg("[pgxinfo.needUpdate] info.complete =  ", info.complete)
    if (check.sigdb) dbg("[pgxinfo.needUpdate] h5.complete =  ", h5.complete)

    dbg("[pgxinfo.needUpdate] nr missing files in allFC : ", length(fc.missing))
    dbg("[pgxinfo.needUpdate] nr missing files in info  : ", length(info.missing))
    if (check.sigdb) dbg("[pgxinfo.needUpdate] nr missing files in sigdb : ", length(h5.missing))
  }

  ## Return checks
  has.files <- (has.fc && has.info)
  if (check.sigdb) has.files <- has.files && has.sigdb
  is.complete <- (fc.complete && info.complete)
  if (check.sigdb) is.complete <- is.complete && h5.complete
  return(!has.files || !is.complete)
}



#' Update
#'
#' @export
pgxinfo.updateDatasetFolder <- function(pgx.dir,
                                        force = FALSE,
                                        delete.old = FALSE,
                                        new.pgx = NULL,
                                        update.sigdb = TRUE,
                                        verbose = TRUE) {
  ## only run pgx.initDatasetFolder if pgx are changed
  if (!dir.exists(pgx.dir)) {
    stop(paste("[pgxinfo.updateDatasetFolder] FATAL ERROR : folder", pgx.dir, "does not exist"))
  }

  pgx.files <- dir(pgx.dir, pattern = "[.]pgx$")
  if (length(pgx.files) == 0) {
    allfc.file1 <- file.path(pgx.dir, allfc.file)
    info.file1 <- file.path(pgx.dir, info.file)
    file.remove(info.file1)
    file.remove(allfc.file1)
    # should return FALSE?
    return(NULL)
  }

  ## all public datasets
  pgx.dir <- pgx.dir[1] ## only one folder!!!
  pgx.files <- dir(pgx.dir, pattern = "[.]pgx$")
  pgx.files <- sub("[.]pgx$", "", pgx.files) ## strip pgx

  ## ----------------------------------------------------------------------
  ## If an allFC file exists
  ## ----------------------------------------------------------------------

  allfc.file <- file.path(pgx.dir, "datasets-allFC.csv")
  has.fc <- file.exists(allfc.file)

  info.file <- file.path(pgx.dir, "datasets-info.csv")
  has.info <- file.exists(info.file)

  sigdb.file <- file.path(pgx.dir, "datasets-sigdb.h5")
  has.sigdb <- file.exists(sigdb.file)

  tsne.file <- file.path(pgx.dir, "datasets-tsne.csv")
  has.tsne <- file.exists(tsne.file)

  if (force) {
    has.fc <- FALSE
    has.info <- FALSE
    has.sigdb <- FALSE
    has.tsne <- FALSE
  }

  ## ----------------------------------------------------------------------
  ## If an allFC file exits, check if it is done for all PGX files
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
    if (verbose) message("[pgxinfo.updateDatasetFolder] checking which pgx files already done in allFC...")
    allFC <- data.table::fread(allfc.file, check.names = FALSE, nrows = 1) ## HEADER!!!
    fc.files <- gsub("^\\[|\\].*", "", colnames(allFC)[-1])
    fc.missing <- setdiff(pgx.files, fc.files)
    fc.delete <- setdiff(fc.files, pgx.files)
    allFC <- NULL
  }

  if (has.info) {
    if (verbose) message("[pgxinfo.updateDatasetFolder] checking which pgx files already in PGX info...")
    ## do not use fread! quoting bug
    pgxinfo <- read.csv(info.file, stringsAsFactors = FALSE, row.names = 1, sep = ",")
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
    if (verbose) message("[pgxinfo.updateDatasetFolder] checking which pgx files already in sigdb...")
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
    sel1 <- which(pgxinfo$dataset %in% new.pgx)
    if (length(sel1)) {
      pgxinfo <- pgxinfo[-sel1, ]
      pgxinfo.changed <- TRUE
    }
  }

  if (verbose) message("[pgxinfo.updateDatasetFolder] length(pgx.missing) = ", length(pgx.missing), "")
  if (verbose) message("[pgxinfo.updateDatasetFolder] length(pgx.delete) =  ", length(pgx.delete), "")

  ## nothing to do???
  if (length(pgx.missing) == 0 && length(h5.missing) == 0) {
    if (verbose) message("[pgxinfo.updateDatasetFolder] no update required. Use FORCE for forced update.")
    return()
  }

  ## pgxinfo and allFC OK, but sigdb not upto date
  if (update.sigdb && length(pgx.missing) == 0 && length(h5.missing) > 0) {
    ## Update SIGDB file using allFC
    allFC <- fread.csv(allfc.file, row.names = 1, check.names = FALSE)
    allFC <- as.matrix(allFC)
    ## update user sigdb or create if not exists
    if (file.exists(sigdb.file)) unlink(sigdb.file)
    if (verbose) message("[updateDatasetFolder] creating sig DB: ", sigdb.file, "...")
    pgx.createSignatureDatabaseH5.fromMatrix(sigdb.file, X = allFC)
    pgx.addEnrichmentSignaturesH5(sigdb.file, X = allFC, methods = "rankcor")
    h5.missing <- NULL
  }

  if (length(pgx.missing) == 0 && length(h5.missing) == 0) {
    if (verbose) message("[updateDatasetFolder] no update required. Use FORCE for forced update.")
    return()
  }

  ## ----------------------------------------------------------------------
  ## Reread allFC file. Before we only read the header.
  ## ----------------------------------------------------------------------
  allFC <- NULL

  if (!force && file.exists(allfc.file) && length(pgx.missing) > 0) {
    allFC <- fread.csv(allfc.file, row.names = 1, check.names = FALSE)
    allFC <- as.matrix(allFC)
  }

  ## these pgx need forced update, so remove
  if (!is.null(allFC) && !is.null(new.pgx)) {
    new.pgx <- sub(".pgx$", "", new.pgx)
    allfc.pgx <- gsub("^\\[|\\].*", "", colnames(allFC))
    sel2 <- which(allfc.pgx %in% new.pgx)
    if (length(sel2)) {
      allFC <- allFC[, -sel2, drop = FALSE]
      pgxfc.changed <- TRUE
    }
  }

  ## ----------------------------------------------------------------------
  ## For all new PGX files, load the PGX file and get the meta FC
  ## matrix.
  ## ----------------------------------------------------------------------
  info.cols <- NULL
  missing.FC <- list()
  message("[updateDatasetFolder] missing pgx = ", paste(pgx.missing, collapse = " "))
  pgxfile <- pgx.missing[1]

  ngs <- NULL
  for (pgxfile in pgx.missing) {
    cat(".")
    pgxfile1 <- file.path(pgx.dir, pgxfile)
    pgxfile1 <- paste0(sub("[.]pgx$", "", pgxfile1), ".pgx")

    pgx <- try(local(get(load(pgxfile1, verbose = 0)))) ## override any name

    if ("try-error" %in% class(pgx)) {
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
      pgxinfo <- pgx.updateInfoPGX(pgxinfo, pgx)
      pgxinfo.changed <- TRUE
    }
  }

  pgx <- NULL

  ## ----------------------------------------------------------------------
  ## Update the INFO meta file
  ## ----------------------------------------------------------------------
  rownames(pgxinfo) <- NULL
  pgxinfo <- data.frame(pgxinfo, check.names = FALSE)

  ## remove unneccessary entries if forced.
  sel.delete <- which(!sub(".pgx$", "", pgxinfo$dataset) %in% pgx.files)
  if (length(sel.delete) && delete.old) {
    pgxinfo <- pgxinfo[-sel.delete, , drop = FALSE]
    pgxinfo.changed <- TRUE
  }

  if (pgxinfo.changed) {
    if (verbose) message("[pgxinfo.updateDatasetFolder] writing updated PGX.INFO file to ", info.file, "...")
    write.csv(pgxinfo, file = info.file)
    Sys.chmod(info.file, "0666")
  }

  ## ----------------------------------------------------------------------
  ## Update the ALL.FC meta file
  ## ----------------------------------------------------------------------

  ## remove unneccessary entries if forced.
  if (!is.null(allFC) && delete.old) {
    fc.files <- gsub("^\\[|\\].*", "", colnames(allFC))
    sel.deleteFC <- which(!fc.files %in% pgx.files)
    if (length(sel.deleteFC)) {
      allFC <- allFC[, -sel.deleteFC, drop = FALSE]
      pgxfc.changed <- TRUE
    }
  }

  if (length(missing.FC) == 0 && !pgxfc.changed) {
    ## no change in info
    dbg("[pgxinfo.updateDatasetFolder] allFC complete. no change needed.")
    return(NULL)
  }

  if (length(missing.FC) > 0) {
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
    allfc.sd <- apply(allFC, 1, sd, na.rm = TRUE)
    allfc.nna <- rowMeans(!is.na(allFC))
    jj <- Matrix::head(order(-allfc.sd * allfc.nna), 20000)
    allFC <- allFC[jj, , drop = FALSE]
    pgxfc.changed <- TRUE
  }

  ## save modified allFC
  if (pgxfc.changed) {
    ## check for duplicates
    if (verbose) message("[pgxinfo.updateDatasetFolder] allFC changed. updating file: ", allfc.file)
    allFC <- allFC[, !duplicated(colnames(allFC)), drop = FALSE]
    allFC <- allFC[, order(colnames(allFC)), drop = FALSE]
    allFC <- round(allFC, digits = 4)
    AA <- data.frame(rownames = rownames(allFC), allFC, check.names = FALSE)

    data.table::fwrite(AA, file = allfc.file)
    Sys.chmod(allfc.file, "0666")
    remove(AA)
  }

  ## update user sigdb or create if not exists
  sigdb <- file.path(pgx.dir, "datasets-sigdb.h5")
  if (update.sigdb) {
    if (!file.exists(sigdb) || pgxfc.changed) {
      ## NEED RETHINK!!!! HERE???
      if (file.exists(sigdb)) unlink(sigdb)
      if (verbose) message("[pgxinfo.updateDatasetFolder] creating signature DB to", sigdb, "...")
      pgx.createSignatureDatabaseH5.fromMatrix(sigdb, X = allFC)
      if (verbose) message("[pgxinfo.updateDatasetFolder] add enrichment signature to", sigdb, "...")
      pgx.addEnrichmentSignaturesH5(sigdb, X = allFC, methods = "rankcor")
    }

    tsne.file <- file.path(pgx.dir, "datasets-tsne.csv")
    if (!file.exists(tsne.file) || pgxfc.changed) {
      rhdf5::h5ls(sigdb)
      cn <- rhdf5::h5read(sigdb, "data/colnames")
      tsne <- rhdf5::h5read(sigdb, "clustering/tsne2d")
      rownames(tsne) <- cn
      colnames(tsne) <- paste0("tsne.", 1:2)
      write.csv(tsne, file = tsne.file)
    }
  }

  return()
}
