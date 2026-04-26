mem.vmrss <- function(digits = 0) {
  mem <- "[? MB]"
  if (Sys.info()["sysname"] %in% c("Linux")) {
    proc <- paste("/proc", Sys.getpid(), "status", sep = "/")
    rss <- gsub("VmRSS:[\t ]+| kB", "", system(paste("grep -i vmrss", proc), intern = TRUE))
    rss <- as.numeric(rss) / (1024) ## MB
    mem <- paste0(round(rss, digits), "MB")
  }
  mem
}

mem.proc <- function(digits = 0) {
  mem <- "[? MB]"
  if (Sys.info()["sysname"] %in% c("Linux")) {
    file <- paste("/proc", Sys.getpid(), "stat", sep = "/")
    what <- vector("character", 52)
    ## In your logging routine
    vsz <- as.numeric(scan(file, what = what, quiet = TRUE)[23])
    vsz <- vsz / (1024**2) ## MB
    mem <- paste0(round(vsz, digits), "MB")
  }
  mem
}

info <- function(..., type = "INFO") {
  dd <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg <- "some message"
  msg <- sapply(list(...), paste, collapse = " ")
  dd <- paste0("[", dd, "]")
  mm <- paste0("[", mem.proc(), "/", mem.vmrss(), "]")
  type <- paste0("[", type, "]")
  message(paste0(type, dd, mm, " --- ", sub("\n$", "", paste(msg, collapse = " "))))
}

dbg <- function(...) info(..., type = "DBUG")

#' @export
capitalize <- function(s) {
  paste(toupper(substring(s, 1, 1)), tolower(substring(s, 2, nchar(s))), sep = "")
}

#' Get Call Stack
#'
#' Returns the call stack up to a specified depth.
#'
#' @param n Number of calls to include in call stack.
#' @param n0 Number of calls to exclude from the call stack.
#' @return Character vector of call stack.
#'
#' @examples
#' \dontrun{
#' f3 <- function(x) {
#'   print("Show whole stack")
#'   message(get_call_stack())
#'   print("Remove message call")
#'   message(get_call_stack(1))
#'   x <- x + 1
#'   x
#' }
#'
#' f2 <- function(x) {
#'   x <- f3(x)
#'   x
#' }
#'
#' f1 <- function(x) {
#'   x
#'   x <- f2(x)
#' }
#' f0 <- function(x) {
#'   x
#'   x <- f1(x)
#' }
#' f0(2)
#' }
get_call_stack <- function(n = 3, n0 = 2) {
  # Get calls, extract last n remove  and remove the last n calls
  calls <- sys.calls()
  fun_names <- lapply(calls, function(x) {
    x <- as.character(x[[1]])
    if (length(x) > 1) {
      x <- x[length(x)]
    }
    return(x)
  })
  n <- ifelse(n > length(fun_names), length(fun_names) - 1, n)
  exclude_last_n <- length(fun_names) - n
  n0 <- ifelse(n0 > exclude_last_n, exclude_last_n - 1, n0)
  get_last_n0 <- exclude_last_n - n0
  fun_names <- fun_names[get_last_n0:exclude_last_n]
  # Concatenate all the calls using the ">" symbol
  call_stack <- paste(fun_names, collapse = ">")

  return(call_stack)
}


#' @title Print Info Message
#'
#' @description Print a formatted info message. By default it includes the last 3 calls.
#'
#' @param msg Message to print.
#' @param n Number of calls to include in call stack.
#' @param n0 Number of calls to exclude from the call stack.
#'
#' @return Prints a formatted message.
#'
#' @examples
#' \dontrun{
#' info_message("Hello")
#' f <- function(x) {
#'   playbase::infomessage("test")
#'   x <- x + 1
#'   x
#' }
#' f1 <- function(x) {
#'   playbase::infomessage("test1")
#'   x <- f(x)
#' }
#' f2 <- function(x) {
#'   playbase::infomessage("test2")
#'   x <- f1(x)
#' }
#'
#' f3 <- function(x) {
#'   playbase::infomessage("test3")
#'   x <- f2(x)
#' }
#'
#' f4 <- function(x) {
#'   playbase::infomessage("test3")
#'   x <- f3(x)
#' }
#'
#' f(1)
#' f4(1)
#' }
#' @export
infomessage <- function(..., n = 3, n0 = 2) {
  stack <- paste0("[", get_call_stack(n, n0), "]: ")
  message(stack, ...)
}


#' Get mini example dataset
#'
#' @description
#' Returns a mini example dataset with sample information,
#' gene expression counts, and contrasts.
#'
#' @param item Select which elements to return c("counts", "samples", "contrats", "all")
#' by default it returns a list with all the data elements to build a pgx object.
#' @return List containing:
#' \itemize{
#'   \item samples - data.frame with sample info
#'   \item counts - matrix of gene counts
#'   \item contrast - data.frame with contrasts
#' }
#'
#' @examples
#' d <- get_mini_example_data()
#' @export
get_mini_example_data <- function(item = "all") {
  if (!item %in% c("counts", "samples", "contrats", "all")) {
    stop("Item not matched. Please use: 'counts', 'samples', 'contrats', or 'all'")
  }
  if (length(item) > 1L) {
    stop("Call only one item at a time. For multiple items call, use 'all'")
  }

  # Prepare input data
  counts <- playbase::COUNTS
  samples <- playbase::SAMPLES
  contrast <- playbase::CONTRASTS

  n_genes <- round(seq(1, nrow(counts), length.out = 500))

  # Subset each data frame to facilitate testing
  mini_counts <- counts[n_genes, c(1:3, 8:9, 11:12)]
  mini_samples <- samples[colnames(mini_counts), ]
  mini_contrast <- contrast[1:3, 1:3]
  mini_data <- list(
    counts = mini_counts,
    samples = mini_samples,
    contrast = mini_contrast
  )

  if (item == "all") {
    return(mini_data)
  } else {
    return(mini_data[[item]])
  }
}


#' Get row names for a data.table
#'
#' @description This function aims to extract the rownames of a data.table object as
#' if it would be a data.frame. The implementation is slow but it helps to maintain back
#' compatibility with data.frame objects.
#'
#' @param dt A data.table object.
#'
#' @return Character vector of row names
#'
#' @examples
#' \dontrun{
#' dt <- data.table(rn = letters[1:5], x = 1:5)
#' rownames(dt)
#' }
#' @export
rownames.data.table <- function(dt) {
  if (is.character(dt[[1]])) {
    return(dt[[1]])
  } else if ("rn" %in% names(dt)) {
    return(dt$rn)
  } else {
    warning("No column contains row names")
    return(NULL)
  }
}

#' Double Center and Scale a Matrix
#'
#' This function performs double centering and scaling on a matrix by:
#' 1. Subtracting row means
#' 2. Subtracting column means
#' 3. Scaling by column standard deviations
#'
#' @param X A numeric matrix to be centered and scaled
#' @return A transposed matrix that has been double centered and scaled
#' @export
double_center_scale_fast <- function(X) {
  X <- X - matrixStats::rowMeans2(X)
  X <- sweep(X, 2, matrixStats::colMeans2(X), "-")
  X <- sweep(X, 2, matrixStats::colSds(X), "/")
  return(t(X))
}

#' Compute coefficient of variation (CV, %)
#' This function computes coefficient of variation (CV) of each row (feature) of data matrix
#' It does so in the linear space.
#' @param counts A numeric matrix
#' @return CV (%) for each row (feature) of the input matrix.
#' @export
compute_CV <- function(counts) {
  if (playbase::is_logged(counts)) counts <- pmax(2**counts - 1, 0)
  sdx <- matrixStats::rowSds(counts, na.rm = TRUE)
  avg <- matrixStats::rowMeans2(counts, na.rm = TRUE)
  cv <- (sdx / avg) * 100
  return(cv)
}

#' Sum up tech replicates of a linear count matrix 
#' This function combines (sums up) features of technical replicates
#' It does so in the linear space.
#' @param counts Numeric matrix
#' @param trep_var Character vector indicating tech rep info for each sample
#' @return counts matrix with technical replicates combined
#' @export
sum_treps <- function(counts, trep_var = "") {

  if (is.null(trep_var)) {
    return(counts)
  } else if (length(trep_var) != ncol(counts)) {
    message("[playbase::sum_treps] trep_var must have same length as ncol of counts")
    return(counts)
  } else {
    if (playbase::is_logged(counts)) {
      prior <- min(counts[which(counts > 0)], na.rm = TRUE)
      counts <- pmax(2**counts - prior, 0)
    }
    trep_var <- as.character(trep_var)
    names(trep_var) <- colnames(counts)
    treps <- unique(trep_var)
    for(i in 1:length(treps)) {
      jj <- which(trep_var == treps[i])
      if (length(jj) <= 1) next
      kk <- names(jj)
      counts1 <- counts[, kk, drop = FALSE]
      mat <- rowSums(counts1, na.rm = TRUE)
      mat[rowSums(is.na(counts1)) == length(kk)] <- NA
      counts[, kk[1]] <- mat
      jx  <- match(kk, colnames(counts))
      counts <- counts[, -jx[-1], drop = FALSE]
    }
  }
  
  return(counts)
  
}

#' Re-order uniprot IDs (when a feature has multiple uniprot ids)
#' Rules:
#' If there are Swiss-Prot entries, Swiss-Prot entries come first, followed by TrEMBL entries.
#' Among Swiss-Prot entries, if there are canonical entries, canonical entries come first (eg P05067),
#' followed by isoforms (P05067-2). In absence of canonical entries, convert isoforms to canonical
#' (eg P05067-2 to P05067).
#' Among Swiss-Prot entries (after step 2), longer proteins come first.
#' Among TrEMBL entires, longer proteins come first.
#' @param feature feature ids, collapsed with ";".
#' @param lengths Character vector of feature lengths, collapsed with ";". Default NULL.
#' @export
rank_uniprots <- function(feature, lengths = NULL, verbose = FALSE) {

  ff <- strsplit(as.character(feature), ";")[[1]]
  ff <- ff[nzchar(trimws(ff))]

  if (all(as.character(feature) %in% c("", "NA", NA))) {
    return(list(feature=feature, lengths=lengths))
  }
  
  if (length(ff) <= 1) {
    sp <- grep("^[OPQ][0-9]", ff)
    if (length(sp) > 0) {
      hh <- grep("[-_][0-9]+$", ff)
      if (length(hh) > 0) ff <- sub("[-_][0-9]+$", "", ff)
    }
    return(list(feature=ff, lengths=lengths))
  }

  ll <- NULL
  if (!is.null(lengths)) {
    lengths <- strsplit(as.character(lengths), ";")[[1]]
    lengths <- lengths[nzchar(trimws(lengths))]
    if (length(ff) != length(lengths)) {
      message("[playbase::rank_uniprots] Protein lengths does not match features. Ignoring length.")
    } else {
      names(lengths) <- ff
      kk <- which(as.character(lengths) %in% c("0", "", "NA", NA))
      if (length(kk) > 0) {
        if (length(kk) == length(lengths)) {
          message("[playbase::rank_uniprots] All protein lengths uninformative. Ignoring length.")
        } else {
          ll <- lengths
          ll[kk] <- "0"
        }
      } else {
        ll <- lengths
      }
    }
  }
  
  ## Swiss-Prot entries come first (if any), followed by TrEMBL entries (if any).
  ff.sp=NULL; ff.trembl=NULL
  sp <- grep("^[OPQ][0-9]", ff)
  if (length(sp) > 0) {
    ## Swiss-Prot canonical entries come first (P05067), followed by isoforms (eg., P05067-2).
    ## In absence of canonical entries, convert isoforms to canonical (eg., P05067-2 to P05067).
    ff.sp <- ff[sp]
    hh <- grep("[-_][0-9]+$", ff.sp)
    if (length(hh) > 0) {      
      if (length(hh) == length(ff.sp)) {
        ff.sp[hh] <- sub("[-_][0-9]+$", "", ff.sp[hh])
        if (!is.null(ll)) {
          names(ll)[sp][hh] <- ff.sp[hh]
          names(lengths)[sp][hh] <- ff.sp[hh]
        }
      } else  {
        canon <- sub("[-_][0-9]+$", "", ff.sp[hh])
        if (length(which(ff.sp %in% canon)) == 0) {
          ff.sp[hh] <- canon
          if (!is.null(ll)) {
            names(ll)[sp][hh] <- ff.sp[hh]
            names(lengths)[sp][hh] <- ff.sp[hh]
          }
        } else {
          ff.sp <- c(ff.sp[-hh], ff.sp[hh])
          if (!is.null(ll)) {
            if (!isTRUE(all.equal(ff.sp, names(ll)))) {
              ll <- ll[match(ff.sp, names(ll))]
              lengths <- lengths[match(ff.sp, names(lengths))]
            }
          }
        }
      }
    }
  }

  if (length(ff.sp) != length(ff)) ff.trembl <- ff[-sp]
  
  ll.sp=NULL; ll.trembl=NULL
  
  ## Among Swiss-Prot entries, longer proteins come first.
  if (!is.null(ff.sp) & !is.null(ll)) {
    hh <- grep("[-_][0-9]+$", ff.sp)
    if (length(hh) > 0) {
      canon <- ff.sp[-hh]
      iso <- ff.sp[hh]
      ll.canon <- as.numeric(ll[match(canon, names(ll))])
      canon <- canon[order(ll.canon, decreasing = TRUE, na.last = TRUE)]
      ll.iso <- as.numeric(ll[match(iso, names(ll))])
      iso <- iso[order(ll.iso, decreasing = TRUE, na.last = TRUE)]
      ff.sp <- c(canon, iso)
    } else {
      if (!isTRUE(all.equal(ff.sp, names(ll)))) {
        ll.sp <- as.numeric(ll[match(ff.sp, names(ll))])
      } else {
        ll.sp <- as.numeric(ll)
      }
      ff.sp <- ff.sp[order(ll.sp, decreasing = TRUE, na.last = TRUE)]
    }
    if (!isTRUE(all.equal(ff.sp, names(ll)))) {
      ll.sp <- as.numeric(ll[match(ff.sp, names(ll))])
      names(ll.sp) <- names(ll)[match(ff.sp, names(ll))]
    } else {
      ll.sp <- as.numeric(ll)
      names(ll.sp) <- names(ll)
    }
  }
  
  ## Among TrEMBL entries, longer proteins come first.
  if (!is.null(ff.trembl) & !is.null(ll)) {
    ll.trembl <- ll[match(ff.trembl, names(ll))]
    ll.trembl <- as.numeric(ll.trembl)
    oo <- order(ll.trembl, decreasing = TRUE)
    ff.trembl <- ff.trembl[oo]
    ll.trembl <- ll.trembl[oo] 
    if (!isTRUE(all.equal(ll.trembl, names(ll)))) {
      names(ll.trembl) <- names(ll)[match(ff.trembl, names(ll))]
    } else {
      names(ll.trembl) <- names(ll)
    }
  }
  
  if (!is.null(ff.sp) | !is.null(ff.trembl)) {
    ff <- paste0(c(ff.sp, ff.trembl), collapse = ";")
  }

  if (!is.null(ll.sp) | !is.null(ll.trembl)) {
    ll <- paste0(c(ll.sp, ll.trembl), collapse = ";")
    names(ll) <- paste0(c(names(ll.sp), names(ll.trembl)), collapse = ";")
  }
  
  ## Keep unique if possible. No deduplication if different lengths present.
  ff1 <- strsplit(ff, ";")[[1]]
  ff1.dup <- which(duplicated(ff1))
  if (any(ff1.dup)) {
    if (!is.null(ll)) {
      ll1 <- strsplit(ll, ";")[[1]]
      nn1 <- strsplit(as.character(names(ll)), ";")[[1]]
      ll1.dup <- which(duplicated(ll1))
      if (isTRUE(all.equal(ff1.dup, ll1.dup))) {
        ff1 <- ff1[-ff1.dup]
        ll1 <- ll1[-ll1.dup]
        nn1 <- nn1[-ll1.dup]
        ff <- paste0(ff1, collapse=";")
        ll <- paste0(ll1, collapse=";")
        names(ll) <- paste0(nn1, collapse=";")
      }
    } else {
      ff <- paste0(ff1[-ff1.dup], collapse=";")
    }
  }
  
  ## Reset lengths to original
  if (!is.null(ll)) {
    vv <- strsplit(as.character(ll), ";")[[1]]
    nn <- strsplit(as.character(names(ll)), ";")[[1]]
    if (!isTRUE(all.equal(names(lengths), nn))) {
      jj <- match(names(lengths), nn)
      vv[jj] <- unname(lengths)
    }
    ll <- paste0(vv, collapse=";")
  }
  
  if (verbose) message("playbase::rank_uniprots] Completed.\n")

  return(list(feature=ff, lengths=ll))
  
}


iconv2utf8 <- function(s) {
  iconv(s, to = "UTF-8//TRANSLIT", sub = "")
}

iconv2ascii <- function(s) {
  iconv(s, to = "ascii//TRANSLIT", sub = "")
}
