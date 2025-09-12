##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##



NORMALIZATION.METHODS <- c("none", "mean.center", "median.center", "sum", "CPM", "TMM", "RLE", "RLE2", "quantile", "maxMedian", "maxSum")

#' @title Normalize count data
#'
#' @description
#' Normalizes a matrix of RNA-seq count data using different methods.
#'
#' @param x Matrix of count data, with genes in rows and samples in columns.
#' @param methods Normalization method(s) to use. Options are "none", "scale", "CPM", "TMM", "RLE".
#' @param keep.zero Logical indicating whether to retain zero counts. Default is TRUE.
#'
#' @details
#' This function normalizes a matrix of RNA-seq count data using different normalization methods:
#'
#' - "none": No normalization
#' - "scale": Scale normalization by column means
#' - "CPM": Counts per million
#' - "TMM": Trimmed mean of M-values (edgeR)
#' - "RLE": Relative log expression (DESeq2)
#'
#' Zero counts can be retained or set to NA after normalization based on the keep.zero parameter.
#'
#' @return
#' The normalized count matrix.
#'
#' @examples
#' \dontrun{
#' counts <- matrix(rnbinom(10000, mu = 10, size = 1), 100, 100)
#' normalized <- pgx.countNormalization(counts, c("TMM", "RLE"))
#' }
#' @export
pgx.countNormalization <- function(x, methods, ref = NULL, rm.zero = FALSE) {
  ## Column-wise normalization (along samples).
  ## x:        counts (linear)
  ## method:   single method

  ## NEED RETHINK: would be better to rewrite Normalization in log2-space (IK)

  methods <- methods[1]

  if (rm.zero) {
    minx <- min(x, na.rm = TRUE)
    which.zero <- which(x == minx, arr.ind = TRUE)
    x[which.zero] <- NA
  }

  for (m in methods) {
    if (m == "none") {
      x <- x
    } else if (m %in% c("scale", "mean.center")) {
      mx <- mean(colMeans(x, na.rm = TRUE))
      x <- t(t(x) / (1e-8 + colMeans(x, na.rm = TRUE))) * mx
    } else if (m == "median.center") {
      mx <- apply(x, 2, median, na.rm = TRUE)
      x <- t(t(x) / (1e-8 + mx)) * median(mx, na.rm = TRUE)
    } else if (m == "median.center.nz") {
      mx <- apply(x, 2, function(x) median(x[x > 0], na.rm = TRUE))
      x <- t(t(x) / (1e-8 + mx)) * median(mx, na.rm = TRUE)
    } else if (m == "sum") {
      sumx <- mean(colSums(x, na.rm = TRUE))
      x <- t(t(x) / (1e-8 + colSums(x, na.rm = TRUE))) * sumx
    } else if (m == "CPM") {
      x <- logCPM(x, log = FALSE)
    } else if (m == "CP10K") {
      x <- logCPM(x, total = 1e4, log = FALSE)
    } else if (m == "TMM") {
      ## normalization on total counts (linear scale)
      x <- normalizeTMM(x, log = FALSE) ## does TMM on counts (edgeR)
    } else if (m == "RLE") {
      ## normalization on total counts (linear scale)
      x <- normalizeRLE(x, log = FALSE, use = "deseq2") ## does RLE on counts (Deseq2)
    } else if (m == "RLE2") {
      ## normalization on total counts (linear scale)
      x <- normalizeRLE(x, log = FALSE, use = "edger") ## does RLE on counts
    } else if (m == "quantile") {
      new.x <- 0.001 * limma::normalizeQuantiles(as.matrix(1000 * x)) ## shift to avoid clipping
      rownames(new.x) <- rownames(x)
      colnames(new.x) <- colnames(x)
      x <- new.x
    } else if (m == "maxMedian") {
      x <- maxMedianNormalization(x, toLog = FALSE)
    } else if (m == "maxSum") {
      x <- maxSumNormalization(x, toLog = FALSE)
    } else if (m == "reference") {
      x <- referenceNormalization(x, ref = ref, toLog = FALSE)
    } else {
      stop("pgx.countNormalization: unknown method")
    }
  }

  if (rm.zero) {
    ## put back zeros as zeros
    x[which.zero] <- minx
  }
  return(x)
}

#' @title Normalize expression data
#'
#' @description
#' Normalizes a matrix of RNA-seq/proteomics data using different methods. Uses log2 data.
#'
#' @param x Matrix of log2-transformed data, with genes in rows and samples in columns.
#' @param method Normalization method(s) to use.
#' @param ref Reference for referenceNormalization
#' @param prior Offset added in X prior to log2
#'
#' @details
#' This function normalizes a log2-transformed data matrix using different normalization methods:
#'
#' @return
#' Normalized, log2-transformed matrix.
#'
#' @export
normalizeExpression <- function(X, method = "CPM", ref = NULL, prior = 1) {
  ## Column-wise normalization (along samples).
  ## x:      log2-transformed counts + prior
  ## method: single method
  ## prior:  prior used for log2-transformation.
  m <- method
  methods <- c(
    "CPM", "quantile", "CPM+quantile", "TMM",
    "maxMedian", "maxSum", "reference"
  )
  if (!m %in% methods) {
    stop("[normalizeExpression]: unknown mormalization method")
  }
  counts <- 2**X - prior
  if (m == "CPM") {
    X <- logCPM(counts = counts, total = 1e+06, prior = prior, log = TRUE)
  } else if (m == "quantile") {
    X <- limma::normalizeQuantiles(X)
  } else if (m == "CPM+quantile") {
    X <- logCPM(counts = counts, total = 1e+06, prior = prior, log = TRUE)
    X <- limma::normalizeQuantiles(X)
  } else if (m == "TMM") {
    tmm <- normalizeTMM(counts, log = FALSE, method = "TMM")
    X <- log2(tmm + prior) ## we do log separate, not by TMM because prior
  } else if (m == "maxMedian") {
    # X <- maxMedianNormalization(counts = counts, prior = prior, toLog = TRUE)
    X <- maxMedianNormalization.logX(X)
  } else if (m == "maxSum") {
    # X <- maxSumNormalization(counts = counts, prior = prior, toLog = TRUE)
    X <- maxSumNormalization.logX(X)
  } else if (m == "reference") {
    # X <- referenceNormalization(counts = counts, ref = ref, prior = prior, toLog = TRUE)
    X <- referenceNormalization.logX(X, ref = ref)
  }
  return(X)
}

#'
#' @export
normalizeMultiOmics <- function(X, method = "median") {
  if (!(method %in% c("median", "combat"))) {
    message("[normalizeMultiOmics] WARNING. skipping normalization. unknown method: ", method)
    return(X)
  }

  ntype <- NA
  if (!all(grepl(":", rownames(X)))) {
    message("[normalizeMultiOmics] WARNING. not multi-omics data.")
    dtype <- rep("gx", nrow(X))
    ntype <- 1
  } else {
    dtype <- mofa.get_prefix(rownames(X))
    ntype <- length(unique(dtype))
  }

  if (ntype == 1 || method == "median") {
    global.median <- median(X[X > min(X, na.rm = TRUE)], na.rm = TRUE)
    for (dt in unique(dtype)) {
      ii <- which(dtype == dt)
      xx <- X[ii, ]
      minx <- min(xx, na.rm = TRUE) ## baseline. usually zero.
      xx[which(xx == minx)] <- NA
      mx <- matrixStats::colMedians(xx, na.rm = TRUE)
      X[ii, ] <- t(t(X[ii, ]) - mx) + global.median
    }
  }

  if (ntype > 1 && method == "combat") {
    X <- t(sva::ComBat(t(X), batch = dtype))
  }

  return(X)
}



## ====================================================================
## Normalization methods
## ====================================================================


#' Log-counts-per-million transformation
#'
#' @param counts Numeric matrix of read counts, with genes in rows and samples in columns.
#' @param total Total count to scale to. Default is 1e6.
#' @param prior Pseudocount to add prior to log transform. Default is 1.
#'
#' @return Matrix of log-transformed values.
#'
#' @details Transforms a matrix of read counts to log-counts-per-million (logCPM).
#' Adds a pseudocount \code{prior} (default 1) before taking the log transform.
#' Values are scaled to \code{total} counts (default 1e6).
#'
#' This stabilizes variance and normalizes for sequencing depth.
#'
#' @examples
#' \dontrun{
#' counts <- matrix(rnbinom(100 * 10, mu = 100, size = 1), 100, 10)
#' logcpm <- logCPM(counts)
#' }
#' @export
logCPM <- function(counts, total = 1e6, prior = 1, log = TRUE) {
  ## Transform to logCPM (log count-per-million) if total counts is
  ## larger than 1e6, otherwise scale to previous avarage total count.
  ##
  ##
  if (is.null(total)) {
    total0 <- mean(Matrix::colSums(counts, na.rm = TRUE)) ## previous sum
    total <- ifelse(total0 < 1e6, total0, 1e6)
    message("[logCPM] setting column sums to = ", round(total, 2))
  }
  if (any(class(counts) == "dgCMatrix")) {
    ## fast/sparse calculate CPM
    cpm <- counts
    cpm[is.na(cpm)] <- 0 ## OK??
    cpm@x <- total * cpm@x / rep.int(Matrix::colSums(cpm), diff(cpm@p)) ## fast divide by columns sum
    if (log) cpm@x <- log2(prior + cpm@x)
    return(cpm)
  } else {
    totcounts <- Matrix::colSums(counts, na.rm = TRUE)
    ## cpm <- t(t(counts) / totcounts * total)
    cpm <- sweep(counts, 2, totcounts, FUN = "/") * total
    if (log) cpm <- log2(prior + cpm)
    return(cpm)
  }
}

## #' @export
## edgeR.normalizeCounts.DEPRECATED <- function(M, method = c("TMM", "TMMwsp", "RLE", "upperquartile", "none")) {
##   method <- method[1]
##   dge <- edgeR::DGEList(M)
##   dge <- edgeR::calcNormFactors(dge, method = method)
##   edgeR::cpm(dge, log = TRUE)
## }


#' @title Normalize counts with TMM method
#'
#' @param counts Count matrix with genes in rows and samples in columns
#' @param log Logical, whether to log transform normalized counts
#' @param method Normalization method, default "TMM"
#'
#' @return Matrix of normalized counts
#'
#' @description
#' Normalizes a count matrix using the TMM (Trimmed Mean of M-values) method.
#'
#' @details
#' This function takes a count matrix and normalizes the counts using the TMM method
#' implemented in the edgeR package.
#'
#' It converts the count matrix to a DGEList object. Then it calculates normalization
#' factors using the calcNormFactors function with default TMM method.
#'
#' The normalized counts are returned as CPM values, optionally log-transformed.
#'
#' @examples
#' \dontrun{
#' counts <- matrix(rnbinom(100, mu = 10, size = 1), ncol = 10)
#' norm_counts <- normalizeTMM(counts)
#' }
#' @export
normalizeTMM <- function(counts, log = FALSE, prior = 1,
                         method = c("TMM", "TMMwsp", "RLE", "upperquartile", "none")[1]) {
  method <- method[1]
  dge <- edgeR::DGEList(as.matrix(counts), group = NULL)
  dge <- edgeR::calcNormFactors(dge, method = method)
  edgeR::cpm(dge, log = log, prior.count = prior)
}


#' Normalize count data using RLE method
#'
#' @description This function normalizes count data using the Relative Log Expression (RLE) method.
#' It supports normalization using either the edgeR or DESeq2 package.
#'
#' @param counts A numeric matrix of count data, where rows represent features and columns represent samples.
#' @param log Logical, whether to return log-transformed normalized counts.
#' @param use Character, the package to use for normalization, either "edger" or "deseq2".
#'
#' @return A numeric matrix of normalized count data.
#'
#' @examples
#' \dontrun{
#' counts <- matrix(rnbinom(100, mu = 10, size = 1), ncol = 10)
#' norm_counts <- normalizeRLE(counts)
#' }
#' @export
normalizeRLE <- function(counts, log = FALSE, use = "deseq2") {
  outx <- NULL
  if (use == "edger") {
    dge <- edgeR::DGEList(as.matrix(counts), group = NULL)
    dge <- edgeR::calcNormFactors(dge, method = "RLE")
    outx <- edgeR::cpm(dge, log = log)
  } else if (use == "deseq2") {
    cts <- round(counts)
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = cts,
      colData = cbind(colnames(cts), 1),
      design = ~1
    )
    disp <- DESeq2::estimateSizeFactors(dds)
    outx <- DESeq2::counts(disp, normalized = TRUE)
  } else {
    stop("unknown method")
  }
  return(outx)
}

#' Scaled log2(x/s+1) function. A generalization of log2(x+1) where we
#' allow a scaling factor s. This is also related to log2(x+p) with
#' custom prior p. However log1s maintains a minimum value of log1s=0
#' for all values of s.
#'
#' @description This function normalizes count data using the scaled
#'   log1p method.
#'
#' @param counts A numeric matrix of count data, where rows represent
#'   features and columns represent samples.
#' @param q numeric value corresponding to quantile to be used for
#'   scaling. Default q=0.20.
#'
#' @return A numeric matrix of log-normalized count data.
#'
#' @examples
#' \dontrun{
#' counts <- matrix(rnbinom(100, mu = 10, size = 1), ncol = 10)
#' norm_counts <- log1s(counts, q = 0.01)
#' }
#' @export
log1s <- function(counts, q = 0.20) {
  min.counts <- min(counts, na.rm = TRUE)
  jj <- which(counts > min.counts) ## only non-zero
  scale <- quantile(counts[jj], probs = q, na.rm = TRUE)
  log2(counts / scale + 1)
}

#' Normalize a multi-omics count matrix with log1s per datatype.
#'
#' @export
mofa.log1s <- function(counts, q = 0.20) {
  if (!is.multiomics(rownames(counts))) {
    X <- log1s(counts = counts, q = q)
    return(X)
  }
  dtype <- mofa.get_prefix(rownames(counts))
  X <- counts * 0
  for (d in unique(dtype)) {
    ii <- which(dtype == d)
    X1 <- counts[ii, ]
    x1.medians <- apply(X1, 2, function(x) median(x[x > 0], na.rm = TRUE))
    X1 <- t(t(X1) / x1.medians)
    X[ii, ] <- log1s(X1, q = q)
  }
  X
}

slog <- function(x, s = 1, q = NULL) {
  if (!is.null(q)) {
    s <- quantile(x[x > 0], probs = q, na.rm = TRUE)[1]
  }
  log2(s + x) - log2(s)
}

## #' @export
## safe.logCPM <- function(x, t = 0.05, prior = 1, q = NULL) {
##   qq <- apply(x, 2, quantile, probs = c(t, 1 - t), na.rm = TRUE)
##   jj <- which(t(t(x) < qq[1, ] | t(x) > qq[2, ]))
##   ax <- x
##   ax[jj] <- NA
##   colSums(x, na.rm = TRUE)
##   totx <- colSums(ax, na.rm = TRUE)
##   meanx <- colMeans(ax, na.rm = TRUE)
##   nnax <- colSums(!is.na(ax))
##   cpm <- sweep(x, 2, totx, FUN = "/") * 1e6
##   colSums(cpm, na.rm = TRUE)
##   ## slog(cx, s=prior, q=q)
##   log2(prior + cpm)
## }

## #' @export
## global_scaling <- function(X, method, shift = "clip") {
##   X[is.infinite(X)] <- NA

##   ## ---logMM & logMS created for MPoC
##   ## ---logCPM conformed to existing deployed master
##   if (method == "maxMedian") {
##     X <- maxMedianNormalization(counts = 2**X - 1)
##   } else if (method == "maxSum") {
##     X <- maxSumNormalization(counts = 2**X - 1)
##   } else if (method == "cpm") {
##     X <- logCPM(counts = 2**X - 1, log = TRUE)
##     ## median.tc <- median(colSums(2**X, na.rm = TRUE), na.rm = TRUE)
##     ## a <- log2(median.tc) - log2(1e6)
##     ## zero.point <- a
##   } else {
##     zero.point <- 0
##     which.zero <- which(X == 0)
##     if (grepl("^m[0-9]", method)) {
##       ## median centering
##       mval <- as.numeric(substring(method, 2, 99))
##       a <- median(X, na.rm = TRUE) - mval
##       zero.point <- a
##     } else if (grepl("^z[0-9]+", method)) {
##       ## zero at z-distance from median
##       zdist <- as.numeric(substring(method, 2, 99))
##       m0 <- mean(apply(X, 2, median, na.rm = TRUE))
##       s0 <- mean(apply(X, 2, sd, na.rm = TRUE))
##       zero.point <- m0 - zdist * s0
##     } else if (grepl("^q[0.][.0-9]+", method)) {
##       ## direct quantile
##       probs <- as.numeric(substring(method, 2, 99))
##       zero.point <- quantile(X, probs = probs, na.rm = TRUE)
##     } else {
##       stop("unknown method = ", method)
##     }

##     message("[normalizeCounts] shifting values to zero: z = ", round(zero.point, 4))
##     if (shift == "slog") {
##       ## smooth log-transform. not real zeros
##       X <- slog(2**X, s = 2**zero.point)
##     } else if (shift == "clip") {
##       ## linear shift and clip. Induces real zeros
##       X <- pmax(X - zero.point, 0)
##     } else {
##       stop("unknown shift method")
##     }
##     ## put back zeros
##     X[which.zero] <- 0
##   }

##   return(X)
## }

#' @export
is.xxl <- function(X, z = 10) {
  ## sdx <- apply(X, 1, function(x) mad(x[x > 0], na.rm = TRUE))
  sdx <- matrixStats::rowSds(X, na.rm = TRUE)
  sdx[is.na(sdx)] <- 0
  sdx0 <- 0.8 * sdx + 0.2 * mean(sdx, na.rm = TRUE) ## moderated SD
  mx <- rowMeans(X, na.rm = TRUE)
  this.z <- (X - mx) / sdx0
  (abs(this.z) > z)
}

#' Normalization for TMT and Silac data: maxMedianNorm
#'
#' @description This function normalizes TMT and Silac data using max median
#'
#' @param counts Numeric matrix of raw intensities: genes in rows, samples in columns.
#' @param prior Pseudocount to add prior to log transform. Default is 1.
#'
#' @return Normalized matrix of log2-transformed values.
#'
#' @details Most used normalization methods for TMT/Silac proteomics data: make the median or sum of each sample column equal to the max median or sum.
#'
#' @examples
#' \dontrun{
#' counts <- matrix(rnbinom(100 * 10, mu = 100, size = 1), 100, 10)
#' norm_counts <- maxMedianNormalization(counts, toLog = TRUE)
#' }
#' @export
maxMedianNormalization <- function(counts, toLog = TRUE, prior = 0) {
  mx <- apply(counts, 2, median, na.rm = TRUE)
  counts <- t(t(counts) / mx) * max(mx, na.rm = TRUE)
  if (toLog) {
    X <- log2(prior + counts)
  } else {
    X <- counts
  }
  return(X)
}

maxMedianNormalization.logX <- function(X) {
  mx <- matrixStats::colMedians(X, na.rm = TRUE)
  X <- t(t(X) - mx) + max(mx, na.rm = TRUE)
  return(X)
}

#' Normalization for TMT and Silac data: maxSumNormalization
#'
#' @description This function normalizes TMT and Silac data using max sum of intensities
#'
#' @param counts Numeric matrix of raw intensities: genes in rows, samples in columns.
#' @param prior Pseudocount to add prior to log transform. Default is 1.
#'
#' @return Normalized matrix of log2-transformed values.
#'
#' @details Most used norm. for TMT/Silac data: make median or sum of each sample equal to max median or sum.
#'
#' @examples
#' \dontrun{
#' counts <- matrix(rnbinom(100 * 10, mu = 100, size = 1), 100, 10)
#' norm_counts <- maxSumNormalization(counts, toLog = TRUE)
#' }
#' @export
maxSumNormalization <- function(counts, toLog = TRUE, prior = 0) {
  mx <- colSums(counts, na.rm = TRUE)
  counts <- t(t(counts) / mx) * max(mx, na.rm = TRUE)
  if (toLog) {
    X <- log2(prior + counts)
  } else {
    X <- counts
  }
  return(X)
}

maxSumNormalization.logX <- function(X) {
  X1 <- X
  X1[is.infinite(X1)] <- NA
  mx <- log2(matrixStats::colSums2(2**X1, na.rm = TRUE))
  X <- t(t(X) - mx) + max(mx, na.rm = TRUE)
  return(X)
}

normalizeMeans.logX <- function(X) {
  X1 <- X
  X1[is.infinite(X1)] <- NA
  mx <- colMeans(X1, na.rm = TRUE)
  X <- t(t(X) - mx) + max(mx, na.rm = TRUE)
  return(X)
}

#' @export
referenceNormalization <- function(counts, ref = 0.01, toLog = TRUE, prior = 0,
                                   bring.back = TRUE) {
  if (is.null(ref) || (length(ref) == 1 && is.numeric(ref[1]))) {
    ## If no reference gene/proteins are given, we take the most
    ## 'stable' genes, i.e. lowest SD in rank.
    rnk.sd <- apply(apply(counts, 2, rank, na.last = "keep"), 1, sd, na.rm = TRUE)
    nref <- ceiling(ref * nrow(counts))
    ref <- head(names(sort(rnk.sd)), nref)
  }
  ref <- intersect(ref, rownames(counts))
  if (length(ref) == 0) {
    return(counts)
  }
  mx <- colMeans(counts[ref, , drop = FALSE], na.rm = TRUE)
  counts <- t(t(counts) / mx)
  if (bring.back) counts <- counts * mean(mx, na.rm = TRUE)
  if (toLog) {
    X <- log2(prior + counts)
  } else {
    X <- counts
  }
  return(X)
}

referenceNormalization.logX <- function(X, ref = 0.01, bring.back = TRUE) {
  if (is.null(ref) || (length(ref) == 1 && is.numeric(ref[1]))) {
    ## If no reference gene/proteins are given, we take the most
    ## 'stable' genes, i.e. lowest SD in rank.
    rnk.sd <- apply(apply(X, 2, rank, na.last = "keep"), 1, sd, na.rm = TRUE)
    nref <- ceiling(ref * nrow(X))
    ref <- head(names(sort(rnk.sd)), nref)
  }
  ref <- intersect(ref, rownames(X))
  if (length(ref) == 0) {
    return(X)
  }
  mx <- colMeans(X[ref, , drop = FALSE], na.rm = TRUE)
  X <- t(t(X) - mx)
  if (bring.back) X <- X + mean(mx, na.rm = TRUE)
  return(X)
}

# Faster quantile normalization. Adapted from
# https://www.spsanderson.com/steveondata/posts/2024-03-28/index.html
fast_qn <- function(.data) {
  .data <- as.matrix(.data) ## Rfast needs dense matrix
  data_sort <- Rfast::colSort(.data)
  row_means <- Rfast::rowmeans(data_sort)
  data_sort <- matrix(row_means,
    nrow = nrow(data_sort),
    ncol = ncol(data_sort),
    byrow = FALSE
  )
  index_rank <- t(matrixStats::colRanks(.data, ties.method = "average"))
  normalized_data <- matrix(nrow = nrow(.data), ncol = ncol(.data))
  for (i in 1:ncol(.data)) {
    normalized_data[, i] <- data_sort[index_rank[, i], i]
  }
  return(normalized_data)
}
