##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## Normalization methods
##
##

NORMALIZATION.METHODS <- c("none", "mean", "scale", "NC", "CPM", "TMM", "RLE", "RLE2", "quantile")

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
normalizeTMM <- function(counts, log = FALSE, method = "TMM") {
  dge <- edgeR::DGEList(as.matrix(counts), group = NULL)
  dge <- edgeR::calcNormFactors(dge, method = method)
  edgeR::cpm(dge, log = log)
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
pgx.countNormalization <- function(x, methods, keep.zero = TRUE) {
  ## Column-wise normalization (along samples).
  ##
  ## x:        counts (linear)
  ## method:   single method

  methods <- methods[1]
  which.zero <- which(x == 0, arr.ind = TRUE)

  for (m in methods) {
    if (m == "none") {
      ## normalization on individual mean
      x <- x
    } else if (m == "scale") {
      ## normalization on individual mean
      mx <- mean(x, na.rm = TRUE)
      x <- t(t(x) / colMeans(x, na.rm = TRUE)) * mx
    } else if (m == "CPM") {
      x <- t(t(x) / Matrix::colSums(x, na.rm = TRUE)) * 1e6
    } else if (m == "TMM") {
      ## normalization on total counts (linear scale)
      x <- normalizeTMM(x, log = FALSE) ## does TMM on counts (edgeR)
    } else if (m == "RLE") {
      ## normalization on total counts (linear scale)
      x <- normalizeRLE(x, log = FALSE, use = "deseq2") ## does RLE on counts (Deseq2)
    } else if (m == "RLE2") {
      ## normalization on total counts (linear scale)
      x <- normalizeRLE(x, log = FALSE, use = "edger") ## does RLE on counts (Deseq2)
      ##        } else if(m %in% c("upperquartile")) {
      ##            ## normalization on total counts (linear scale)
    } else if (m == "quantile") {
      new.x <- 0.01 * limma::normalizeQuantiles(as.matrix(100 * x)) ## shift to avoid clipping
      rownames(new.x) <- rownames(x)
      colnames(new.x) <- colnames(x)
      x <- new.x
    }
  } ## end of for method

  x <- pmax(x, 0) ## prevent negative values
  ## put back zeros as zeros
  if (keep.zero && nrow(which.zero) > 0) {
    x[which.zero] <- 0
  }

  return(x)
}


#' @export
edgeR.normalizeCounts <- function(M, method = c("TMM", "TMMwsp", "RLE", "upperquartile", "none")) {
  method <- method[1]
  dge <- edgeR::DGEList(M)
  dge <- edgeR::calcNormFactors(dge, method = method)
  edgeR::cpm(dge, log = TRUE)
}


## --------------------------- new stuff -----------------------------------------

#' @export
normalizeCounts <- function(counts,
                            q.zero = 0.01,
                            method = "sdc",
                            shift.method = "slog",
                            weighted.qn = TRUE,
                            svd.init = NULL,
                            plot = FALSE) {
  which.zero <- which(counts == 0)
  which.missing <- which(is.na(counts))
  nzero <- sum(length(which.zero))
  nmissing <- sum(length(which.missing))
  message("[normalizeCounts] ", nzero, " zero values")
  message("[normalizeCounts] ", nmissing, " missing values")

  hist2 <- function(x1, ...) {
    dx <- min(x1, na.rm = TRUE) + 0.05 * diff(range(x1, na.rm = TRUE))
    hist(x1[x1 > dx], ...)
  }
  main <- ""
  boxplothist <- function(x, main) {
    boxplot(x, main = main)
    x[is.infinite(x)] <- NA
    if (length(which.missing)) {
      hist2(x[-which.missing], breaks = 200)
      if (any(!is.na(x[which.missing]))) {
        hist(x[which.missing], breaks = 200, col = "red", add = TRUE, border = NA)
      }
    } else {
      hist2(x, breaks = 200)
    }
  }

  if (plot) {
    x <- log2(1e-10 + counts)
    boxplothist(log2(1e-10 + counts), main = "log(counts)")
  }

  xmedians <- apply(counts, 2, function(x) median(x[x > 0], na.rm = TRUE))
  xcounts <- t(t(counts) / xmedians * mean(xmedians))

  if (plot) {
    x <- log2(xcounts)
    boxplothist(log2(xcounts), main = "median centered (MC)")
  }

  S <- scaled.log1p(xcounts, q = 0.0001)
  ##  S <- log2(xcounts)  ## zeros become -Inf
  if (plot) {
    boxplothist(S, main = "quantile scaled (MC+QS)")
  }

  ## detect outliers (values)
  mdx <- apply(S, 1, median, na.rm = TRUE)
  sdx <- apply(S, 1, mad, na.rm = TRUE)
  sdx0 <- mean(sdx, na.rm = TRUE)
  sdx[is.na(sdx)] <- sdx0
  sdx <- (0.8 * sdx + 0.2 * sdx0) ## bayesian
  Z <- (S - mdx) / sdx
  smax <- mdx + 10 * sdx
  smin <- mdx - 10 * sdx
  smax[is.na(smax)] <- mean(smax, na.rm = TRUE)
  smin[is.na(smin)] <- mean(smin, na.rm = TRUE)
  is.outlier <- which(S > smax)
  noutlier <- length(is.outlier)
  message("[normalizeCounts] detected ", noutlier, " outlier values")
  if (noutlier > 0) {
    S[is.outlier] <- NA
  }

  ## impute missing values (and/or previous outliers)
  any(is.na(S))
  if (any(is.na(S))) {
    ## nv=floor(max(ncol(counts)/5,3));
    ## svd.init=NULL
    nmissing <- sum(is.na(S))
    message("[normalizeCounts] imputing ", nmissing, " missing values")
    S1 <- svdImpute2(S, nv = 3, init = svd.init, fill.empty = "sample")
    if (plot) {
      boxplothist(S1, main = "missing imputed (MC+QS+MI)")
    }
    S <- S1
  }

  ## clip range of values
  message("[normalizeCounts] clipping values")
  S <- pmin(pmax(S, smin), smax)
  if (plot) {
    boxplothist(S, main = "clipped range (MC+QS+MI+CR)")
  }

  ## regress on total counts. here or later??
  regress.tc <- FALSE
  if (regress.tc) {
    message("[normalizeCounts] regressing on total counts")
    totcounts <- scale(log2(colSums(2**S) + 1))
    S <- limma::removeBatchEffect(S, batch = NULL, covariates = totcounts)
  }

  ## Do quantile normalization. weighted will retain a bit of
  ## variation in the values.
  message("[normalizeCounts] apply quantile normalization")
  if (weighted.qn) {
    S <- 0.1 * S + 0.9 * limma::normalizeQuantiles(S)
  } else {
    S <- limma::normalizeQuantiles(S)
  }

  if (plot) {
    boxplothist(S, main = "quantile normalized (MC+QS+MI+CR+QN)")
  }

  ## shift to 1%
  smin <- min(S, na.rm = TRUE) + 0.05 * diff(range(S))
  smin
  zero.point <- 0
  S1 <- S
  S1[S <= smin] <- NA
  S1[which.missing] <- NA
  if (method == "sd") {
    ## SD is calculated on global feature values
    m0 <- mean(apply(S1, 2, median, na.rm = TRUE))
    s0 <- mean(apply(S1, 2, sd, na.rm = TRUE))
    zdist <- qnorm(1 - q.zero) ## SD distance
    zero.point <- m0 - zdist * s0
  } else if (method == "sdc") {
    ## SD is calculated on centered feature values
    C1 <- S1 - rowMeans(S1, na.rm = TRUE)
    m0 <- median(S1, na.rm = TRUE)
    s0 <- mean(apply(C1, 2, sd, na.rm = TRUE))
    zdist <- qnorm(1 - q.zero) ## SD distance
    zero.point <- m0 - zdist * s0
  } else if (method == "cpm") {
    median.tc <- median(colSums(2**S, na.rm = TRUE), na.rm = TRUE)
    median.tc
    a <- log2(median.tc / 1e6)
    median(colSums(2**(S - a), na.rm = TRUE), na.rm = TRUE)
    zero.point <- a
  } else if (grepl("^m[0-9]", method)) {
    ## median centering
    mval <- as.numeric(substring(method, 2, 99))
    a <- median(S1, na.rm = TRUE) - mval
    median(S1 - a, na.rm = TRUE)
    median(S - a, na.rm = TRUE)
    zero.point <- a
  } else if (method == "quantile") {
    ## determine by direct quantile
    zero.point <- quantile(S1, probs = q.zero, na.rm = TRUE)
  } else {
    stop("unknown method = ", method)
  }

  if (0) {
    colSums(2**S1, na.rm = TRUE)
    colSums(2**(S1 - zero.point), na.rm = TRUE)
    colSums(2**S, na.rm = TRUE)
    colSums(2**(S - zero.point), na.rm = TRUE)
    median(S[S > 0.1], na.rm = TRUE)
    median(S[S > 0.1] - zero.point, na.rm = TRUE)
    hist(S, breaks = 200)
    hist((S - zero.point), breaks = 200)
  }

  message("[normalizeCounts] shifting values to zero: z = ", round(zero.point, 4))
  if (shift.method == "slog") {
    ## smooth log-transform. not real zeros
    S <- slog(2**S, s = 2**zero.point)
  } else {
    ## linear shift. Induces real zeros
    S <- pmax(S - zero.point, 0)
  }

  if (plot) {
    boxplothist(S, main = "quantile shifted (MC+QS+MI+CR+QN+SH)")
  }

  S
}

scaled.log1p <- function(counts, q = 0.01) {
  min.counts <- min(counts, na.rm = TRUE)
  ii <- which(counts == min.counts) ## left censored values
  jj <- which(counts > min.counts) ## only non-zero
  scale <- 1 / quantile(counts[jj], probs = q, na.rm = TRUE)
  log2(scale * counts + 1)
}

q <- 0.01
slog <- function(x, s = 1, q = NULL) {
  if (!is.null(q)) {
    s <- quantile(x[x > 0], probs = q, na.rm = TRUE)[1]
  }
  log2(s + x) - log2(s)
}

#' @export
safe.logCPM <- function(x, t = 0.05, prior = 1, q = NULL) {
  qq <- apply(x, 2, quantile, probs = c(t, 1 - t), na.rm = TRUE)
  jj <- which(t(t(x) < qq[1, ] | t(x) > qq[2, ]))
  ax <- x
  ax[jj] <- NA
  colSums(x, na.rm = TRUE)
  totx <- colSums(ax, na.rm = TRUE)
  meanx <- colMeans(ax, na.rm = TRUE)
  nnax <- colSums(!is.na(ax))
  cpm <- sweep(x, 2, totx, FUN = "/") * 1e6
  colSums(cpm, na.rm = TRUE)
  ## slog(cx, s=prior, q=q)
  log2(prior + cpm)
}

