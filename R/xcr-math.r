##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##



#' Calculate signed ranks. Rank on absolute value and add sign.
#'
#' @export
signedRank <- function(x) {
  sign(x) * rank(abs(x), na.last = "keep", ties.method = "random") / length(x)
}

#' Calculate columnwise signed ranks for matrices. Rank on absolute
#' value and add sign. Values are rescaled to [-1,1].
#'
#' @export
colSignedRanks <- function(x) {
  sign(x) * t(matrixStats::colRanks(abs(x),
    na.last = "keep",
    ties.method = "random"
  )) / nrow(x)
}

#' @title Calculate cosine similarity
#'
#' @param X Numeric matrix
#' @param Y Optional second numeric matrix to compare against X
#' @param method Optional specification for method to calculate cosine similarity
#'
#' @return Matrix of cosine similarity values
#'
#' @description Calculates the cosine similarity between two matrices.
#'
#' @details This function takes two numeric matrices \code{X} and \code{Y} and calculates the cosine similarity between them.
#' It can handle sparse matrices and missing values.
#'
#' Cosine similarity is defined as the cosine of the angle between two vectors.
#' It is calculated as the dot product of the vectors divided by the L2 norms of the vectors.
#'
#' If only \code{X} is provided, the cosine similarity between columns of \code{X} is calculated.
#' If both \code{X} and \code{Y} are provided, the cosine similarity between columns of \code{X} and \code{Y} is calculated.
#'
#' The \code{method} parameter allows specifying alternate methods to calculate the cosine similarity.
#' By default the cosine similarity formula is used directly.
#'
#' @export
tcosine_similarity <- function(X, Y = NULL, method = NULL) {
  if (!is.null(Y)) {
    return(cosine_similarity(Matrix::t(X), Y = Matrix::t(Y), method = method))
  }
  return(cosine_similarity(Matrix::t(X), Y = NULL, method = method))
}


#' @describeIn tcosine_similarity calculates the cosine similarity between two matrices X and Y.
#' @export
cosine_similarity <- function(X, Y = NULL, method = NULL) {
  ## sparse cosine: A.B / (|A||B|)
  ## handles sparse matrix and missing values
  X <- methods::as(X, "dgCMatrix")
  if (is.null(Y)) {
    ii <- which(Matrix::colMeans(is.na(X)) != 1)
    zX <- X[, ii]
    zX[is.na(zX)] <- 0
    s1 <- crossprod(!is.na(X[, ii]), zX**2)
    ss <- (s1 * t(s1))**0.5
    M <- crossprod(zX) / (ss + 1e-20)
    S <- matrix(NA, nrow = ncol(X), ncol = ncol(X))
    S[ii, ii] <- as.matrix(M)
    return(S)
  } else {
    Y <- methods::as(Y, "dgCMatrix")
    ii <- which(Matrix::colMeans(is.na(X)) != 1)
    jj <- which(Matrix::colMeans(is.na(Y)) != 1)
    zX <- X[, ii]
    zX[is.na(zX)] <- 0
    zY <- Y[, jj]
    zY[is.na(zY)] <- 0
    s1 <- crossprod(!is.na(X[, ii]), zY**2)
    s2 <- crossprod(!is.na(Y[, jj]), zX**2)
    ss <- (s1 * t(s2))**0.5
    M <- crossprod(zX, zY) / (ss + 1e-20)
    S <- matrix(NA, nrow = ncol(X), ncol = ncol(Y))
    S[ii, jj] <- as.matrix(M)
    return(S)
  }
}


#' @describeIn tcosine_similarity calculates a sparse cosine similarity matrix for the data matrix X
#' @param k: number of nearest neighbors (default 100)
#' @param th: similarity threshold (default 0.01)
#' @param block: block size for processing (default 100)
#' @param ties.method: method for breaking ties (default "random")
#' @param gpu: use GPU acceleration (default FALSE)
#' @export
tcosine.sparse <- function(X, k = 100, th = 0.01, block = 100, ties.method = "random",
                           gpu = FALSE) {
  dim(X)
  X <- X / (1e-20 + Matrix::rowSums(X**2, na.rm = TRUE)**0.5)
  nblock <- ceiling(nrow(X) / block)
  nblock
  idx <- c()
  i <- 1
  for (i in 1:nblock) {
    j0 <- (i - 1) * block + 1
    j1 <- min((i - 1) * block + block, nrow(X))
    jj <- j0:j1

    if (gpu == TRUE) {
      rho <- gpuTcrossprod(X[jj, ], X[, ])
    } else {
      rho <- tcrossprod(X[jj, ], X[, ])
    }
    if (ties.method == "random") {
      idx.j <- apply(rho, 1, function(x) Matrix::head(order(x, stats::rnorm(length(x)), decreasing = TRUE), k))
    } else {
      idx.j <- apply(rho, 1, function(x) Matrix::head(order(x, decreasing = TRUE), k))
    }
    idx.j <- as.vector(idx.j)
    idx.i <- as.vector(sapply(1:nrow(rho), rep, k))
    x <- rho[cbind(idx.i, idx.j)]
    idx.i <- idx.i + j0 - 1
    x1 <- cbind(idx.i, idx.j, x = x)
    if (th > 0) x1 <- x1[which(x1[, 3] > th), , drop = FALSE]
    if (nrow(x1) > 0) {
      ## add symmetric part
      x2 <- cbind(x1[, 2], x1[, 1], x = x1[, 3])
      x1 <- rbind(x1, x2)
      idx <- rbind(idx, x1)
    }
  }
  S <- Matrix::sparseMatrix(
    i = idx[, 1], j = idx[, 2], x = idx[, 3],
    dims = c(nrow(X), nrow(X)), use.last.ij = TRUE
  )
  dim(S)
  rownames(S) <- colnames(S) <- rownames(X)

  sp <- round(100 * mean(S != 0, na.rm = TRUE), digits = 2)
  sp.na <- round(100 * mean(is.na(S), na.rm = TRUE), digits = 2)
  cat("tcosine.sparse: sparsity=", sp, "%\n")
  if (sp.na > 0) cat("tcosine.sparse: NA=", sp.na, "%\n")

  return(S)
}


#' @title Encode values to length scales
#'
#' @description
#' Encodes a numeric vector to a length scale based on the logarithm of the values.
#'
#' @param x A numeric vector to encode.
#' @param r The resolution of the length scale. Default is 0.1.
#' @param a The scale parameter for length scale. Default is 0.25.
#'
#' @details
#' This function takes a numeric vector \code{x} and encodes the values to a length scale.
#' It first takes the log of the values in \code{x}, excluding any NA values.
#' It then constructs a set of length bins based on the min and max log values, with bin width determined by the \code{r} parameter.
#' Values in \code{x} are assigned to bins based on their log values.
#' The \code{a} parameter controls the scale of the length encoding.
#'
#' The output is a sparse matrix with rows corresponding to the input \code{x} and columns for the length bins.
#'
#' @return
#' A sparse matrix encoding the input values to lengths.
#'
#' @examples
#' \dontrun{
#' x <- rpois(100)
#' len <- length_encode(x)
#' }
#' @export
length_encode <- function(x, r = 0.1, a = 0.25) {
  x0 <- x
  x <- x[which(!is.na(x))]
  logx <- log(x)
  maxlen <- max(logx, na.rm = TRUE)
  minlen <- min(logx, na.rm = TRUE)

  dx <- log(1 + a * r)
  brks <- seq(minlen - dx, maxlen + dx, dx)
  ix <- 2 + as.integer(cut(logx, breaks = brks))
  endx <- max(ix, na.rm = TRUE) + 2
  M0 <- stats::model.matrix(~ 0 + factor(ix, levels = 1:endx))
  M <- matrix(NA, nrow = length(x0), ncol = ncol(M0))
  M[which(!is.na(x0)), ] <- M0
  M <- Matrix::Matrix(M, sparse = TRUE)
  colnames(M) <- NULL
  rownames(M) <- NULL
  if (r == 0) {
    S <- cosine_similarity(t(M))
    return(S)
  }
  n <- ceiling(r / dx)
  n <- 2 ## alway 2 for now..
  X <- Matrix::Matrix(M)
  i <- 1
  for (i in 1:n) {
    dM <- cbind(M[, -i:-1], matrix(0, nrow(M), ncol = i))
    X <- X + dM
  }
  for (i in 1:n) {
    dM <- cbind(matrix(0, nrow(M), ncol = i), M[, 1:(ncol(M) - i)])
    X <- X + dM
  }
  return(X)
}

if (0) {
  x[1] <- NA
  length_encode(x, r = 0.1, a = 0.5)
  length_encode(x, r = 0.1, a = 0.33)
  length_encode(x, r = 0.1, a = 0.25)
}


## ===================================================================================
## ===================================================================================
## ===================================================================================


#' @title Calculate tagged Hamming distance
#'
#' @description Calculate the Hamming distance between two sequences, taking into account tags.
#'
#' @param aa Character vector. The first sequence.
#' @param bb Character vector. The second sequence.
#' @param align Logical. Whether to align the sequences before calculating Hamming distance. Default is TRUE.
#'
#' @details This function calculates the Hamming distance between two character vector
#' sequences \code{aa} and \code{bb}. It first parses any tags in the sequences using \code{parse.tags}.
#' It then calculates the Hamming distance between the sequences indicated by any matching tags (e.g.
#' \code{cdr1}, \code{cdr2}). For tagged sequences, it will align the sequences before calculating Hamming
#' distance if \code{align=TRUE}. For non-tagged sequences, it calculates simple Hamming distance.
#'
#' @return Named numeric vector of Hamming distances for each matched tag.
#'
#' @export
tagged.hamming <- function(aa, bb, align = TRUE) {
  aligned.dist <- function(seq1, seq2) {
    seq1 <- gsub("[ ]", "", seq1)
    seq2 <- gsub("[ ]", "", seq2)
    seq.distance(c(seq1, seq2), verbose = 0)$distance[1, 2]
  }
  hamming.dist <- function(s, t) {
    s1 <- sapply(s, strsplit, split = " ")[[1]]
    t1 <- sapply(t, strsplit, split = " ")[[1]]
    n <- max(length(s1), length(t1))
    s1 <- c(s1, rep(" ", n))[1:n]
    t1 <- c(t1, rep(" ", n))[1:n]

    sum(s1 != t1)
  }
  tag.hamming0 <- function(a, b) {
    a1 <- parse.tags(a)
    b1 <- parse.tags(b)
    tags <- intersect(names(a1), names(b1))
    if (length(tags) == 0) {
      return(NA)
    }
    tags
    dist <- rep(NA, length(tags))
    names(dist) <- tags
    for (i in 1:length(tags)) {
      k <- tags[i]
      is.seq <- grepl("cdr|seq", k)
      if (is.seq && align) {
        dist[i] <- aligned.dist(a1[k], b1[k])
      } else if (is.seq && !align) {
        dist[i] <- hamming.dist(a1[k], b1[k])
      } else {
        dist[i] <- 1 * (a1[k] != b1[k])
      }
    }
    dist
  }
  a <- aa[1]
  b <- bb[1]
  aa.tags <- unique(unlist(lapply(aa, function(a) names(parse.tags(a)))))
  bb.tags <- unique(unlist(lapply(bb, function(b) names(parse.tags(b)))))
  all.tags <- sort(intersect(aa.tags, bb.tags))

  if (length(all.tags) == 0 || is.null(all.tags)) {
    return(NA)
  }
  if (length(aa) == 1 && length(bb) > 1) aa <- rep(aa, length(bb))
  if (length(bb) == 1 && length(aa) > 1) bb <- rep(bb, length(aa))
  D <- lapply(1:length(aa), function(i) tag.hamming0(aa[i], bb[i]))
  D <- t(sapply(D, function(e) e[match(all.tags, names(e))]))

  if (length(all.tags) == 1) {
    D <- t(D)
    rownames(D) <- NULL
  }
  colnames(D) <- all.tags
  D
}

#' Rescale values to 0-1 range
#'
#' @title Rescale values to 0-1 range
#'
#' @param x A numeric vector to rescale.
#' @param symm Logical indicating if rescaled values should be symmetrized around 0. Default is FALSE.
#'
#' @return A numeric vector with rescaled values.
#'
#' @description Rescales a numeric vector to a 0-1 range.
#'
#' @details This function takes a numeric vector \code{x} and rescales the values to lie between 0 and 1.
#' It subtracts the minimum value and divides by the range.
#' This transforms the values to a 0-1 range while maintaining relative differences.
#'
#' If \code{symm=TRUE}, the rescaled values are further transformed to be symmetrized around 0
#' by subtracting 0.5 and multiplying by 2.
#'
#' The rescaled vector is returned.
#'
#' @examples
#' \dontrun{
#' x <- rnorm(100)
#' x_scaled <- uscale(x)
#' }
#' @export
uscale <- function(x, symm = FALSE) {
  uscale.func <- function(x) (x - min(x)) / (max(x) - min(x) + 1e-99)
  if (is.matrix(x)) {
    y <- apply(x, 2, uscale.func)
    if(nrow(x)==1) {
      y <- matrix(y,nrow=1)
      dimnames(y) <- dimnames(x)
    }
  } else {
    y <- uscale.func(x)
  }
  y[is.na(y)] <- NA
  if (symm) y <- (y - 0.5) * 2
  return(y)
}


#' Test
#'
#' @export
duplicated.dgCMatrix <- function (dgCMat, MARGIN) {
  MARGIN <- as.integer(MARGIN)
  n <- nrow(dgCMat)
  p <- ncol(dgCMat)
  J <- rep(1:p, diff(dgCMat@p))
  I <- dgCMat@i + 1
  x <- dgCMat@x
  if (MARGIN == 1L) {
    ## check duplicated rows
    names(x) <- J
    RowLst <- split(x, I)
    is_empty <- setdiff(1:n, I)
    result <- duplicated.default(RowLst)
  } else if (MARGIN == 2L) {
    ## check duplicated columns
    names(x) <- I
    ColLst <- split(x, J)
    is_empty <- setdiff(1:p, J)
    result <- duplicated.default(ColLst)
  } else {
    warning("invalid MARGIN; return NULL")
    result <- NULL
  }
  
  if(any(is_empty)){
    out <- logical(if(MARGIN == 1L) n else p)
    out[-is_empty] <- result
    if(length(is_empty) > 1)
      out[is_empty[-1]] <- TRUE
    result <- out
  }
  
  result
}

#' Fast test if rows are duplicated in sparse matrix M using random
#' vector and fast matrix multiplication.
#'
#' @export
rowsduplicated.dgCMatrix <- function(M) {
  rowMeans(apply(M %*% matrix(rnorm(2*ncol(M)),ncol=2),2,duplicated))==1
}


#' Calculate sparse correlation matrix handling missing values
#'
#' @param G Sparse matrix containing gene sets
#' @param rank_matrix Matrix of ranked values
#' @return Correlation matrix between G and rank_matrix
#' @details If rank_matrix has no missing values, calculates correlation directly using corSparse.
#' Otherwise computes column-wise correlations only using non-missing values.
#' @export
sparse_cor_matrix <- function(G, rank_matrix) {
  if (sum(is.na(rank_matrix)) == 0) {
    cor_matrix <- qlcMatrix::corSparse(G, rank_matrix)
  } else {
    message("rank matrix has missing values: computing column-wise reduced rankcor")
    rankcorSparse.vec <- function(X, y) {
      y <- y[!is.na(y)]
      gg <- intersect(rownames(X), names(y))
      y <- rank(y, na.last = "keep")
      qlcMatrix::corSparse(X[gg, , drop = FALSE], cbind(y[gg]))
    }
    cor_matrix <- lapply(1:ncol(rank_matrix), function(i) rankcorSparse.vec(G, rank_matrix[, i]))
    cor_matrix <- do.call(cbind, cor_matrix)
  }
  return(cor_matrix)
}

## ===================================================================================
## ============================== END OF FILE ========================================
## ===================================================================================
