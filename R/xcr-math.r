##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##







sparse <- NULL



#' @export
tcosine_similarity <- function(X, Y = NULL, method = NULL) {
  if (!is.null(Y)) {
    return(cosine_similarity(t(X), Y = t(Y), method = method))
  }
  return(cosine_similarity(t(X), Y = NULL, method = method))
}

#' @export
cosine_similarity <- function(X, Y = NULL, method = NULL) {
  ## sparse cosine: A.B / (|A||B|)
  ## handles sparse matrix and missing values
  X <- as(X, "dgCMatrix")
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
    Y <- as(Y, "dgCMatrix")
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


#' @export
tcosine.sparse <- function(X, k = 100, th = 0.01, block = 100, ties.method = "random",
                           gpu = FALSE) {
  dim(X)
  X <- X / (1e-20 + Matrix::rowSums(X**2)**0.5)
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
      idx.j <- apply(rho, 1, function(x) Matrix::head(order(x, rnorm(length(x)), decreasing = TRUE), k))
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








#' @export
length_encode <- function(x, r = 0.1, a = 0.25) {
  x0 <- x
  x <- x[which(!is.na(x))]
  logx <- log(x)
  maxlen <- max(logx, na.rm = TRUE)
  minlen <- min(logx, na.rm = TRUE)
  c(minlen, maxlen)
  dx <- log(1 + a * r)
  brks <- seq(minlen - dx, maxlen + dx, dx)
  ix <- 2 + as.integer(cut(logx, breaks = brks))
  endx <- max(ix, na.rm = TRUE) + 2
  M0 <- model.matrix(~ 0 + factor(ix, levels = 1:endx))
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
  all.tags
  if (length(all.tags) == 0 || is.null(all.tags)) {
    return(NA)
  }
  if (length(aa) == 1 && length(bb) > 1) aa <- rep(aa, length(bb))
  if (length(bb) == 1 && length(aa) > 1) bb <- rep(bb, length(aa))
  D <- lapply(1:length(aa), function(i) tag.hamming0(aa[i], bb[i]))
  D <- t(sapply(D, function(e) e[match(all.tags, names(e))]))
  D
  if (length(all.tags) == 1) {
    D <- t(D)
    rownames(D) <- NULL
  }
  colnames(D) <- all.tags
  D
}


#' @export
uscale <- function(x, symm = FALSE) {
  uscale.func <- function(x) (x - min(x)) / (max(x) - min(x))
  if (NCOL(x) == 1) {
    y <- uscale.func(x)
  } else {
    y <- apply(x, 2, uscale.func)
  }
  y[is.na(y)] <- NA
  if (symm) y <- (y - 0.5) * 2
  return(y)
}








## NEED RETHINK!! for non centereds, non-scaled values!!!
## fun="mean";m=1000,n=1000

## ===================================================================================
## ===================================================================================
## ===================================================================================
