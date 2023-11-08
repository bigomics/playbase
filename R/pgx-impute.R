##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## ----------------------------------------------------------------------
## -------------------------- subroutines -------------------------------
## ----------------------------------------------------------------------

#' @title Impute missing values
#' @description Generic function to impute missing (NA) value. Input is a matrix.
#'
#' @param X A matrix containing the input data
#' @param methods Methods for imputation. Default is c("LLS", "bpca", "msImpute")
#' @param plot Boolean for generating histogram plots
#'
#' @return An updated matrix containing imputed values
#'
#' @export
imputeMissing <- function(X,
                          method = c(
                            "LLS", "bpca", "msImpute", "SVD", "SVD2", "RF",
                            "knn", "QRILC", "MLE", "MinDet", "MinProb",
                            "min", "zero", "nbavg"
                          )[1:3],
                          rf.ntree = 100,
                          plot = FALSE) {
  impX <- list()

  ## ------------ msImpute --------------
  if ("msImpute" %in% method) {
    sel <- which(rowSums(!is.na(X)) >= 4)
    msi <- try(msImpute::msImpute(X[sel, ], method = "v2", group = NULL))
    if (!"try-error" %in% class(msi)) {
      cx <- X
      cx[sel, ] <- msi
      impX[["msImpute"]] <- cx
    }
  }

  ## ------------ LLSimpute --------------
  if ("LLS" %in% method) {
    ii <- which(rowMeans(is.na(X)) < 1)
    jj <- which(colMeans(is.na(X)) < 1)
    X1 <- X[ii, jj]
    X1 <- X1 + 1e-2 * matrix(rnorm(length(X1)), nrow = nrow(X1), ncol = ncol(X1))
    k <- min(10, dim(X1) - 1)
    result <- try(pcaMethods::llsImpute(t(X1), k = k))
    if (!"try-error" %in% class(result)) {
      resX <- t(pcaMethods::completeObs(result))
      resX[!is.na(X1)] <- X1[!is.na(X1)]
      cx <- X
      cx[ii, jj] <- resX
      impX[["LLS"]] <- cx
    }
  }

  ## ------------ promor --------------
  ## methods: 'minProb', 'minDet', 'RF', and 'SVD'.
  if ("SVD" %in% method) {
    ## pcaMethods::pca(..., method = 'svdImpute').
    ## SVD does not allow rows with all NA
    ii <- which(rowMeans(is.na(X)) < 1)
    jj <- which(colMeans(is.na(X)) < 1)
    res <- try(pcaMethods::pca(t(X[ii, jj]), method = "svdImpute", nPcs = 3))
    if (!"try-error" %in% class(res)) {
      cx <- X * NA
      cx[ii, jj] <- t(pcaMethods::completeObs(res))
      impX[["SVD"]] <- cx
    }
  }

  if ("SVD2" %in% method) {
    impX[["SVD2"]] <- svdImpute2(X, nv = 3)
  }

  if ("RF" %in% method) {
    ## missForest
    res <- missForest::missForest(as.data.frame(t(X)),
      maxiter = 10,
      ## parallelize = 'variables',
      ntree = rf.ntree
    )
    impX[["RF"]] <- t(as.matrix(res$ximp))
  }

  ## ------------ MSnbase --------------

  msn.impute <- function(X, method) {
    cx <- X * NA
    sel <- 1:nrow(X)
    if (method[1] == "bpca") sel <- which(rowMeans(is.na(X)) < 1)
    rd <- data.frame(name = rownames(X), ID = rownames(X))[sel, ]
    cd <- data.frame(sample = colnames(X))
    selX <- X[sel, ]
    rownames(selX) <- rownames(rd)
    colnames(selX) <- rownames(cd)
    mset <- MSnbase::MSnSet(selX, pData = cd, fData = rd)
    res <- try(MSnbase::impute(mset, method = method[1]))
    if (!"try-error" %in% class(res)) {
      cx[sel, ] <- MSnbase::exprs(res)
    }
    if ("try-error" %in% class(res)) {
      cx <- NULL
    }
    cx
  }

  ms.methods <- c(
    "bpca", "knn", "QRILC", "MLE", "MinDet", "MinProb",
    "min", "zero", "nbavg"
  )
  ms.methods2 <- intersect(method, ms.methods)
  m <- "MLE" ## good and fast
  m <- "bpca" ## the best?
  for (m in ms.methods2) {
    impX[[m]] <- msn.impute(X, method = m)
  }

  names(impX)
  if (length(impX) == 0) {
    return(NULL)
  }

  ## ------------ meta --------------
  metaX <- lapply(impX, as.vector)
  metaX <- do.call(cbind, metaX)
  metaX <- rowMeans(metaX, na.rm = TRUE)
  metaX <- matrix(metaX,
    nrow = nrow(X), ncol = ncol(X),
    dimnames = dimnames(X)
  )

  ## any remaining NA we fill with col/row medians
  if (any(is.na(metaX))) {
    missing <- which(is.na(metaX), arr.ind = TRUE)
    row.mx <- apply(metaX, 1, median, na.rm = TRUE)
    col.mx <- apply(metaX, 2, median, na.rm = TRUE)
    mx <- rowMeans(cbind(row.mx[missing[, 1]], col.mx[missing[, 2]]), na.rm = TRUE)
    metaX[missing] <- mx
  }

  if (plot) {
    impX[["metaX"]] <- resX
    length(impX)
    nr <- ceiling(sqrt(length(impX)))
    nr
    par(mfrow = c(nr, nr), mar = c(5, 5, 4, 1))
    xlim <- c(0, 1.2 * max(X, na.rm = TRUE))
    jj <- which(is.na(X))
    for (n in names(impX)) {
      h <- hist(impX[[n]][, ],
        breaks = 100, xlim = xlim, xlab = "log2(x)",
        main = n, col = "grey80", border = "grey80", cex.main = 1.4
      )
      hist(impX[[n]][jj], add = TRUE, col = "red", border = "red", breaks = h$breaks)
    }
  }

  metaX
}

#' @export
svdImpute2 <- function(M, nv = 3, threshold = 0.001, init = NULL,
                       maxSteps = 100, fill.empty = "median",
                       verbose = FALSE) {
  ## nv=3;threshold=0.001;init=NULL;maxSteps=100;fill.empty="median";verbose=FALSE

  ind.missing <- which(is.na(M), arr.ind = TRUE)
  empty.rows <- which(rowMeans(is.na(M)) == 1)
  empty.cols <- which(colMeans(is.na(M)) == 1)

  if (!is.null(init)) {
    ## initialize missing values with fixed value
    M[ind.missing] <- init
  } else {
    ## initialize missing values with col/row medians
    row.mx <- apply(M, 1, median, na.rm = TRUE)
    col.mx <- apply(M, 2, median, na.rm = TRUE)
    mx <- rowMeans(cbind(row.mx[ind.missing[, 1]], col.mx[ind.missing[, 2]]), na.rm = TRUE)
    M[ind.missing] <- mx
  }

  ## do SVD iterations
  count <- 0
  error <- Inf
  Mold <- M
  nv <- min(nv, dim(M))
  while ((error > threshold) && (count < maxSteps)) {
    res <- irlba::irlba(M, nv = nv)
    imx <- res$u %*% (diag(res$d) %*% t(res$v))
    M[ind.missing] <- imx[ind.missing]
    count <- count + 1
    if (count > 0) {
      error <- sqrt(sum((Mold - M)^2) / sum(Mold^2))
      if (verbose) {
        cat(count, ": change in estimate: ", error, "\n")
      }
    }
    Mold <- M
  }

  ## extra corrections (refill empty columns or rows)
  has.empty <- (length(empty.rows) > 0 || length(empty.cols) > 0)
  if (has.empty && fill.empty == "NA") {
    if (length(empty.rows)) M[empty.rows, ] <- NA
    if (length(empty.cols)) M[, empty.cols] <- NA
  }
  if (has.empty && fill.empty == "sample") {
    ii <- which(
      (!ind.missing[, 1] %in% empty.rows) &
        (!ind.missing[, 2] %in% empty.cols)
    )
    mm <- M[ind.missing[ii, ]]
    if (length(empty.rows) && length(mm)) {
      n1 <- length(M[empty.rows, ])
      message("[svdImpute2] warning: empty rows : n1 = ", n1)
      M[empty.rows, ] <- sample(mm, n1, replace = TRUE)
    }
    if (length(empty.cols) && length(mm)) {
      n2 <- length(M[, empty.cols])
      message("[svdImpute2] warning: empty cols : n2 = ", n2)
      M[, empty.cols] <- sample(mm, n2, replace = TRUE)
    }
  }

  return(M)
}




#' @describeIn knnImputeMissing Impute missing values with non-negative matrix factorization
#' @export
nmfImpute <- function(x, k = 5) {
  ## Impute missing values with NMF
  ##

  k <- min(k, dim(x))
  nmf <- NNLM::nnmf(x, k = k, check.k = FALSE, rel.tol = 1e-2, verbose = 0)
  xhat <- with(nmf, W %*% H)
  x[is.na(x)] <- xhat[is.na(x)]
  if (sum(is.na(x)) > 0) {
    nmf1 <- NNLM::nnmf(x, k = 1, check.k = FALSE, rel.tol = 1e-2, verbose = 0)
    xhat1 <- with(nmf1, W %*% H)
    x[is.na(x)] <- xhat1[is.na(x)]
  }
  x
}


#' @title Impute Missing Values with k-Nearest Neighbors
#'
#' @description This function imputes missing values in a vector using k-nearest neighbors.
#' @param x A numeric vector containing missing values to be imputed.
#' @param pos A matrix of positions for each element in `x`.
#' @param missing An optional value specifying the value used to represent missing values in `x`.
#' The default value is `NA`.
#' @param k An optional numeric value specifying the number of nearest neighbors to use for imputation.
#' The default value is 10.
#'
#' @details This function takes a numeric vector `x` containing missing values, a matrix of positions `pos`
#' for each element in `x`, and an optional value `missing` representing the missing values in `x` as input.
#' The function uses the k-nearest neighbors algorithm to impute the missing values in `x` based on their positions in `pos`.
#' The number of nearest neighbors used for imputation is specified by the `k` parameter.
#' The imputed values are returned as a numeric vector of the same length as `x`.
#'
#' @return A numeric vector of the same length as `x`, containing the imputed values.
#'
#' @export
knnImputeMissing <- function(x, pos, missing = NA, k = 10) {
  k0 <- which(x == missing)
  k1 <- which(x != missing)
  if (length(k0) == 0) {
    return(x)
  }
  pos0 <- pos[k0, ]
  pos1 <- pos[k1, ]
  nb <- FNN::get.knnx(pos1, pos0, k = k)$nn.index
  fx <- factor(x[k1])
  mx <- matrix(fx[as.vector(nb)], nrow = nrow(nb), ncol = ncol(nb))
  x.imp <- apply(mx, 1, function(x) names(which.max(table(x))))
  x[which(x == missing)] <- x.imp
  x
}



## ----------------------------------------------------------------------
## -------------------------- end of file -------------------------------
## ----------------------------------------------------------------------
