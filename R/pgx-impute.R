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
                            "LLS", "bpca", "msImpute", "SVD", "SVD2", "NMF", "NMF2",
                            "RF", "knn", "QRILC", "MLE", "MinDet", "MinProb",
                            "min", "zero", "nbavg", "rowmeans", "Perseus"
                          )[1:3],
                          rf.ntree = 100, nv = 5, keep.limits = FALSE,
                          infinite.na = TRUE, plot = FALSE) {
  ## set infinite as NA
  if (infinite.na) X[is.infinite(X)] <- NA

  impX <- list()

  ## ------------ simple rowmeans -----------
  if ("rowmeans" %in% method) {
    cx <- X
    ii <- which(is.na(cx), arr.ind = TRUE)
    cx[ii] <- rowMeans(cx, na.rm = TRUE)[ii[, 1]]
    ii <- which(is.na(cx), arr.ind = TRUE)
    cx[ii] <- colMeans(cx, na.rm = TRUE)[ii[, 2]]
    ii <- which(is.na(cx))
    cx[ii] <- mean(cx, na.rm = TRUE)
    impX[["rowmeans"]] <- cx
  }

  ## ------------ Perseus like --------------
  if ("Perseus" %in% method) {
    impX[["Perseus"]] <- perseusImpute(X, shift = 1.8, width = 0.3, method = "sample", seed = NULL)
  }

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
    result <- try(pcaMethods::llsImpute(t(X1), k = k, allVariables = TRUE))
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
    res <- try(pcaMethods::pca(t(X[ii, jj]), method = "svdImpute", nPcs = nv))
    if (!"try-error" %in% class(res)) {
      cx <- X * NA
      cx[ii, jj] <- t(pcaMethods::completeObs(res))
      impX[["SVD"]] <- cx
    }
  }

  if ("SVD2" %in% method) {
    impX[["SVD2"]] <- svdImpute2(X, nv = nv, init = "5%")
  }

  if ("NMF" %in% method) {
    minx <- min(X, na.rm = TRUE)
    X1 <- X - minx
    impX[["NMF"]] <- log(nmfImpute(exp(X1), k = nv)) + minx
  }

  if ("NMF2" %in% method) {
    minx <- min(X, na.rm = TRUE)
    X1 <- X - minx
    impX[["NMF2"]] <- log(nmfImpute2(exp(X1), k = nv)) + minx
  }

  if ("RF" %in% method) {
    ## missForest
    df <- as.data.frame(t(X), check.names = FALSE)
    colnames(df) <- paste0("rfvar.", 1:ncol(df))
    res <- missForest::missForest(df,
      maxiter = 10,
      ## parallelize = 'variables',
      ntree = rf.ntree
    )
    ximp <- as.matrix(t(res$ximp))
    cx <- X
    ii <- match(rownames(ximp), colnames(df))
    cx[ii, ] <- ximp
    impX[["RF"]] <- cx
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

  ## constrain limits??
  if (keep.limits) {
    minx <- min(X, na.rm = TRUE)
    maxx <- max(X, na.rm = TRUE)
    impX <- lapply(impX, function(x) pmin(pmax(x, minx), maxx))
  }

  metaX <- NULL
  if (length(impX) > 1) {
    ## ------------ meta --------------
    metaX <- lapply(impX, as.vector)
    metaX <- do.call(cbind, metaX)
    metaX <- rowMeans(metaX, na.rm = TRUE)
    metaX <- matrix(metaX,
      nrow = nrow(X), ncol = ncol(X),
      dimnames = dimnames(X)
    )
  } else {
    metaX <- impX[[1]]
  }

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

#' @title Impute missing values separately for each omics type in multi-omics data.
#' @description Generic function to impute missing (NA) value. Input is a matrix.
#' @description It uses the imputeMissing function.
#' @param X A multiomics matrix containing the input data.
#' @return Imputed multi-omics matrix
#' @export
imputeMissing.mox <- function(X,
                              method = c("LLS", "bpca", "msImpute", "SVD", "SVD2",
                                "NMF", "RF", "knn", "QRILC", "MLE", "MinDet", "MinProb",
                                "min", "zero", "nbavg", "rowmeans", "Perseus")[1:3],
                              rf.ntree = 100,
                              nv = 5,
                              keep.limits = FALSE,
                              infinite.na = TRUE,
                              plot = FALSE) {


  impX <- NULL
  is.mox <- playbase::is.multiomics(rownames(X))
  
  if (is.mox) {

    II <- list()
    dtypes <- unique(sub(":.*", "", rownames(X)))

    for (i in 1:length(dtypes)) {
      ii <- grep(paste0("^", dtypes[i], ":"), rownames(X))              
      message("[playbase::imputeMissing.mox]: Multi-omics data. Imputing ", dtypes[i])
      message("[playbase::imputeMissing.mox]: ", length(ii), " features..\n") 
      II[[dtypes[i]]] <- playbase::imputeMissing(X = X[ii, ], method = method,
        rf.ntree = rf.ntree, nv = nv, keep.limits = keep.limits,
        infinite.na = infinite.na, plot = plot)
    }

    impX <- do.call(rbind, II)
    rm(II); gc()
    
  } else {

    impX <- playbase::imputeMissing(X = X, method = method,
      rf.ntree = rf.ntree, nv = nv, keep.limits = keep.limits,
      infinite.na = infinite.na, plot = plot)

  }

  return(impX)

}


#' @export
svdImpute2 <- function(X, nv = 10, threshold = 0.001, init = NULL,
                       maxSteps = 100, fill.empty = "median",
                       infinite.na = TRUE,
                       ## randomize.init = FALSE,
                       verbose = FALSE) {
  ## nv=3;threshold=0.001;init=NULL;maxSteps=100;fill.empty="median";verbose=FALSE

  if (infinite.na) X[is.infinite(X)] <- NA
  ind.missing <- which(is.na(X), arr.ind = TRUE)
  empty.rows <- which(rowMeans(is.na(X)) == 1)
  empty.cols <- which(colMeans(is.na(X)) == 1)

  if (is.null(nv)) {
    nv <- max(1, round(mean(is.na(X)) * min(dim(X)))) ## heuristic..
    message(paste0("setting nv to ", nv))
  }
  nv <- min(nv, round(min(dim(X)) / 3))

  init.methods <- c("MinDet","MinProb","QRILC","min")
  
  if (!is.null(init) && is.character(init) && grepl("%", init)) {
    ## initialize missing values with quantile fixed value
    q <- as.numeric(sub("%", "", init))
    init <- quantile(X[!is.na(X)], probs = q * 0.01)[1]
    message(paste0("setting initial values to ", q, "%. init=", round(init, 4)))
    X[ind.missing] <- init
  } else if (!is.null(init) && is.character(init) &&
               init %in% init.methods ) {
    ## initialize missing values with other impute method
    message(paste("setting initial values using", init))
    initX <- imputeMissing(X, method=init)
    X[ind.missing] <- initX[ind.missing]
  } else {
    ## initialize missing values with col/row medians
    row.mx <- apply(X, 1, median, na.rm = TRUE)
    col.mx <- apply(X, 2, median, na.rm = TRUE)
    row.mx[is.nan(row.mx)] <- NA
    col.mx[is.nan(col.mx)] <- NA
    X[ind.missing] <- row.mx[ind.missing[, 1]]
    ind.missing2 <- which(is.na(X), arr.ind = TRUE)
    X[ind.missing2] <- col.mx[ind.missing2[, 2]]
  }

  ## randomize initial values?
  ## if (randomize.init) {
  ##   sdx <- mean(matrixStats::rowSds(X, na.rm = TRUE), na.rm = TRUE)
  ##   rr <- rnorm(nrow(ind.missing), mean = 0, sd = sdx)
  ##   X[ind.missing] <- X[ind.missing] + rr
  ## }

  ## do SVD iterations
  count <- 0
  error <- Inf
  Xold <- X
  nv <- min(nv, dim(X) - 1)
  while ((error > threshold) && (count < maxSteps)) {
    if (nv < min(dim(X)) / 5) {
      res <- irlba::irlba(X, nv = nv, nu = nv)
    } else {
      res <- svd(X, nv = nv, nu = nv)
      res$d <- res$d[1:nv]
    }
    if (nv == 1) {
      imx <- res$d * (res$u %*% t(res$v))
    } else {
      imx <- res$u %*% (diag(res$d) %*% t(res$v))
    }
    X[ind.missing] <- imx[ind.missing]
    count <- count + 1
    if (count > 0) {
      error <- sqrt(sum((Xold - X)^2, na.rm = TRUE) / sum(Xold^2, na.rm = TRUE))
      if (verbose) {
        cat(count, ": change in estimate: ", error, "\n")
      }
    }
    Xold <- X
  }

  ## extra corrections (refill empty columns or rows)
  has.empty <- (length(empty.rows) > 0 || length(empty.cols) > 0)
  if (has.empty && fill.empty == "NA") {
    if (length(empty.rows)) X[empty.rows, ] <- NA
    if (length(empty.cols)) X[, empty.cols] <- NA
  }
  if (has.empty && fill.empty == "sample") {
    ii <- which(
      (!ind.missing[, 1] %in% empty.rows) &
        (!ind.missing[, 2] %in% empty.cols)
    )
    mm <- X[ind.missing[ii, ]]
    if (length(empty.rows) && length(mm)) {
      n1 <- length(X[empty.rows, ])
      message("[svdImpute2] warning: empty rows : n1 = ", n1)
      X[empty.rows, ] <- sample(mm, n1, replace = TRUE)
    }
    if (length(empty.cols) && length(mm)) {
      n2 <- length(X[, empty.cols])
      message("[svdImpute2] warning: empty cols : n2 = ", n2)
      X[, empty.cols] <- sample(mm, n2, replace = TRUE)
    }
  }

  return(X)
}


#' @title Impute missing values with non-negative matrix factorization
#'
#' @description Imputes missing values in matrix non-negative matrix factorization
#'
#' @export
nmfImpute <- function(x, k = 5) {
  ## Impute missing values with NMF
  ##
  k <- min(k, dim(x) - 1)
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


#' @title Impute missing values with iterative NMF
#'
#' @description Imputes missing values in matrix with iterative
#'    non-negative matrix factorization
#' 
#' @export
nmfImpute2 <- function(x, k=5, niter=10, init=0.05) {
  k <- min(k, dim(x) - 1)
  ii <- which(is.na(x), arr.ind=TRUE)
  impx <- NULL
  if(is.numeric(init) && length(init)==1) {
    impx <- x
    minx <- min(x,na.rm=TRUE)  
    impx[ii] <- quantile(x[x > minx], probs=init, na.rm=TRUE)
  } else if(is.character(init) && length(init)==1) {
    impx <- imputeMissing(x, method=init, nv = k)
  } else if(is.matrix(init) && all(dim(init)==dim(x))) {
    impx <- init
  } 
  if(is.null(impx)) stop("[nmfImpute2] invalid init:", init)
  
  i=1
  for(i in 1:niter) {
    m1 <- RcppML::nmf(impx, k=k, verbose=FALSE)
    if(all(c("w","d","h") %in% names(m1))) {
      wh <- m1$w %*% diag(m1$d) %*% m1$h
    } else if(all(c("w","d","h") %in% slotNames(m1))) {
      wh <- m1@w %*% diag(m1@d) %*% m1@h
    } else {
      stop("[nmfImpute2] Fatal error in RcppML::nmf")
    }
    impx[ii] <- wh[ii]
  }
  impx
}


#' @title Impute missing values using Bayesian PCA
#'
#' @description Imputes missing values using BPCA
#'
#' @export
BPCAimpute <- function(X, k = 2) {
  k <- min(k, dim(X) - 1)
  ii <- which(rowMeans(is.na(X)) < 1)
  pc <- pcaMethods::pca(t(X[ii, ]), method = "bpca", nPcs = k)
  obsX <- t(pcaMethods::completeObs(pc))
  impX <- X
  impX[ii, ] <- obsX
  impX
}


## https://www.biorxiv.org/content/10.1101/2020.08.12.248963v1.full
## Perseus, by default, impute for each sample separately.



#' @title Perseus-style imputation
#' @description Imputes missing values using Perseus-style imputation.
#' Draw from a narrow Gaussian distribution, which is down-shifted toward
#' lower intensities to simulate values below the detection limit.
#' Default Perseus parameters are width=0.3 and down shift = 1.8 (log2 scale)
#' from the mean of the observed data distribution. This is commonly used in
#' proteomics where missing values often reflect low abundance. By imputing
#' from a left-shifted distribution (mu-shift*sigma), the function reflects
#' the assumption that missing values tend to be lower than observed ones.
#' @export
perseusImpute <- function(X,
                          shift = 1.8,
                          width = 0.3,
                          method = c("sample", "global"),
                          seed = NULL) {
  if (sum(is.na(X)) == 0) {
    return(X)
  }

  method <- method[1]
  if (!method %in% c("sample", "global")) {
    stop("[playbase::imputeMissing perseusImpute: method must be 'sample' or 'global']")
  }

  message("[perseusImpute: Performing Perseus-style imputation")

  impX <- X
  if (!is.null(seed)) set.seed(seed)

  min.pos <- min(impX, na.rm = TRUE)

  if (method == "sample") {
    mu0 <- colMeans(impX, na.rm = TRUE)
    sdx0 <- matrixStats::colSds(impX, na.rm = TRUE)
    mu <- mu0 - shift * sdx0
    sdx <- width * sdx0
    na.samples <- apply(impX, 2, function(x) sum(is.na(x)))
    na.samples <- names(na.samples[which(na.samples > 0)])
    mu <- unname(mu[match(na.samples, names(mu))])
    sdx <- unname(sdx[match(na.samples, names(sdx))])
    i <- 1
    for (i in 1:length(na.samples)) {
      nas <- which(is.na(impX[, na.samples[i]]))
      impX[nas, na.samples[i]] <- pmax(rnorm(length(nas), mean = mu[i], sd = sdx[i]), min.pos)
    }
  } else if (method == "global") {
    mu0 <- mean(as.numeric(impX), na.rm = TRUE)
    sdx0 <- sd(as.numeric(impX), na.rm = TRUE)
    mu <- mu0 - shift * sdx0
    sdx <- width * sdx0
    impX[which(is.na(impX))] <- pmax(rnorm(sum(is.na(impX)), mean = mu, sd = sdx), min.pos)
  }

  return(impX)
}


## ----------------------------------------------------------------------
## -------------------------- end of file -------------------------------
## ----------------------------------------------------------------------
