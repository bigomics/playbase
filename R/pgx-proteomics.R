##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


########################################################################
##
## PROTEOMICS AND SILAC FUNCTIONS
##
########################################################################








## ======================================================================
## =================== ROGER FUNCTIONS ==================================
## ======================================================================
##
##



k <- 10
#' @export
prot.nmfImpute <- function(X, groups, k = 10, r = 0.5) {
  setZERO <- function(x, y) {
    b <- tapply(x, y, function(a) {
      if (mean(is.na(a)) >= r) a[is.na(a)] <- 0
      a
    })
    names(b) <- NULL
    unlist(b)[names(x)]
  }
  X[X == 0] <- NA
  sum(is.na(X))
  impX <- t(apply(X, 1, setZERO, y = groups))
  sum(is.na(impX))
  out <- NNLM::nnmf(impX[, ], k = k, check.k = 0)
  hatX <- (out$W %*% out$H)
  impX[is.na(impX)] <- hatX[is.na(impX)]
  return(impX)
}

#' @export
prot.medianImpute <- function(X, groups) {
  logmedian <- function(x) {
    exp(median(log(x)))
  }
  medianImputeZERO <- function(x, y) {
    b <- tapply(x, y, function(a) {
      a[a == 0] <- logmedian(a)
      a
    })
    names(b) <- NULL
    unlist(b)[names(x)]
  }
  X[which(is.na(X))] <- 0
  impX <- t(apply(X, 1, medianImputeZERO, y = groups))
  return(impX)
}

### IMPUTATION by downshifted Gaussian
#' @export
.RGimputation <- function(x, width = 0.3, downshift = 1.8, bycolumn = T) {
  if (bycolumn == T) {
    for (i in 1:ncol(x)) {
      x[, i][is.na(x[, i])] <- stats::rnorm(sum(is.na(x[, i])), mean(as.numeric(x[, i]), na.rm = T) - downshift * stats::sd(as.numeric(x[, i]), na.rm = T), width * stats::sd(as.numeric(x[, i]), na.rm = T))
    }
  } else {
    x[is.na(x)] <- stats::rnorm(sum(is.na(x)), mean(as.numeric(x), na.rm = T) - downshift * stats::sd(as.numeric(x), na.rm = T), width * stats::sd(as.numeric(x), na.rm = T))
  }
  return(x)
}



## ======================================================================
## =================== SPECIAL SILAC FUNCTIONS ==========================
## ======================================================================
##
##







#' @export
silac.calcCopyNumber <- function(data, mol.weight, y) {
  ## data should be a numeric matrix
  ## mol is the molecular weight
  ## y is the mass of the cell in PICOGRAM (pg)
  total.intensity <- Matrix::colSums(data, na.rm = TRUE)
  mass <- t(t(data) * y / total.intensity) * 1e-12

  ## calculate moles
  cn <- mass / (mol.weight * 1000) * 6.022140857e23
  return(cn)
}


###
# FIT WIBULL
###

#' @export
fit.weibull2 <- function(x, y) {
  y <- as.numeric(y)
  x <- as.numeric(x)
  if (all(y == 0)) {
    xfit <- seq(0, 48, 1)
    return(list(x = xfit, y = xfit * 0, t50 = Inf))
  }
  jj <- which(x == 0 | y > 0) ## skip zero measurements (outlier?) if not x==0
  x <- x[jj]
  y <- y[jj]
  cdf.weibull <- function(x, lambda, k) (1 - exp(-(x / lambda)^k))
  fit <- stats::nls(y ~ cdf.weibull(x, lambda, k),
    start = list(lambda = 50, k = 1),
    lower = list(lambda = 0.01, k = 0.001),
    algorithm = "port"
  )
  xfit <- seq(0, 48, 1)
  lambda0 <- stats::coef(fit)["lambda"]
  k0 <- stats::coef(fit)["k"]
  yfit <- cdf.weibull(xfit, lambda0, k0)
  t50 <- lambda0 * log(2)**(1 / k0) ## time at 50%
  list(x = xfit, y = yfit, t50 = t50)
}

# ;samples=NULL;protein=p



# ;samples="NaiveRest_6h";main=NULL;minq=3
