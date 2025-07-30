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


## #' @title Impute missing values with NMF
## #'
## #' @param X Numeric matrix of proteomics data
## #' @param groups Factor indicating groups for samples
## #' @param k Number of NMF factors. Default 10.
## #' @param r Cutoff proportion of missing values to impute within group. Default 0.5.
## #'
## #' @return Matrix with imputed values
## #'
## #' @description Imputes missing values in a proteomics matrix using non-negative matrix factorization (NMF).
## #'
## #' @details This function takes a numeric matrix \code{X} containing proteomics data, with samples in columns.
## #' It uses \code{groups} to define sample groups.
## #'
## #' It first replaces missing values (0s) with NAs.
## #' For each sample group, any rows with more than \code{r} proportion of NAs are set to 0.
## #'
## #' NMF with \code{k} factors is then applied to impute the remaining NAs.
## #' The imputed matrix is returned.
## #'
## #' @export
## prot.nmfImpute <- function(X, groups, k = 10, r = 0.5) {
##   setZERO <- function(x, y) {
##     b <- tapply(x, y, function(a) {
##       if (mean(is.na(a)) >= r) a[is.na(a)] <- 0
##       a
##     })
##     names(b) <- NULL
##     unlist(b)[names(x)]
##   }
##   X[X == 0] <- NA
##   sum(is.na(X))
##   impX <- t(apply(X, 1, setZERO, y = groups))
##   sum(is.na(impX))
##   out <- NNLM::nnmf(impX[, ], k = k, check.k = 0)
##   hatX <- (out$W %*% out$H)
##   impX[is.na(impX)] <- hatX[is.na(impX)]
##   return(impX)
## }


## #' @title Median Imputation of Proteomics Data
## #'
## #' @description This function performs median imputation on proteomics data, replacing missing or zero values with the median of the non-missing values in the same group.
## #'
## #' @param X A numeric matrix of proteomics data, where rows represent proteins and columns represent samples.
## #' @param groups A factor or character vector indicating the group membership of each sample (column) in X.
## #'
## #' @details The function first replaces any missing values in X with zeros.
## #' Then, for each protein (row), the function calculates the median of the non-missing values in each group, using a log-transformation to handle skewed data.
## #' The function then replaces any zero values in X with the calculated median for their respective group.
## #'
## #' @return A numeric matrix of the same dimensions as X, with missing and zero values replaced by the median of their respective group.
## #'
## #' @export
## prot.medianImpute <- function(X, groups) {
##   logmedian <- function(x) {
##     exp(stats::median(log(x)))
##   }
##   medianImputeZERO <- function(x, y) {
##     b <- tapply(x, y, function(a) {
##       a[a == 0] <- logmedian(a)
##       a
##     })
##     names(b) <- NULL
##     unlist(b)[names(x)]
##   }
##   X[which(is.na(X))] <- 0
##   impX <- t(apply(X, 1, medianImputeZERO, y = groups))
##   return(impX)
## }


## #' @title Impute Missing Values by Downshifted Gaussian
## #'
## #' @description Imputes missing values in a numeric matrix by drawing random values from a downshifted Gaussian distribution.
## #'
## #' @param x Numeric matrix or vector containing missing values to impute.
## #' @param width Width of Gaussian distribution as fraction of standard deviation. Default is 0.3.
## #' @param downshift Number of standard deviations to downshift the mean. Default is 1.8.
## #' @param bycolumn Logical indicating whether to impute by column (default) or across entire matrix.
## #'
## #' @details This function imputes missing values by drawing from a Gaussian distribution with mean shifted down by \code{downshift} standard deviations, and width equal to \code{width} times the standard deviation.
## #'
## #' If \code{bycolumn} is TRUE, this is done separately for each column. The mean and standard deviation are calculated from the non-missing values in each column.
## #'
## #' If \code{bycolumn} is FALSE, the parameters are calculated across the entire matrix and missing values are imputed randomly.
## #'
## #' @return The numeric matrix \code{x} with missing values imputed.
## #'
## #' @examples
## #' \dontrun{
## #' mat <- matrix(rnorm(100), ncol = 10)
## #' mat[sample(100, 20)] <- NA
## #' imputed <- .RGimputation(mat)
## #' }
## #' @export
## .RGimputation <- function(x, width = 0.3, downshift = 1.8, bycolumn = TRUE) {
##   if (bycolumn == T) {
##     for (i in 1:ncol(x)) {
##       x[, i][is.na(x[, i])] <- stats::rnorm(sum(is.na(x[, i])), mean(as.numeric(x[, i]), na.rm = T) - downshift * stats::sd(as.numeric(x[, i]), na.rm = T), width * stats::sd(as.numeric(x[, i]), na.rm = T))
##     }
##   } else {
##     x[is.na(x)] <- stats::rnorm(sum(is.na(x)), mean(as.numeric(x), na.rm = T) - downshift * stats::sd(as.numeric(x), na.rm = T), width * stats::sd(as.numeric(x), na.rm = T))
##   }
##   return(x)
## }


## ======================================================================
## =================== SPECIAL SILAC FUNCTIONS ==========================
## ======================================================================


## #' @title Calculate Copy Number from SILAC Data
## #'
## #' @description This function calculates the copy number of proteins from SILAC data, using the molecular weight and mass of the cell.
## #'
## #' @param data A numeric matrix of SILAC data, where rows represent proteins and columns represent samples.
## #' @param mol.weight A numeric vector of molecular weights for each protein.
## #' @param y The mass of the cell in picograms (pg).
## #'
## #' @details The function first calculates the total intensity for each sample by summing the intensities across all proteins, ignoring missing values.
## #' Then, for each protein and sample, the function calculates the mass by multiplying the intensity by the mass of the cell and dividing by the total intensity.
## #' The function then converts the mass to moles by dividing by the molecular weight and Avogadro's number.
## #'
## #' @return A numeric matrix of the same dimensions as data, containing the calculated copy number for each protein and sample.
## #'
## #' @export
## silac.calcCopyNumber <- function(data, mol.weight, y) {
##   ## data should be a numeric matrix
##   ## mol is the molecular weight
##   ## y is the mass of the cell in PICOGRAM (pg)
##   total.intensity <- Matrix::colSums(data, na.rm = TRUE)
##   mass <- t(t(data) * y / total.intensity) * 1e-12

##   ## calculate moles
##   cn <- mass / (mol.weight * 1000) * 6.022140857e23
##   return(cn)
## }


## #' Fit a Weibull distribution to bivariate data
## #'
## #' @title Fit a bivariate Weibull distribution
## #'
## #' @param x A numeric vector for the x variable
## #' @param y A numeric vector for the y variable
## #'
## #' @return A list with the estimated Weibull parameters
## #'
## #' @description Fits a bivariate Weibull distribution to two input variables.
## #'
## #' @details This function takes two numeric vectors \code{x} and \code{y} and fits a bivariate Weibull distribution.
## #' The Weibull distribution has shape and scale parameters for each variable.
## #'
## #' Maximum likelihood estimation is used to estimate the 4 Weibull parameters from the input data.
## #' The estimated parameters are returned as a named list.
## #'
## #' @export
## fit.weibull2 <- function(x, y) {
##   y <- as.numeric(y)
##   x <- as.numeric(x)
##   if (all(y == 0)) {
##     xfit <- seq(0, 48, 1)
##     return(list(x = xfit, y = xfit * 0, t50 = Inf))
##   }
##   jj <- which(x == 0 | y > 0) ## skip zero measurements (outlier?) if not x==0
##   x <- x[jj]
##   y <- y[jj]
##   cdf.weibull <- function(x, lambda, k) (1 - exp(-(x / lambda)^k))
##   fit <- stats::nls(y ~ cdf.weibull(x, lambda, k),
##     start = list(lambda = 50, k = 1),
##     lower = list(lambda = 0.01, k = 0.001),
##     algorithm = "port"
##   )
##   xfit <- seq(0, 48, 1)
##   lambda0 <- stats::coef(fit)["lambda"]
##   k0 <- stats::coef(fit)["k"]
##   yfit <- cdf.weibull(xfit, lambda0, k0)
##   t50 <- lambda0 * log(2)**(1 / k0) ## time at 50%
##   list(x = xfit, y = yfit, t50 = t50)
## }

## ================================================================================
## =============================== END OF FILE ====================================
## ================================================================================
