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
                            "min", "zero", "nbavg")[1:3],
                          rf.ntree = 100,
                          plot = FALSE) {
  impX <- list()

  ## ------------ msImpute --------------
  if('msImpute' %in% method) {
    sel <- which( rowSums(!is.na(X)) >= 4 )
    msi <- try(msImpute::msImpute(X[sel,], method="v2", group=NULL))
    if(!"try-error" %in% class(msi)) {
      cx <- X
      cx[sel,] <- msi
      impX[['msImpute']] <- cx
    }
  }

  ## ------------ LLSimpute --------------
  if('LLS' %in% method) {  
    ii <- which( rowMeans(is.na(X)) < 1)
    jj <- which( colMeans(is.na(X)) < 1)    
    X1 <- X[ii,jj]
    X1 <- X1 + 1e-2*matrix(rnorm(length(X1)),nrow=nrow(X1),ncol=ncol(X1))    
    k <- min(10,dim(X1)-1)
    result <- try( pcaMethods::llsImpute(t(X1), k=k) )
    if(!"try-error" %in% class(result)) {
      resX <- t(pcaMethods::completeObs(result))
      resX[!is.na(X1)] <- X1[!is.na(X1)]
      cx <- X
      cx[ii,jj] <- resX
      impX[['LLS']] <- cx
    }
  }
  
  ## ------------ promor --------------
  ## methods: 'minProb', 'minDet', 'RF', and 'SVD'.
  if('SVD' %in% method) {
    ## pcaMethods::pca(..., method = 'svdImpute').
    ## SVD does not allow rows with all NA
    ii <- which( rowMeans(is.na(X)) < 1)
    jj <- which( colMeans(is.na(X)) < 1)    
    res <- try(pcaMethods::pca( t(X[ii,jj]), method = "svdImpute", nPcs=3))
    if(!"try-error" %in% class(res)) {    
      cx <- X * NA
      cx[ii,jj] <- t(pcaMethods::completeObs(res))
      impX[['SVD']] <- cx
    }
  }
  if('SVD2' %in% method) {
    impX[['SVD2']] <- svdImpute2(X, nv=3)
  }
  if('RF' %in% method) {
    ## missForest
    res <- missForest::missForest( as.data.frame(t(X)),
                                  maxiter = 10,
                                  ## parallelize = 'variables',
                                  ntree = rf.ntree )
    impX[['RF']] <- t(as.matrix(res$ximp))
  }
  
  ## ------------ MSnbase --------------

  msn.impute <- function(X, method) {
    cx <- X * NA
    sel <- 1:nrow(X)
    if(method[1] == 'bpca') sel <- which(rowMeans(is.na(X)) < 1)
    rd <- data.frame(name=rownames(X), ID=rownames(X))[sel,]
    cd <- data.frame(sample=colnames(X))
    rownames(cd) <- colnames(X)
    rownames(rd) <- rownames(X)[sel]  
    mset <- MSnbase::MSnSet(X[sel,], pData=cd, fData=rd)
    res <- try( MSnbase::impute(mset, method = method[1]) )    
    if(!"try-error" %in% class(res)) {
      cx[sel,] <- MSnbase::exprs(res)
    }
    if("try-error" %in% class(res)) {
      cx <- NULL
    }
    cx
  }

  ms.methods = c("bpca", "knn", "QRILC", "MLE", "MinDet", "MinProb",
    "min", "zero", "nbavg")  
  ms.methods2 <- intersect(method, ms.methods)
  m = "MLE"  ## good and fast
  m = "bpca" ## the best?
  for(m in ms.methods2) {
    impX[[m]] <- msn.impute(X, method=m)
  }  

  names(impX)
  if(length(impX)==0) return(NULL)
  
  ## ------------ meta --------------
  metaX <- lapply(impX, as.vector)
  metaX <- do.call(cbind, metaX)
  metaX <- rowMeans(metaX, na.rm=TRUE)    
  metaX <- matrix(metaX, nrow = nrow(X), ncol = ncol(X),
    dimnames = dimnames(X))

  ## any remaining NA we fill with col/row medians
  if(any(is.na(metaX))) {
    missing <- which(is.na(metaX),arr.ind=TRUE)
    row.mx <- apply(metaX, 1, median, na.rm=TRUE)
    col.mx <- apply(metaX, 2, median, na.rm=TRUE)  
    mx <- rowMeans(cbind(row.mx[missing[,1]], col.mx[missing[,2]]),na.rm=TRUE)
    metaX[missing] <- mx
  }
  
  if(plot) {
    impX[['metaX']] <- resX
    length(impX)
    nr <- ceiling(sqrt(length(impX)))
    nr
    par(mfrow=c(nr,nr), mar=c(5,5,4,1))
    xlim = c(0, 1.2*max(X,na.rm=TRUE))    
    jj <- which(is.na(X))
    for(n in names(impX)) {
      h <- hist(impX[[n]][,], breaks=100, xlim=xlim, xlab="log2(x)",
        main=n, col="grey80", border="grey80", cex.main=1.4)
      hist(impX[[n]][jj], add=TRUE, col="red", border="red", breaks=h$breaks)
    }
  }

  metaX
}

#' @export
svdImpute2 <- function (M, nv = 3, threshold = 0.001,
                        maxSteps = 100, verbose=FALSE) {
  missing <- which(is.na(M), arr.ind=TRUE)

  ## initialize missing values with col/row medians
  row.mx <- apply(M, 1, median, na.rm=TRUE)
  col.mx <- apply(M, 2, median, na.rm=TRUE)  
  mx <- rowMeans(cbind(row.mx[missing[,1]], col.mx[missing[,2]]),na.rm=TRUE)
  M[missing] <- mx

  ## do SVD iterations
  count <- 0
  error <- Inf
  MOld <- M
  nv <- min(nv,dim(M))
  while ((error > threshold) && (count < maxSteps)) {
    res <- irlba::irlba( M, nv = nv )
    imx <- res$u %*% (diag(res$d) %*% t(res$v))
    M[missing] <- imx[missing]
    count <- count + 1
    if (count > 0) {
      error <- sqrt(sum((MOld - M)^2)/sum(MOld^2))
      if (verbose) {
        cat("change in estimate: ", error, "\n")
      }
    }
    MOld <- M
  }
  return(M)
}















## ----------------------------------------------------------------------
## -------------------------- end of file -------------------------------
## ----------------------------------------------------------------------
