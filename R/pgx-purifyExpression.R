##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' @export
pgx.purifyExpression <- function(tumorX, normalX,
                                 method = c("nnlm", "nnmf", "isopurer", "demixt", "undo")) {
  alpha <- list()
  xhat <- list()

  method <- tolower(method)
  if ("nnlm" %in% method) {
    ## ----------------------------------------------------------------------
    ## NNLM (BigOmics's own method...)
    ## ----------------------------------------------------------------------

    pen <- rep(0, 3)
    res <- NNLM::nnlm(normalX, tumorX, alpha = pen)
    cf <- res$coefficients
    cf
    normal.frac <- (normalX %*% cf)
    alpha0 <- (1 - Matrix::colSums(normal.frac) / Matrix::colSums(tumorX))
    alpha0

    xhat[["nnlm"]] <- pmax(tumorX - normal.frac, 0)
    dim(x.hat)
    alpha[["nnlm"]] <- alpha0
  }

  if ("nnmf" %in% method) {
    ## ----------------------------------------------------------------------
    ## NNMF (as in vignette)
    ## ----------------------------------------------------------------------


    ## compute proportion of contaminant content using NNMF
    k <- 10
    res.nmf <- NNLM::nnmf(tumorX, k = k, init = list(W0 = normalX), check.k = FALSE)

    x.hat <- res.nmf$W[, 1:k, drop = FALSE] %*% res.nmf$H[1:k, , drop = FALSE]
    nnlm.alpha <- with(res.nmf, Matrix::colSums(x.hat) / Matrix::colSums(W %*% H))
    round(nnlm.alpha, 2)

    xhat[["nnmf"]] <- x.hat
    alpha[["nnmf"]] <- nnlm.alpha
  }

  if ("isopurer" %in% method) {
    ## ----------------------------------------------------------------------
    ## IsoPureR (quite slow...)
    ## https://cran.r-project.org/web/packages/ISOpureR/vignettes/ISOpureRGuide.pdf
    ## ----------------------------------------------------------------------

    ISOpureS1model <- ISOpureR::ISOpure.step1.CPE(tumorX, normalX)
    ISOpureS2model <- ISOpureR::ISOpure.step2.PPE(tumorX, normalX, ISOpureS1model)
    isopurer.alpha <- ISOpureS2model$alphapurities
    isopurer.alpha

    x.hat <- ISOpureS2model$cc_cancerprofiles
    dim(x.hat)
    alpha[["isopurer"]] <- isopurer.alpha
    xhat[["isopurer"]] <- x.hat
  }

  if ("demixt" %in% method) {
    ## ----------------------------------------------------------------------
    ## DeMixT (crashes often...)
    ## https://bioinformatics.mdanderson.org/main/DeMixT
    ## ----------------------------------------------------------------------
    ## ?DeMixT
    res <- DeMixT::DeMixT(data.Y = tumorX, data.comp1 = normalX, if.filter = FALSE)
    res$pi

    Matrix::head(res$decovExprT, 3) ## purified tumor data
    Matrix::head(res$decovExprN1, 3) ## normal contiminant profile
    Matrix::head(res$decovMu, 3)
    Matrix::head(res$decovSigma, 3)

    x.hat <- res$decovExprT

    demixt.alpha <- (1 - res$pi[1, ])
    alpha[["demixt"]] <- demixt.alpha
    xhat[["demixt"]] <- x.hat
  }


  if ("undo" %in% method) {
    ## ----------------------------------------------------------------------
    ## UNDO
    ## ----------------------------------------------------------------------

    ## load tumor stroma mixing tissue samples
    ## two_source_deconv(
    ##    X, lowper=0.4, highper=0.1, epsilon1=0.01,
    ##    epsilon2=0.01, A, S[,1], S[,2], return=0)

    res <- UNDO::two_source_deconv(
      tumorX,
      lowper = 0.8, highper = 0.1, epsilon1 = 0.4,
      epsilon2 = 0, A = NULL, S1 = normalX, S2 = NULL, return = 1
    )
    str(res)
    res
    res$Estimated_Mixing_Matrix
    undo.alpha <- res$Estimated_Mixing_Matrix[, 2]
    x.hat <- NULL

    alpha[["undo"]] <- undo.alpha
    xhat[["undo"]] <- NULL
  }

  return(list(alpha = alpha, xhat = xhat))
}
