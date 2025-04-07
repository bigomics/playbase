##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2025 BigOmics Analytics SA. All rights reserved.
##

#' Create feature sets (like genesets but collections of
#' features/rownames) based on gene families and custom class columns
#' in gene annation. This allows to create sets where probes have no
#' symbol (like in metabolomics).
#' 
#' @export
pgx.getFeatureSets <- function(pgx, min.size=1) {
  
  ## Get default classification from gene families
  ftclass <- lapply( playdata::FAMILIES, function(x) {
    map_probes( pgx$genes, x, column='human_ortholog', target='rownames')
  })

  ## Update 'all' group with all symbols
  ftclass[['<all>']] <- pgx$genes$feature

  ## Add any class columns in gene annotation data
  if(any(grepl("class", tolower(colnames(pgx$genes))))) {
    sel <- grep("class", tolower(colnames(pgx$genes)))
    i <- sel[1]
    for(i in sel) {
      aa <- pgx$genes[,i]
      aa[ aa=='' ] <- NA
      ii <- which(!is.na(aa))
      aa <- paste0(colnames(pgx$genes)[i],":",aa)
      aa.gmt <- tapply( pgx$genes$feature[ii], aa[ii], function(x) setdiff(x,c(NA,"","-")) )
      ftclass <- c(ftclass, aa.gmt)
    }
  }

  ## filter
  ftclass <- ftclass[ sapply(ftclass,length) >= min.size]
  return(ftclass)
}


#' Create feature sets (like genesets but collections of features)
#' based on gene families and custom class columns in gene annation.
#' 
#' @export
pgx.featureSetScore <- function(pgx) {
  
  Y <- pgx$samples
  X <- pgx$X

  if("families" %in% names(pgx)) {
    sets <- pgx$families
  } else {
    sets <- pgx.getFeatureSets(pgx, min.size=3) 
  }
  
  ## the code below overwrittes user input, and should be removed
  cvar <- playbase::pgx.getCategoricalPhenotypes(Y, max.ncat = 999)
  ss <- "sample|patient|years|days|months|gender|repl"
  cvar <- grep(ss, cvar, invert = TRUE, value = TRUE, ignore.case=TRUE)
  Y <- Y[colnames(X), cvar, drop = FALSE]
  kk <- which(apply(Y, 2, function(y) length(unique(y)) > 1))
  Y <- Y[, kk, drop = FALSE]
  dim(Y)
  
  method <- "meta"

  sdtop=1000
  sdx <- matrixStats::rowSds(X, na.rm = TRUE)
  names(sdx) <- rownames(X)

  S <- matrix(NA, nrow = length(sets), ncol = ncol(Y))
  rownames(S) <- names(sets)
  colnames(S) <- colnames(Y)
  i=1
  for (i in 1:ncol(Y)) {

    grp <- as.character(Y[, i])
    score <- rep(NA, length(sets))
    names(score) <- names(sets)
    j=1
    for (j in 1:length(sets)) {
      pp <- sets[[j]]
      pp <- head(pp[order(-sdx[pp])], sdtop)      
      pp <- intersect(pp, rownames(X))
      X1 <- X[pp, , drop = FALSE]
      if (nrow(X1) == 0) { score[j]=NA; next }
      
      s1 <- s2 <- 1
      if (method %in% c("correlation", "meta")) {
        mx <- t(apply(X1, 1, function(x) tapply(x, grp, mean, na.rm = TRUE)))
        if (nrow(mx) == 0 || ncol(mx) == 0) next
        D <- 1 - cor(mx, use = "pairwise")
        diag(D) <- NA
        s1 <- mean(D, na.rm = TRUE)
      }
      
      if (method %in% c("p-value", "meta")) {
        jj <- which(!is.na(grp))
        design <- model.matrix(~ grp[jj])
        fit <- tryCatch(
        {
          suppressWarnings(limma::eBayes(limma::lmFit(X1[, jj, drop = FALSE], design)))
        },
        error = function(w) { NA }
        )
        if (all(is.na(fit))) { score[j]=NA; next }        
        suppressWarnings(suppressMessages(top <- limma::topTable(fit)))
        s2 <- mean(-log10(1e-99 + top$adj.P.Val), na.rm = TRUE)
      }
      f <- 1
      f <- (1 - exp(-(length(pp) / 10)**2)) ## penalize smaller sets
      score[j] <- f * (s1 * s2) ** ifelse(method == "meta", 0.5, 1)
    }    
    S[, i] <- score    
  }
  
  S[is.na(S)] <- 0
  
  return(S)
}


#' Compute featureset score using F-test on phenotypes.
#' 
pgx.featureSetScore.Ftest <- function(pgx) {
  
  Y <- pgx$samples
  X <- pgx$X

  if("families" %in% names(pgx)) {
    sets <- pgx$families
  } else {
    sets <- pgx.getFeatureSets(pgx, min.size=3) 
  }
  
  ## Get categorical phenotypes
  cvar <- playbase::pgx.getCategoricalPhenotypes(Y, max.ncat = 999)
  ss <- "sample|patient|years|days|months|gender|repl"
  cvar <- grep(ss, cvar, invert = TRUE, value = TRUE, ignore.case=TRUE)
  Y <- Y[colnames(X), cvar, drop = FALSE]
  kk <- which(apply(Y, 2, function(y) length(unique(y)) > 1))
  Y <- Y[, kk, drop = FALSE]
  dim(Y)

  sdtop=1000
  sdx <- matrixStats::rowSds(X, na.rm = TRUE)
  names(sdx) <- rownames(X)

  ## Calculate score for each featureset, for each phenotype
  S <- matrix(NA, nrow = length(sets), ncol = ncol(Y))
  rownames(S) <- names(sets)
  colnames(S) <- colnames(Y)
  i=1
  for (i in 1:ncol(Y)) {

    y <- as.character(Y[, i])
    jj <- which(!is.na(y))

    score <- rep(NA, length(sets))
    names(score) <- names(sets)
    j=1
    for (j in 1:length(sets)) {
      pp <- sets[[j]]
      pp <- head(pp[order(-sdx[pp])], sdtop)      
      pp <- intersect(pp, rownames(X))
      X1 <- X[pp,jj, drop = FALSE]
      if (nrow(X1) == 0) { score[j]=NA; next }

      res <- gx.limmaF(X1, y[jj], fdr=1, lfc=0)
      score[j] <- mean(-log10(1e-99 + res$P.Value))
    }    
    S[, i] <- score    
  }
  
  S[is.na(S)] <- 0
  S <- S[order(-rowMeans(S)),,drop=FALSE]
  
  return(S)
}

#' Compute featureset score using limma on contrasts
#' 
pgx.featureSetScore.byContrast <- function(pgx) {
  

  if("families" %in% names(pgx)) {
    sets <- pgx$families
  } else {
    sets <- pgx.getFeatureSets(pgx, min.size=3) 
  }

  Y <- pgx$contrasts
  X <- pgx$X
  
  sdtop=1000
  sdx <- matrixStats::rowSds(X, na.rm = TRUE)
  names(sdx) <- rownames(X)

  ## Calculate score for each featureset, for each phenotype
  S <- matrix(NA, nrow = length(sets), ncol = ncol(Y))
  rownames(S) <- names(sets)
  colnames(S) <- colnames(Y)
  i=1
  for (i in 1:ncol(Y)) {
    y <- as.character(Y[, i])
    jj <- which(!is.na(y))
    score <- rep(NA, length(sets))
    names(score) <- names(sets)
    j=1
    for (j in 1:length(sets)) {
      pp <- sets[[j]]
      pp <- head(pp[order(-sdx[pp])], sdtop)      
      pp <- intersect(pp, rownames(X))
      X1 <- X[pp,jj, drop = FALSE]
      if (nrow(X1) == 0) { score[j]=NA; next }
      res <- gx.limma(X1, y[jj], fdr=1, lfc=0)
      score[j] <- mean(-log10(1e-99 + res$P.Value))
    }    
    S[, i] <- score    
  }
  
  S[is.na(S)] <- 0
  S <- S[order(-rowMeans(S)),,drop=FALSE]
  
  return(S)
}


