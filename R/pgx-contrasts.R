##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## ----------------------------------------------------------------------
## Auto-detection helper functions
## ----------------------------------------------------------------------


##' Get contrast groups from PGX object
#'
#' @title Get contrast groups
#'
#' @description Retrieves the group labels for a given contrast from a PGX object.
#'
#' @param pgx A PGX object.
#' @param contrast The index or name of the contrast.
#' @param as.factor Whether to return as a factor.
#'
#' @details This function extracts the experimental design matrix from the PGX object,
#' selects the specified contrast column, and returns the corresponding group labels.
#'
#' @return A vector or factor of group labels for the given contrast.
#'
#' @examples
#' \dontrun{
#' # TODO
#' }
#' @export
pgx.getContrastGroups <- function(pgx, contrast, as.factor = TRUE) {
  exp.matrix <- pgx$model.parameters$exp.matrix
  grp <- contrastAsLabels(exp.matrix[, contrast, drop = FALSE], as.factor = as.factor)
  if (NCOL(grp) == 1) {
    grp <- grp[, 1]
    names(grp) <- rownames(exp.matrix)
  }
  return(grp)
}






#' Get sample conditions from expression matrix
#'
#' @param exp.matrix Expression matrix with samples in columns
#' @param nmax Maximum number of groups to return
#'
#' @return Character vector of sample conditions
#'
#' @description
#' Extracts sample conditions from an expression matrix by pasting column names.
#'
#' @details
#' This function takes an expression matrix and generates a sample condition vector
#' by pasting the column names with an underscore.
#'
#' If there are more than nmax unique conditions, they are renamed to "group1", "group2", etc.
#'
#' @export
pgx.getConditions <- function(exp.matrix, nmax = 3) {
  group <- apply(exp.matrix, 1, paste, collapse = "_")
  group <- factor(group)
  if (ncol(exp.matrix) > nmax) {
    ngroup <- length(unique(group))
    levels(group) <- paste0("group", 1:ngroup)
  }
  return(as.character(group))
}

#' Generate design matrix for expression data
#'
#' @param pheno Data frame containing sample metadata
#' @param contr.matrix Contrast matrix defining comparisons of interest
#'
#' @return Design matrix for expression data analysis
#'
#' @description
#' Generates a design matrix for expression data based on sample metadata and contrasts.
#'
#' @details
#' This function takes a data frame of sample metadata and a contrast matrix as input.
#' It constructs a model matrix with columns for the metadata variables, followed by
#' columns encoding the contrasts of interest.
#'
#' The result is a design matrix that can be used directly in differential expression
#' analysis with the expression data.
#'
#' @export
pgx.expMatrix <- function(pheno, contr.matrix) {
  ctx <- rownames(contr.matrix)
  ## already an experiment contrast
  if (length(ctx) == nrow(pheno) && all(ctx %in% rownames(pheno))) {
    contr.matrix <- contr.matrix[match(rownames(pheno), ctx), , drop = FALSE]
    return(contr.matrix)
  }
  group.col <- names(which(apply(pheno, 2, function(x) all(ctx %in% x))))[1]
  if (is.null(group.col) || length(group.col) == 0 || is.na(group.col)) {
    message("[pgx.expMatrix] WARNING: could not resolve group column. No exp.matrix\n")
    return(NULL)
  }

  grp <- pheno[, group.col]
  exp.matrix <- contr.matrix[match(grp, ctx), , drop = FALSE]

  rownames(exp.matrix) <- rownames(pheno)

  return(exp.matrix)
}

## -----------------------------------------------------------------------------
## Contrast creation functions
## -----------------------------------------------------------------------------









#' Make direct contrasts from a factor
#'
#' @param Y A factor vector
#' @param ref The reference level to contrast against
#' @param na.rm Remove missing values? Default is TRUE
#' @param warn Show warnings about missing levels? Default is FALSE.
#'
#' @return A contrast matrix
#'
#' @description
#' Creates a contrast matrix directly from a factor by contrasting each level against the reference.
#'
#' @details
#' This function takes a factor Y and generates a contrast matrix with columns encoding comparisons of each factor level versus the reference level specified by ref.
#'
#' It first expands Y into dummy variables using model.matrix.
#' Then it subtracts the reference level column to create +1/-1 contrasts.
#' Columns with all 0s are dropped.
#'
#' Levels missing in Y will create columns of all 0s. Set na.rm=TRUE to remove these.
#' Warnings about missing levels can be shown with warn=TRUE.
#'
#' @export
makeDirectContrasts <- function(Y, ref, na.rm = TRUE) {
  ## check enough levels
  nlevel <- apply(Y, 2, function(y) length(unique(y)))
  if (any(nlevel < 2)) {
    notlevely <- colnames(Y)[which(nlevel < 2)]
    message("warning:: not enough levels: ", notlevely)
  }
  ii <- which(nlevel > 1)
  ii
  if (length(ii) == 0) {
    message("warning:: no valid phenotypes")
    return(NULL)
  }


  Y <- Y[, ii, drop = FALSE]
  ref <- ref[ii]

  ## make contrast
  exp.matrix <- makeDirectContrasts000(Y = Y, ref = ref, na.rm = na.rm, warn = FALSE)
  exp.matrix <- sign(exp.matrix)
  no.vs <- grep("_vs_|_VS_", colnames(exp.matrix), invert = TRUE)
  no.vs
  if (length(no.vs) > 0) {
    colnames(exp.matrix)[no.vs] <- paste0(colnames(exp.matrix)[no.vs], ":Y_vs_N")
  }
  exp.matrix0 <- exp.matrix
  if (all(grepl("_vs_|_VS_", colnames(exp.matrix0)))) {
    exp.matrix0 <- contrastAsLabels(exp.matrix0)
  }
  group <- pgx.getConditions(exp.matrix0)
  if (length(levels(group)) > 0.5 * nrow(exp.matrix)) {
    cat("WARNING:: contrast matrix looks degenerate. consider removing a contrast.\n")
  }

  contr.matrix <- exp.matrix[which(!duplicated(group)), , drop = FALSE]
  rownames(contr.matrix) <- group[which(!duplicated(group))]

  out <- list(contr.matrix = contr.matrix, group = group, exp.matrix = exp.matrix)
  return(out)
}


#' Make direct contrasts from a factor
#'
#' @param Y A factor vector
#' @param ref The reference level to contrast against
#' @param na.rm Remove missing values? Default is TRUE
#' @param warn Show warnings about missing levels? Default is FALSE.
#'
#' @return A contrast matrix
#'
#' @description
#' Creates a contrast matrix directly from a factor by contrasting each level against the reference.
#'
#' @details
#' This function takes a factor Y and generates a contrast matrix with columns encoding comparisons of each factor level versus the reference level specified by ref.
#'
#' It first expands Y into dummy variables using model.matrix.
#' Then it subtracts the reference level column to create +1/-1 contrasts.
#' Columns with all 0s are dropped.
#'
#' Levels missing in Y will create columns of all 0s. Set na.rm=TRUE to remove these.
#' Warnings about missing levels can be shown with warn=TRUE.
#'
#' @export
makeDirectContrasts000 <- function(Y, ref, na.rm = TRUE, warn = FALSE) {
  if (NCOL(Y) == 1) Y <- data.frame(Y = Y)

  ## check
  all <- c("all", "other", "others", "rest")
  full <- c("*", "full")
  has.ref <- rep(NA, ncol(Y))
  for (i in 1:ncol(Y)) has.ref[i] <- (ref[i] %in% Y[, i] || ref[i] %in% c(all, full))
  has.ref
  if (!all(has.ref)) {
    stop("ERROR:: reference ", which(!has.ref), " not in phenotype matrix\n")
    return(NULL)
  }

  contr.matrix <- c()
  if (length(ref) < ncol(Y)) ref <- Matrix::head(rep(ref, 99), ncol(Y))
  ref.pattern <- "wt|contr|ctr|untreat|normal|^neg|ref|^no$|^0$|^0h$|scrambl|none|dmso|vehicle"
  i <- 1
  for (i in 1:ncol(Y)) {
    m1 <- NULL
    ref1 <- ref[i]
    x <- as.character(Y[, i])
    x[is.na(x) | x == "NA"] <- "_"
    detect.ref <- any(grepl(ref.pattern, x, ignore.case = TRUE))

    if (is.na(ref1) && detect.ref) {
      ref1 <- grep(ref.pattern, x, ignore.case = TRUE, value = TRUE)
      ref1 <- sort(ref1)[1]
      cat("reference auto-detected:", ref1, "\n")
    }
    cref <- as.character(ref1)
    m1 <- stats::model.matrix(~ 0 + x)
    colnames(m1) <- sub("^x", "", colnames(m1))
    if (ref1 %in% full) {
      levels <- names(table(x))
      levels <- setdiff(levels, c(NA, "NA"))

      if (length(levels) > 1) {
        cc <- makeFullContrasts(levels)
        m1 <- m1[, rownames(cc)] %*% cc
      }
    } else if (!is.na(ref1) && !(ref1 %in% all)) {
      m1 <- m1 - m1[, cref] ## +1/-1 encoding
      m1 <- m1[, which(colnames(m1) != cref), drop = FALSE] ## remove refvsref...
      m1 <- m1[, !colnames(m1) %in% c("NA", "_"), drop = FALSE]
      colnames(m1) <- paste0(colnames(m1), "_vs_", ref1)
    } else if (!is.na(ref1) && (ref1 %in% all)) {
      m1 <- t(t(m1 == 1) / Matrix::colSums(m1 == 1) - t(m1 == 0) / Matrix::colSums(m1 == 0))
      m1 <- m1[, !colnames(m1) %in% c("NA", "_"), drop = FALSE]
      colnames(m1) <- paste0(colnames(m1), "_vs_others")
    } else {
      stop("[makeDirectContrasts000] FATAL")
    }
    if (!is.null(m1)) {
      mm <- gsub("[: ]", "_", colnames(Y)[i])
      colnames(m1) <- paste0(mm, ":", colnames(m1))
      contr.matrix <- cbind(contr.matrix, m1)
    }
  }

  ## take out any empty comparisons
  contr.matrix <- contr.matrix[, which(Matrix::colSums(contr.matrix != 0) > 0), drop = FALSE]

  ## normalize to zero mean and symmetric sum-to-one. Any NA to zero.
  for (i in 1:ncol(contr.matrix)) {
    m <- contr.matrix[, i]
    m[is.na(m)] <- 0
    contr.matrix[, i] <- 1 * (m > 0) / sum(m > 0) - 1 * (m < 0) / sum(m < 0)
  }
  rownames(contr.matrix) <- rownames(Y)
  return(sign(contr.matrix))
}


#' Generate contrasts between clusters
#'
#' @param clusters Vector of cluster assignments
#' @param min.freq Minimum cluster frequency to include in contrasts
#' @param full Logical for whether to generate all pairwise contrasts
#' @param by.sample Logical for sample-wise contrasts instead of cluster-wise
#'
#' @return Contrast matrix
#'
#' @description
#' Generates a contrast matrix for comparing clusters.
#'
#' @details
#' This function takes a vector of cluster assignments and generates a contrast matrix
#' encoding comparisons between each pair of clusters above the min.freq threshold.
#'
#' The contrast matrix has clusters as rows and comparison columns with names like
#' "cluster1_vs_cluster2". The matrix values are 1/-1 for the clusters being compared,
#' and 0 otherwise.
#'
#' If full=TRUE, all pairwise comparisons are generated instead of only adjacent clusters.
#' If by.sample=TRUE, the contrasts are encoded at the sample level instead of cluster level.
#'
#' @export
makeFullContrasts <- function(labels, by.sample = FALSE) {
  levels <- sort(unique(as.character(labels)))
  cc <- t(utils::combn(levels, 2))
  contr.matrix <- c()
  for (i in nrow(cc):1) {
    ctr <- 1 * (levels == cc[i, 1]) - 1 * (levels == cc[i, 2])
    contr.matrix <- cbind(ctr, contr.matrix)
    colnames(contr.matrix)[1] <- paste(cc[i, ], collapse = "_vs_")
  }
  rownames(contr.matrix) <- levels
  if (by.sample) {
    design <- stats::model.matrix(~ 0 + labels)
    colnames(design) <- sub("^labels", "", colnames(design))
    rownames(design) <- names(labels)
    design <- design[, rownames(contr.matrix)]
    contr.matrix <- design %*% contr.matrix
  }
  return(contr.matrix)
}





#' Automatically generate contrasts within strata
#'
#' @param df Data frame with sample metadata
#' @param strata.var Column name for the stratification variable
#' @param mingrp Minimum group size to include
#' @param slen Maximum contrast name length
#' @param ref Reference level for factors
#' @param fix.degenerate Fix degenerate designs?
#' @param skip.hidden Skip hidden metadata columns?
#'
#' @return Data frame with automatically generated contrasts
#'
#' @description
#' This function automatically generates contrasts for differential analysis from a
#' data frame of sample metadata, stratified by a grouping variable. It extracts factor
#' columns within each stratum, generates all pairwise comparisons between groups, and
#' outputs a clean contrast matrix.
#'
#' @details
#' The sample metadata is split by the stratification variable into separate strata.
#' Categorical variables with at least mingrp groups are expanded into dummy variables
#' within each stratum. All pairwise comparisons between groups are generated, with
#' contrast names cleaned and truncated to slen characters.
#'
#' Degenerate designs are checked for and resolved if fix.degenerate=TRUE.
#' Metadata columns starting with "." are skipped if skip.hidden=TRUE.
#'
#' The final output is a data frame with samples in rows, contrast columns, and values
#' of 0/1/-1 indicating the sample groups being compared within each stratum.
#'
#' @export
pgx.makeAutoContrastsStratified <- function(df, strata.var, mingrp = 3, slen = 20, ref = NULL,
                                            fix.degenerate = FALSE, skip.hidden = TRUE) {
  df1 <- df[, -match(strata.var, colnames(df)), drop = FALSE]
  strata <- df[, strata.var]
  strata.levels <- unique(strata)
  s <- strata.levels[1]
  ct.all <- NULL
  for (s in strata.levels) {
    sel <- which(strata == s)
    if (length(sel) < mingrp) next
    suppressWarnings(
      ct1 <- pgx.makeAutoContrasts(df1[sel, , drop = FALSE],
        mingrp = mingrp, slen = slen,
        ref = ref, fix.degenerate = fix.degenerate,
        skip.hidden = skip.hidden
      )
    )
    if (is.null(ct1)) next
    ct1x <- ct1$exp.matrix
    colnames(ct1x) <- paste0(colnames(ct1x), "@", s)
    ss <- rownames(df1)[sel]
    if (is.null(ct.all)) {
      ct.all <- data.frame(sample = ss, ct1x, check.names = FALSE)
    } else {
      df2 <- data.frame(sample = ss, ct1x, check.names = FALSE)
      ct.all <- dplyr::full_join(ct.all, df2, by = "sample")
    }
  }

  if (is.null(ct.all)) {
    message("[pgx.makeAutoContrastsStratified] WARNING : no valid contrasts")
    return(NULL)
  }

  ## sample-wise contrasts
  rownames(ct.all) <- ct.all[, "sample"]
  ct.all <- as.matrix(ct.all[, -1, drop = FALSE]) ## drop sample column
  ct.all <- ct.all[match(rownames(df), rownames(ct.all)), ]
  rownames(ct.all) <- rownames(df)
  ct.all[is.na(ct.all)] <- 0

  return(ct.all)
}


#' Automatically generate contrasts from sample metadata
#'
#' @param df Data frame with sample metadata
#' @param mingrp Minimum group size to include
#' @param slen Maximum contrast name length
#' @param ref Reference level for factors
#' @param fix.degenerate Fix degenerate designs?
#' @param skip.hidden Skip hidden metadata columns?
#'
#' @return A data frame with automatically generated contrasts
#'
#' @description
#' This function automatically generates contrasts for differential analysis from a
#' data frame of sample metadata. It extracts factor columns, generates comparisons between
#' groups, handles dummy variable creation, and outputs a clean contrast matrix.
#'
#' @details
#' The sample metadata is scanned to find categorical variables with at least mingrp
#' groups. These are expanded into dummy variables and all pairwise comparisons between
#' groups are generated. Contrast names are cleaned and truncated to slen characters.
#'
#' Degenerate designs are checked for and resolved if fix.degenerate=TRUE.
#' Metadata columns starting with "." are skipped if skip.hidden=TRUE.
#'
#' The final output is a data frame with samples in rows, contrast columns, and values
#' of 0/1/-1 indicating the sample groups being compared. This defines the contrasts
#' for differential analysis.
#'
#' @export
pgx.makeAutoContrasts <- function(df, mingrp = 3, slen = 20, ref = NULL,
                                  fix.degenerate = FALSE, skip.hidden = TRUE) {
  shortestunique <- function(xx, slen = 3) {
    dup <- sapply(
      1:max(nchar(xx)),
      function(i) any(duplicated(substring(xx, 1, i)))
    )
    if (!any(!dup)) {
      return(xx)
    }
    k <- min(which(!dup))
    substring(xx, 1, max(k, slen))
  }

  autoContrast1 <- function(x, ref1, slen, mingrp) {
    ## Automatically create contrast. If 2 levels, create A-vs-B,
    ## otherwise create A-vs-others.
    if (is.null(ref1)) ref1 <- NA
    x <- as.character(x)
    x <- iconv(x, "latin1", "ASCII", sub = "")
    nx <- table(x)
    too.small <- names(which(nx < mingrp))
    if (length(too.small)) x[which(x %in% too.small)] <- NA
    nx <- table(x)
    if (length(nx) < 2) {
      return(NULL)
    }
    x <- factor(x)
    if (!is.na(ref1)) x <- stats::relevel(x, ref = ref1)

    xlevels <- gsub("[^[:alnum:]+-]", "", levels(x))
    levels(x) <- shortestunique(xlevels, slen = slen)
    xref <- gsub("[^[:alnum:]]", "", levels(x)[1])
    nn <- length(nx)
    nn

    if (length(levels(x)) != nn) {
      ## something got wrong...
      return(NULL)
    } else if (nn < 2) {
      return(NULL)
    } else if (nn == 2) {
      ct <- stats::model.matrix(~x)[, 2, drop = FALSE]
      colnames(ct) <- paste0(levels(x)[2], "_vs_", levels(x)[1])
    } else if (nn >= 3) {
      if (is.na(ref1)) {
        ct <- stats::model.matrix(~ 0 + x)
        colnames(ct) <- paste0(levels(x), "_vs_other")
      } else {
        ct <- stats::model.matrix(~ 0 + x)
        colnames(ct) <- paste0(levels(x), "_vs_", xref)
        i <- 1
        for (i in 1:ncol(ct)) {
          j <- which(!(x %in% levels(x)[c(1, i)]))
          j <- intersect(as.character(j), rownames(ct))
          ct[j, i] <- NA
        }
        ct <- ct[, 2:ncol(ct), drop = FALSE] ## remove REFvsREF
      }
    }
    ## NA can make ct smaller than full
    ct <- ct[match(1:length(x), rownames(ct)), , drop = FALSE]
    rownames(ct) <- 1:length(x)
    ct
  }

  ## repeat ref if too short
  if (!is.null(ref) && length(ref) < ncol(df)) ref <- Matrix::head(rep(ref, 99), ncol(df))

  ## filter out 'internal/hidden' and 'group' parameters
  not.used <- grepl("^[.]", colnames(df))
  if (skip.hidden && sum(not.used) > 0 && sum(!not.used) > 0) {
    df <- df[, !not.used, drop = FALSE]
  }

  ## first all to characters
  df.rownames <- rownames(df)
  df <- data.frame(df, check.names = FALSE)
  df <- apply(df, 2, as.character)

  ## trim leading/end parts that are equal
  ##  df <- apply(df, 2, trimsame)  ## need rethink!

  ## try detect (fluffy) comment fields (and remove)
  countSpaces <- function(s) {
    sapply(gregexpr(" ", s), function(p) {
      sum(p >= 0)
    })
  }
  justComment <- function(x) {
    x <- iconv(x, "latin1", "ASCII", sub = "")
    (nchar(x) > 50 || countSpaces(x) >= 4)
  }
  is.comment <- sapply(df[1, ], justComment)
  is.comment
  sel <- which(!is.comment)
  sel
  if (length(sel) == 0) {
    cat("WARNING:: could not auto-find variables...\n")
    return(NULL)
  }
  df <- df[, sel, drop = FALSE]
  if (!is.null(ref)) ref <- ref[sel]
  df[df == ""] <- NA
  df[df == " "] <- NA

  ## ----------- careful cleanup
  df <- apply(df, 2, function(s) gsub("[`'\"]","",s))
  df <- apply(df, 2, function(s) trimws(gsub("[ ]+"," ",s)))
  
  ## ----------- use type.convert to infer parameter type
  df <- utils::type.convert(data.frame(df, check.names = FALSE), as.is = TRUE)

  ## ----------- convert numeric variables into bins
  ii <- which(sapply(df, class) %in% c("integer", "numeric"))
  if (length(ii)) {
    for (i in ii) {
      x <- as.numeric(as.character(df[, i]))
      x <- c("low", "high")[1 + 1 * (x > stats::median(x, na.rm = TRUE))]
      df[, i] <- factor(x, levels = c("low", "high"))
    }
  }  
  
  ## ----------- remove phenotype with too many levels
  level.num   <- apply(df,2,function(x) length(unique(x[!is.na(x)])))  
  level.ratio <- level.num / apply(df,2,function(x) length(x[!is.na(x)]))
  sel <- ( level.num <= 10 | level.ratio < 0.20 )
  df <- df[, sel, drop=FALSE]
  
  ## emergency bail out...
  if (ncol(df) == 0) {
    return(NULL)
  }
  
  ## For each phenotype parameter we 'automagically' try to create a
  ## contrast
  K <- NULL
  i <- 1
  for (i in 1:ncol(df)) {
    ref1 <- NA
    if (!is.null(ref)) ref1 <- ref[i]
    x <- df[, i]
    too.small <- (x %in% names(which(table(x) < mingrp)))
    x[too.small] <- NA
    x <- iconv(x, "latin1", "ASCII", sub = "") ## ???
    if (!(ref1 %in% x)) ref1 <- NA
    ref.pattern <- "wt|contr|ctr|untreat|normal|^neg|ref|^no$|^0$|^0h$|scrambl|none|dmso|vehicle|low|^lo$|null|zero|^not"
    detect.ref <- any(grepl(ref.pattern, x, ignore.case = TRUE))
    if (is.na(ref1) & detect.ref) {
      ref1 <- grep(ref.pattern, x, ignore.case = TRUE, value = TRUE)
      ref1 <- sort(ref1)[1]
      cat("reference auto-detected:", ref1, "\n")
    }
    ct <- autoContrast1(x, ref = ref1, slen = slen, mingrp = mingrp)
    if (!is.null(ct)) {
      colnames(ct) <- paste0(colnames(df)[i], ":", colnames(ct))
      K <- cbind(K, ct)
    }
  }

  if (is.null(K)) {
    warning("[pgx.makeAutoContrasts] non valid contrasts")
    return(NULL)
  }

  rownames(K) <- df.rownames

  ## Now try to infer the underlying "conditions"
  K1 <- contrastAsLabels(K - 0.5)
  kcode <- apply(K1, 1, paste, collapse = "_")
  xc <- factor(kcode, levels = unique(kcode)) ## experimental condition
  if (ncol(K1) > 10) levels(xc) <- paste0("condition", 1:length(levels(xc))) ## too long...
  jj <- which(!duplicated(kcode))
  K2 <- K[jj, colnames(K1), drop = FALSE]
  rownames(K2) <- xc[jj]
  is.degenerate <- (length(jj) > 0.9 * nrow(K1) || mean(table(xc) == 1) > 0.5)
  is.degenerate

  ## THIS IS EXPERIMENTAL: remove
  if (fix.degenerate && is.degenerate) {
    cat("WARNING:: contrast matrix looks degenerate. trying to remove some contrasts...\n")
    is.degenerate <- TRUE
    iter <- 0
    while (is.degenerate && iter < 100) {
      ## Now try to infer the underlying "conditions"
      kcode <- apply(K1, 1, paste, collapse = "_")
      xc <- factor(kcode, levels = unique(kcode)) ## experimental condition
      if (ncol(K1) > 10) levels(xc) <- paste0("condition", 1:length(levels(xc))) ## too long...
      jj <- which(!duplicated(kcode))

      ## SPECIAL CASE!!! if comparisons are degenerate (no valid
      ## condition groups). LIMMA does not like that. Then delete
      ## phenotype with most levels one by one
      is.degenerate <- (length(jj) == nrow(K1) || mean(table(xc) == 1) > 0.5)
      is.degenerate
      if (is.degenerate) {
        ptype <- sub("[:].*", "", colnames(K1))
        del.ptype <- names(which.max(table(ptype)))
        del.ptype
        del <- which(ptype %in% del.ptype)
        K1 <- K1[, -del, drop = FALSE]
      }
      iter <- iter + 1
    }
    iter
    K2 <- K[jj, colnames(K1), drop = FALSE]
    rownames(K2) <- xc[jj]
  } else if (!fix.degenerate && is.degenerate) {
    cat("WARNING:: contrast matrix looks degenerate. going for NULL design...\n")
    ## Go for zero design (no-replicates)
    K2 <- NULL
  }

  ## Translate coding 0/NA/1 to -1/0/+1 coding of contrast
  K[K == 0] <- -1
  K[is.na(K)] <- 0
  if (!is.null(K2)) {
    K2[K2 == 0] <- -1
    K2[is.na(K2)] <- 0
  }
  rownames(K) <- df.rownames

  out <- list(group = xc, contr.matrix = K2, exp.matrix = K)

  return(out)
}


#' Convert contrast matrix to group labels
#'
#' @title Convert contrast matrix to group labels
#'
#' @description Converts a contrast matrix to a data frame of group labels.
#'
#' @param contr.matrix The contrast matrix.
#' @param as.factor Whether to return as factor. Default is FALSE.
#'
#' @details This function takes a contrast matrix and returns a data frame
#' with a column of group labels for each contrast. The labels are generated
#' from the contrast names.
#'
#' @return Data frame of group labels.
#'
#' @examples
#' \dontrun{
#' contrast <- playbase::CONTRASTS
#' z <- playbase::contrastAsLabels(contrast)
#' z
#' }
#' @export
contrastAsLabels <- function(contr.matrix, as.factor = FALSE) {
  contrastAsLabels.col <- function(contr, contr.name) {
    grp1 <- gsub(".*[:]|_vs_.*", "", contr.name)
    grp0 <- gsub(".*_vs_|@.*", "", contr.name)
    x <- rep(NA, length(contr))
    x[which(contr < 0)] <- grp0
    x[which(contr > 0)] <- grp1
    if (as.factor) x <- factor(x, levels = c(grp0, grp1))
    x
  }
  K <- data.frame(contr.matrix[, 0])
  i <- 1
  for (i in 1:ncol(contr.matrix)) {
    contr <- contr.matrix[, i]
    contr.name <- colnames(contr.matrix)[i]
    k1 <- contrastAsLabels.col(contr, contr.name)
    K <- cbind(K, k1)
  }
  colnames(K) <- colnames(contr.matrix)
  rownames(K) <- rownames(contr.matrix)

  return(K)
}


#' Make contrast matrix from label matrix
#'
#' @title Make contrast matrix from label matrix
#'
#' @param lab.matrix Matrix of sample labels
#'
#' @return Contrast matrix
#'
#' @description Creates a contrast matrix from a matrix of sample labels.
#'
#' @details This function takes a matrix of sample labels as input, where the column
#' names indicate the sample groups being compared (e.g. "Group1_vs_Group2"). It parses
#' the column names to extract the two groups, then constructs a contrast matrix by
#' assigning +1/-1 weights to samples belonging to each group.
#'
#' The resulting contrast matrix has rows corresponding to samples, and columns
#' corresponding to the label matrix column names. This encodes the contrasts between
#' each pair of groups.
#'
#' @export
makeContrastsFromLabelMatrix <- function(lab.matrix) {
  ct.names <- colnames(lab.matrix)
  main.grp <- sapply(strsplit(ct.names, split = "_vs_"), "[", 1)
  ctrl.grp <- sapply(strsplit(ct.names, split = "_vs_"), "[", 2)
  main.grp <- sub(".*:", "", main.grp)
  ctrl.grp <- sub("@.*", "", ctrl.grp)

  contr.mat <- matrix(0, nrow(lab.matrix), ncol(lab.matrix))
  rownames(contr.mat) <- rownames(lab.matrix)
  colnames(contr.mat) <- colnames(lab.matrix)
  for (i in 1:ncol(lab.matrix)) {
    lab1 <- trimws(lab.matrix[, i])
    j1 <- which(lab1 == main.grp[i])
    j0 <- which(lab1 == ctrl.grp[i])
    contr.mat[j1, i] <- +1 / length(j1)
    contr.mat[j0, i] <- -1 / length(j0)
  }

  return(contr.mat)
}
