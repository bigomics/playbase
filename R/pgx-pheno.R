##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

##
## Phenotype (samples.csv, pgx$samples) related functions.
##


#' @title Create a Phenotype Matrix
#'
#' @description This function creates a phenotype matrix from a PGX object.
#'
#' @param pgx A list representing a PGX object containing the data to be analyzed.
#' @param phenotype A character string specifying the name of the phenotype to use for creating the matrix.
#'
#' @details This function takes a PGX object `pgx` and a character string `phenotype` as input and
#' creates a phenotype matrix using the specified phenotype. The phenotype matrix is created by extracting
#' the specified phenotype from the `samples` data frame in the `pgx` object and converting it into a model matrix.
#' The resulting phenotype matrix is returned as a numeric matrix.
#'
#' @return A numeric matrix representing the phenotype matrix created from the specified phenotype.
#'
#' @export
pgx.phenoMatrix <- function(pgx, phenotype) {
  y <- pgx$samples[, phenotype]
  mm <- t(stats::model.matrix(~ 0 + y))
  rownames(mm) <- sub("^y", "", rownames(mm))
  colnames(mm) <- rownames(pgx$samples)
  as.matrix(mm)
}






#' @title Check if a Variable is Categorical
#'
#' @description This function checks if a variable is categorical based on its class,
#' number of unique values, and proportion of non-missing values.
#'
#' @param x A vector representing the variable to be checked.
#' @param max.ncat An optional numeric value specifying the maximum number of unique
#' values for a variable to be considered categorical.
#' If not provided, the default value is one less than the length of `x`.
#' @param min.ncat An optional numeric value specifying the minimum number of unique
#' values for a variable to be considered categorical. The default value is 2.
#'
#' @details The function takes a vector `x` as input and checks if it is a factor or character vector.
#' If `x` is not a factor or character vector, the function returns `FALSE`.
#' Otherwise, the function calculates the number of unique values in `x` (excluding missing values) and the proportion of non-missing values in `x`.
#' If the number of unique values is greater than 80% of the proportion of non-missing values, the function returns `FALSE`.
#' Otherwise, if the number of unique values is within the range specified by the `min.ncat` and `max.ncat` parameters, the function returns `TRUE`.
#' Otherwise, the function returns `FALSE`.
#'
#' @return A logical value indicating whether the input variable is categorical.
#'
#' @export
is.categorical <- function(x, max.ncat = NULL, min.ncat = 2) {
  max.ncat <- length(x) - 1
  is.factor <- any(class(x) %in% c("factor", "character"))
  n.unique <- length(unique(setdiff(x, NA)))
  n.notna <- length(x[!is.na(x)])
  is.id <- (n.unique > 0.8 * n.notna)
  is.factor2 <- (is.factor & !is.id & n.unique >= min.ncat & n.unique <= max.ncat)
  return(is.factor2)
}


#' Discretize phenotype matrix
#'
#' @param df Data frame containing phenotype data
#' @param min.ncat Minimum number of categories for discretization. Default is 2.
#' @param max.ncat Maximum number of categories for discretization. Default is 20.
#' @param remove.dup Logical indicating whether to remove duplicate factor levels. Default is FALSE.
#'
#' @return Discretized phenotype data frame
#'
#' @details This function discretizes the phenotype data frame \code{df} into categorical variables.
#' Numerical variables are converted into "high" and "low" values based on the median.
#' Categorical variables are kept as factors.
#'
#' The number of categories is controlled by \code{min.ncat} and \code{max.ncat} parameters.
#' Duplicate factor levels can be removed by setting \code{remove.dup = TRUE}.
#'
#' @export
pgx.discretizePhenotypeMatrix <- function(df, min.ncat = 2, max.ncat = 20, remove.dup = FALSE) {
  df <- as.data.frame(df)
  catpheno <- pgx.getCategoricalPhenotypes(
    df,
    max.ncat = max.ncat, min.ncat = min.ncat, remove.dup = remove.dup
  )
  numpheno <- pgx.getNumericalPhenotypes(df)

  numpheno <- setdiff(numpheno, catpheno) ## already in categories?
  df.num <- c()
  if (length(numpheno)) {
    df.num <- as.matrix(df[, numpheno, drop = FALSE])
    is.high <- t(t(df.num) > apply(df.num, 2, stats::median, na.rm = TRUE))
    df.num[is.high] <- "high"
    df.num[!is.high] <- "low"
  }
  df1 <- df[, 0]
  if (length(catpheno)) df1 <- cbind(df1, df[, catpheno, drop = FALSE])
  if (length(numpheno)) df1 <- cbind(df1, df.num)
  rownames(df1) <- rownames(df)
  df1
}


#' Binarize numerical columns in phenotype matrix
#'
#' @param S Data frame containing phenotype data
#'
#' @return Discretized phenotype data frame
#'
#' @details This function binarizes numerical columns of phenotype
#'   data frame \code{df} into categorical variables.  Numerical
#'   variables are converted into "high" and "low" values based on the
#'   median.  Categorical variables are kept as factors.
#'
#' The number of categories is controlled by \code{min.ncat} and \code{max.ncat} parameters.
#' Duplicate factor levels can be removed by setting \code{remove.dup = TRUE}.
#'
#' @export
binarizeNumericalColumns <- function(S) {
  S <- as.data.frame(S)
  typedS <- type.convert(S, as.is = TRUE)
  is.num <- sapply(typedS, class) %in% c("integer", "numeric")
  numlevels <- apply(S, 2, function(s) length(unique(s[!is.na(s)])))
  is.cat <- (!is.num | numlevels <= 3)
  which.num <- which(is.num & !is.cat)
  which.num
  i <- which.num[1]
  for (i in which.num) {
    mx <- round(median(typedS[, i], na.rm = TRUE), digits = 2)
    grp <- paste0(c(">=", "<"), mx)[1 + 1 * (typedS[, i] >= mx)]
    S[, i] <- grp
  }
  S
}


#' @title Get Numerical Phenotypes
#'
#' @description This function retrieves the column names of numerical
#' phenotypes from a data frame.
#'
#' @param df A data frame containing the data to be analyzed.
#'
#' @details The function takes a data frame `df` as input and searches for columns that
#' represent numerical phenotypes.
#' Numerical phenotypes are identified as columns that are numeric, have a majority of
#' unique values, and do not contain certain keywords in their names or values (e.g., "sample",
#' "id", "replicate", "patient", "donor", "individual", "year", "month", "day", "survival").
#' The names of the identified columns are returned as a character vector.
#'
#' @return A character vector representing the column names of the numerical phenotypes in the input data frame.
#'
#' @export
pgx.getNumericalPhenotypes <- function(df) {
  is.bad <- 0
  is.bad1 <- grepl("^sample$|[_.]id$|replic|rep|patient|donor|individ", tolower(colnames(df)))
  is.bad2 <- grepl("year|month|day|^efs|^dfs|surv|follow", tolower(colnames(df)))
  is.bad3 <- apply(df, 2, function(x) any(grepl("^sample|patient|replicate|donor|individ", x, ignore.case = TRUE)))
  is.bad <- (is.bad1 | is.bad2 | is.bad3)


  numratio <- apply(df, 2, function(x) length(unique(x))) / nrow(df)

  numpheno <- (apply(df, 2, is.num) & !is.bad & numratio > 0.5)

  names(which(numpheno == TRUE))
}


#' @title Get Categorical Phenotypes
#'
#' @description
#' Identifies categorical phenotype columns in a data frame.
#'
#' @param df Data frame containing phenotype data.
#' @param min.ncat Minimum number of categories required. Default 2.
#' @param max.ncat Maximum number of categories allowed. Default 20.
#' @param remove.dup Logical indicating whether to remove duplicate factor
#' levels. Default FALSE.
#'
#' @details
#' This function examines a phenotype data frame \code{df} and identifies
#' columns that represent categorical variables.
#' It first excludes any columns that appear to contain sample IDs.
#' It then selects columns that are factors, or character/numeric with
#' between \code{min.ncat} and \code{max.ncat} unique values.
#'
#' The function returns a character vector of column names representing putative categorical phenotypes.
#' Setting \code{remove.dup = TRUE} will remove any duplicate factor levels prior to counting number of categories.
#'
#' @return
#' A character vector of categorical phenotype column names
#'
#' @export
pgx.getCategoricalPhenotypes <- function(df, min.ncat = 2, max.ncat = 20, remove.dup = FALSE) {
  is.bad <- 0

  ## ... exclude sample IDs

  is.bad1 <- grepl("^sample$|[_.]id$|patient|donor|individ", tolower(colnames(df)))

  ## ... exclude numerical dates/age/year
  is.bad2 <- grepl("ratio|year|month|day|^age$|^efs|^dfs|surv|follow", tolower(colnames(df)))

  is.num <- sapply(df, class) == "numeric"
  is.bad2 <- (is.bad2 & is.num) ## no numeric

  ## ... exclude any sample ID coded in columns...
  is.bad3 <- apply(df, 2, function(x) {
    mean(grepl("^sample|patient|donor|individ",
      x[!is.na(x)],
      ignore.case = TRUE
    )) > 0.8
  })

  is.bad <- (is.bad2 | is.bad3)



  n.unique <- apply(df, 2, function(x) length(unique(setdiff(x, c(NA, "NA", "")))))
  n.notna <- apply(df, 2, function(x) length(x[!is.na(x)]))

  is.factor2 <- (!is.bad & n.unique >= min.ncat & n.unique <= max.ncat)
  is.factor2

  ## take reduced matrix
  df1 <- df[, which(is.factor2), drop = FALSE]
  nlevel <- apply(df1, 2, function(x) length(unique(x)))
  nchars <- apply(df1, 2, function(x) max(nchar(iconv(x, "latin1", "ASCII", sub = ""))))
  df1 <- df1[, order(nlevel, -nchars), drop = FALSE]


  if (remove.dup && ncol(df1) > 1) {
    i <- 1
    j <- 2
    is.dup <- rep(FALSE, ncol(df1))
    for (i in 1:(ncol(df1) - 1)) {
      is.dup[i] <- FALSE
      for (j in (i + 1):ncol(df1)) {
        is.dup[i] <- is.dup[i] || all(rowSums(table(df1[, i], df1[, j]) != 0) == 1)
      }
    }
    df1 <- df1[, which(!is.dup), drop = FALSE]
  }
  return(colnames(df1))
}



#' @title Get Levels of Group Variables
#'
#' @description This function retrieves the levels of group variables in a data frame.
#'
#' @param Y A data frame containing the data to be analyzed.
#'
#' @details The function takes a data frame `Y` as input and searches for
#' columns that represent group variables. Group variables are identified as
#' columns that have a name that does not contain "title", "name", "sample", or
#' "patient" and have a majority of non-unique values.
#' Numeric columns are excluded from the search.
#' The levels of the identified group variables are then extracted and returned as a character vector.
#'
#' @return A character vector representing the levels of the group variables in the input data frame.
#'
#' @export
getLevels <- function(Y) {
  yy <- Y[, grep("title|name|sample|patient", colnames(Y), invert = TRUE), drop = FALSE] ## NEED RETHINK!!!!
  is.grpvar <- apply(yy, 2, function(y) length(table(y)) >= 2)
  is.numeric <- apply(yy, 2, function(y) is.numeric(type.convert(y, as.is = TRUE)))
  is.grpvar <- is.grpvar & !is.numeric
  yy <- yy[, is.grpvar, drop = FALSE]


  levels <- lapply(1:ncol(yy), function(i) unique(paste0(colnames(yy)[i], "=", yy[, i])))
  levels <- sort(unlist(levels))
  return(levels)
}


#' @title Select samples from selected levels
#'
#' @description Selects rows from a data frame corresponding to selected factor levels.
#'
#' @param Y Data frame containing sample data.
#' @param levels Character vector of selected factor level descriptions in format 'factor=level'.
#'
#' @details This function takes a data frame Y containing sample data, and a character vector
#' levels defining a subset of factor levels to select. It parses the level descriptions,
#' identifies matching samples in Y, and returns the rownames of selected samples.
#'
#' The levels should be strings in format 'factor=level' defining the
#' factor name and level value.  Samples in Y matching any of the
#' provided levels will be selected. Note: across factors samples as
#' intersected ('and' operator), within a factor samples are combined
#' ('or' operator). Not very logical but this is what people in practice want.
#'
#' @return A character vector of selected sample names.
#'
#' @export
selectSamplesFromSelectedLevels <- function(Y, levels) {
  if (is.null(levels) || all(levels == "")) {
    return(rownames(Y))
  }

  # fill ="" will (unfortunately) still return NA when level is "NA"... which crashes when phenotype is ""
  pheno <- data.table::tstrsplit(levels, "=", keep = 1) %>%
    unlist()
  ptype <- data.table::tstrsplit(levels, "=", keep = 2) %>%
    unlist()
  # force replace "NA" by NA
  ptype[ptype == "NA"] <- NA
  ##  sel <- rep(FALSE, nrow(Y))
  sel <- rep(TRUE, nrow(Y))
  for (ph in unique(pheno)) {
    # ph = pheno[1]
    k <- which(pheno == ph)
    ##    sel <- sel | (Y[, ph] %in% ptype[k])
    sel <- sel & (Y[, ph] %in% ptype[k])
  }

  return(rownames(Y)[which(sel)])
}


#' @title Expand annotation matrix
#'
#' @description Expands a phenotype annotation matrix into dummy variables.
#'
#' @param A Data frame containing the annotation variables.
#'
#' @details This function takes an annotation data frame and expands any categorical variables into
#' dummy variables using model.matrix.
#'
#' For each column, it determines if the variable is numeric or categorical.
#' Numeric variables are ranked. Categorical variables are expanded into dummy variables.
#'
#' @return An expanded annotation matrix with dummy variables.
#'
#' @export
expandAnnotationMatrix <- function(A) {
  expandPhenoMatrix(A)
}


#' @title Expand phenotype matrix
#'
#' @description Expands a phenotype data matrix into dummy variables while optionally
#' dropping the reference level.
#'
#' @param pheno Data frame containing the phenotype variables.
#' @param drop.ref Logical indicating whether to drop the reference level for each factor.
#'
#' @details This function takes a phenotype data matrix and expands any categorical variables into
#' dummy variables
#'
#' For each column, it determines if the variable is numeric or categorical. Numeric variables are
#' dichotomized. Categorical variables are expanded into dummy variables using \code{model.matrix}.
#'
#' @return An expanded phenotype matrix with dummy variables suitable for regression modeling.
#' @export
expandPhenoMatrix <- function(M, drop.ref = TRUE, keep.numeric = FALSE, check = TRUE) {
  ## get expanded annotation matrix
  a1 <- tidy.dataframe(M)
  nlevel <- apply(a1, 2, function(x) length(setdiff(unique(x), NA)))
  nterms <- colSums(!is.na(a1))
  nratio <- nlevel / nterms
  if (inherits(a1, "data.frame")) {
    a1.typed <- utils::type.convert(a1, as.is = TRUE)
    y.class <- sapply(a1.typed, function(a) class(a)[1])
  } else {
    ## matrix??
    a1.typed <- utils::type.convert(a1, as.is = TRUE)
    y.class <- apply(a1.typed, 2, function(a) class(a)[1])
  }

  ## these integers are probably factors... (mostly...)
  is.fac <- rep(FALSE, ncol(a1))
  is.int <- (y.class == "integer")
  ii <- which(is.int)
  is.fac[ii] <- apply(a1[, ii, drop = FALSE], 2, function(x) {
    nlev <- length(unique(x[!is.na(x)]))
    max(x, na.rm = TRUE) %in% c(nlev, nlev - 1) ## boolean
  })
  is.fac2 <- (y.class == "integer" & nlevel <= 3 & nratio < 0.66)
  y.class[is.fac | is.fac2] <- "character"

  ## select allowed columns: numeric or with "sensible" levels
  y.isnum <- (y.class %in% c("numeric", "integer"))
  kk <- which(y.isnum | (!y.isnum & nlevel > 1 & nratio < 0.66))
  if (length(kk) == 0) {
    return(NULL)
  }
  a1 <- a1[, kk, drop = FALSE]
  a1.isnum <- y.isnum[kk]

  i <- 1
  m1 <- list()
  for (i in 1:ncol(a1)) {
    if (a1.isnum[i]) {
      suppressWarnings(x <- as.numeric(a1[, i]))
      if (keep.numeric) {
        m0 <- matrix(x, ncol = 1)
        colnames(m0) <- "#"
      } else {
        if (drop.ref) {
          m0 <- matrix((x > stats::median(x, na.rm = TRUE)), ncol = 1)
          colnames(m0) <- "high"
        } else {
          mx <- stats::median(x, na.rm = TRUE)
          m0 <- matrix(cbind(x <= mx, x > mx), ncol = 2)
          colnames(m0) <- c("low", "high")
        }
      }
    } else if (drop.ref && nlevel[i] == 2) {
      x <- as.character(a1[, i])
      x1 <- Matrix::tail(sort(x), 1)
      m0 <- matrix(x == x1, ncol = 1)
      colnames(m0) <- x1
    } else {
      x <- as.character(a1[, i])
      x[is.na(x) | x == "NA" | x == " "] <- "_"
      m0 <- stats::model.matrix(~ 0 + x)
      colnames(m0) <- sub("^x", "", colnames(m0))
    }
    rownames(m0) <- rownames(a1)
    ## remove "_"
    if ("_" %in% colnames(m0)) {
      m0 <- m0[, -which(colnames(m0) == "_")]
    }
    m1[[i]] <- m0
  }

  ## create level names
  names(m1) <- colnames(a1)
  for (i in 1:length(m1)) {
    colnames(m1[[i]]) <- paste0(names(m1)[i], "=", colnames(m1[[i]]))
  }
  m1 <- do.call(cbind, m1)
  colnames(m1) <- sub("=#", "", colnames(m1))
  rownames(m1) <- rownames(M)
  return(m1)
}


## =================================================================================
## ========================= END OF FILE ===========================================
## =================================================================================
