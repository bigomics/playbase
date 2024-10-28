##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


## ----------------------------------------------------------------------
## Partial correlation functions
## ----------------------------------------------------------------------


#' @title Compute partial correlations around a gene
#'
#' @description
#' Estimates partial correlations between a target gene and other
#' genes using glasso.
#'
#' @param X Expression matrix, genes in rows and samples in columns.
#' @param gene Gene name or index for target gene.
#' @param nmax Maximum number of genes to use.
#'
#' @return List containing:
#' \itemize{
#'   \item cor: Correlation matrix
#'   \item pcor: Partial correlation matrix
#'   \item wi: Inverse covariance matrix
#'   \item sampleSize: Number of samples used
#' }
#'
#' @details
#' This function calculates the correlation between a target gene and
#' the top nmax correlated genes.
#'
#' It uses the glasso algorithm to estimate a sparse inverse covariance
#' matrix and derive partial correlations.
#'
#' The correlation matrix, partial correlation matrix, inverse covariance
#' matrix (precision matrix), and number of samples used are returned.
#'
#' @export
pgx.computeGlassoAroundGene <- function(X, gene, nmax = 100) {
  rho <- stats::cor(t(X), t(X[gene, , drop = FALSE]))
  jj <- Matrix::head(order(-rowMeans(rho**2, na.rm = TRUE)), nmax)
  tX <- t(X[jj, ])

  vX <- stats::var(tX)
  res <- glasso::glasso(vX, 0.1)
  res$cor <- Matrix::cov2cor(vX)
  res$sampleSize <- nrow(tX)

  res$pcor <- -stats::cov2cor(res$wi)
  diag(res$pcor) <- 1
  rownames(res$pcor) <- rownames(res$cor)
  colnames(res$pcor) <- colnames(res$cor)

  return(res)
}


#' Compute partial correlation around a gene
#'
#' @title Compute partial correlation around a gene
#'
#' @description Computes partial correlation between a given gene and all other genes.
#'
#' @param pgx pgx object
#' @param gene Character string specifying the gene name to compute partial correlation around.
#' @param rho.min Numeric threshold for including genes in partial correlation graph.
#' @param pcor.min Numeric threshold for including genes in partial correlation graph.
#' @param nsize Integer specifying number of top correlated genes to include.
#' @param main Character string for plot main title.
#' @param what Character vector specifying which plots to generate. Options are
#' "cor", "pcor", "graph".
#' @param edge.width Numeric line width for graph edges.
#' @param layout Layout algorithm for graph. Options are "fr", "kk", "lgl", etc.
#'
#' @details This function computes the partial correlation between a given
#' gene and all other genes in the dataset.
#' It returns plots of the top correlated genes, a correlation graph, and
#' optionally a matrix of all gene-gene correlations.
#'
#' @return List with components:
#' \itemize{
#'  \item rho: List of partial correlation matrices
#'  \item meta.pcor: Matrix of mean partial correlations
#'  \item timings: Matrix of timings
#' }
#'
#' @export
pgx.plotPartialCorrelationGraph <- function(res, gene, rho.min = 0.1, nsize = -1, main = "",
                                            vsize = 10, edge.width = 10, label.cex = 0.8,
                                            radius = -1, plot = TRUE, layout = "fr") {
  ## GLASSO object
  nn <- rownames(res$cor)
  R <- stats::cov2cor(res$w)
  diag(R) <- 0
  utils::tail(sort(abs(R)), 20)
  P <- stats::cov2cor(res$wi)
  diag(P) <- 0
  utils::tail(sort(abs(P)), 20)
  rownames(R) <- colnames(R) <- nn
  rownames(P) <- colnames(P) <- nn


  G <- igraph::graph_from_adjacency_matrix(
    abs(P),
    mode = "undirected", diag = FALSE, weighted = TRUE
  )

  ee <- igraph::get.edges(G, igraph::E(G))
  ew <- (abs(P[ee]) / max(abs(P[ee])))
  igraph::E(G)$rho <- P[ee]
  igraph::V(G)$name
  igraph::V(G)$size <- vsize * c(1, 1.8)[1 + 1 * (igraph::V(G)$name == gene)]
  igraph::V(G)$label.cex <- label.cex * c(1, 1.4)[1 + 1 * (igraph::V(G)$name == gene)]
  igraph::V(G)$color <- c("lightskyblue", "salmon")[1 + 1 * (igraph::V(G)$name == gene)]

  igraph::E(G)$width <- edge.width * ew**1.5
  igraph::E(G)$color <- c("red", "darkgreen")[1 + 1 * (sign(R[ee]) > 0)]
  igraph::E(G)$color <- c("#FF000033", "#44444433")[1 + 1 * (sign(R[ee]) > 0)]
  igraph::E(G)$color <- c(paste0(omics_colors("red"), "33"), paste0(omics_colors("brand_blue"), "33"))[1 + 1 * (sign(R[ee]) > 0)]

  ## delete weak edges
  G1 <- igraph::delete_edges(G, which(abs(igraph::E(G)$weight) < rho.min))
  if (length(igraph::V(G1)) <= 1) {
    return()
  }

  ## delete orphans
  isolated <- which(igraph::degree(G1) == 0 & igraph::V(G)$name != gene)
  G1 <- igraph::delete.vertices(G1, isolated)

  if (length(igraph::V(G1)) <= 1) {
    return()
  }

  ## get only nodes within radius
  dist <- igraph::distances(G1, v = gene, mode = "all", weights = NA)[1, ]
  vv <- igraph::V(G1)
  if (radius > 0) {
    vv <- names(which(dist <= radius))
  }
  G1 <- igraph::induced_subgraph(G1, vv)


  L1 <- igraph::layout_with_fr(G1)

  if (length(igraph::V(G1)) == 0) {
    return()
  }

  if (plot == TRUE) {
    graphics::par(mar = c(1, 1, 1, 1) * 0)
    plot(G1,
      vertex.color = "lightskyblue1",
      vertex.frame.color = "skyblue"
    )
  }

  return(G1)
}


PCOR.FAST <- c("cor", "pcor", "pcor.shrink", "SILGGM", "FastGGM")
PCOR.METHODS <- c(
  "cor", "pcor", "pcor.shrink", "glasso", "huge",
  ## "clime", "BigQuic", "QUIC",
  "fastclime", "FastGGM", "SILGGM"
)


#' Computes partial correlation between a given gene and all other genes in a data matrix X.
#'
#' @param X Numeric gene expression matrix with genes as rows.
#' @param gene Gene name to compute partial correlation around.
#' @param method Partial correlation method(s) to use. Default is PCOR.METHODS.
#' @param nmax Maximum number of genes to return correlations for. Default 100.
#' @param fast Use fast approximation methods if TRUE. Default FALSE.
#'
#' @details
#' This function calculates the partial correlation between a specified gene and
#' all other genes in the expression matrix X.
#' It returns a list containing the partial correlation matrix and p-values.
#'
#' The \code{method} argument specifies which partial correlation algorithm to use.
#' The default \code{PCOR.METHODS} uses both Pearson and Spearman correlation.
#'
#' The number of top genes to return correlations for can be controlled with \code{nmax}.
#' Setting \code{fast=TRUE} will use faster approximations for large matrices.
#'
#' @return
#' A named list with components:
#' \itemize{
#'   \item rho - matrix of partial correlations
#'   \item p.value - matrix of p-values
#' }
#'
#' @export
pgx.computePartialCorrelationAroundGene <- function(X, gene, method = PCOR.METHODS, nmax = 100, fast = FALSE) {
  rho <- stats::cor(t(X), t(X[gene, , drop = FALSE]))
  jj <- Matrix::head(order(-rowMeans(rho**2, na.rm = TRUE)), nmax)
  tX <- t(X[jj, ])
  res <- pgx.computePartialCorrelationMatrix(
    tX,
    method = method, fast = fast
  )
  lapply(res$rho, dim)
  names(res)
  return(res)
}


#' @title Compute partial correlation matrix
#'
#' @description
#' Computes a partial correlation matrix between genes adjusting for covariates.
#'
#' @param pgx A PGX object containing gene expression data.
#' @param covariates Character vector of covariate names to adjust for.
#' @param nsize Number of top connected genes to extract network. Default 0.
#' @param what Character vector specifying output elements. Default c("matrix", "graph").
#' @param gene Character vector of one or more genes to highlight. Default NULL.
#'
#' @details
#' This function calculates the partial correlation between all genes in the PGX expression matrix,
#' adjusting for the specified covariates. It removes the linear effect of the covariates before
#' computing the gene-gene partial correlations.
#'
#' The main output is a partial correlation matrix. Optionally a graph can be extracted, showing the
#' top connected genes. Specific genes can be highlighted in the graph.
#'
#' @return
#' A list with components:
#' \itemize{
#'   \item rho - Matrix of partial correlations
#'   \item graph - igraph object of top connected genes
#' }
#'
#' @export
pgx.computePartialCorrelationMatrix <- function(tX, method = PCOR.METHODS, fast = FALSE) {
  if (fast || ncol(tX) > 1000) {
    method <- intersect(method, PCOR.FAST)
    method
  }
  method <- unique(c("cor", method)) ## always cor
  timings <- list()
  rho <- list()
  timings

  if ("cor" %in% method) {
    timings[["cor"]] <- system.time(
      rho[["cor"]] <- stats::cor(tX)
    )
  }

  if ("pcor" %in% method) {
    timings[["pcor"]] <- system.time(
      suppressWarnings(rho[["pcor"]] <- ppcor::pcor(tX)$estimate)
    )
  }

  if ("pcor.shrink" %in% method) {
    timings[["pcor.shrink"]] <- system.time(
      suppressWarnings(rho[["pcor.shrink"]] <- corpcor::pcor.shrink(tX))
    )
  }

  if ("huge" %in% method) {
    timings[["huge.glasso"]] <- system.time(
      out <- huge(tX, method = "glasso")
    )

    c2 <- -out$icov[[length(out$icov)]]
    diag(c2) <- -diag(c2)
    rho[["huge"]] <- Matrix::cov2cor(c2)
  }

  if ("FastGGM" %in% method) {
    timings[["FastGGM"]] <- system.time(
      # The f FastGGM_Parallel likely comes from wt2015-github/FastGGM but better confirm
      suppressWarnings(res <- FastGGM_Parallel(tX))
    )
    rho[["FastGGM"]] <- res$partialCor
  }

  if ("SILGGM" %in% method) {
    timings[["SILGGM"]] <- system.time({
      ## these are quite fast
      rho[["SILGGM"]] <- SILGGM::SILGGM(tX)$partialCor ## default: "D-S_NW_SL"
    })
  }

  if ("fastclime" %in% method) {
    timings[["fastclime"]] <- system.time(
      # The f fastclime.selector likely comes from fastclime but better confirm
      out1 <- fastclime(tX)
    )
    # The f fastclime.selector likely comes from fastclime but better confirm
    out2 <- fastclime.selector(out1$lambdamtx, out1$icovlist, 0.1)
    out2$adaj
    rho[["fastclime"]] <- out2$adaj
  }

  if ("QUIC" %in% method) {
    timings[["QUIC"]] <- system.time(
      rho[["QUIC"]] <- {
        # The f QUIC likely comes from QUIC but better confirm
        r <- -QUIC(stats::cov(tX), 1e-1)$X
        diag(r) <- -diag(r)
        r
      }
    )
  }

  if ("glasso" %in% method) {
    timings[["glasso"]] <- system.time(
      res <- glasso::glasso(stats::cov(tX), 1e-2)
    )
    r <- -stats::cov2cor(res$wi)
    diag(r) <- -diag(r)
    rho[["glasso"]] <- r
  }

  if ("clime" %in% method) {
    timings[["clime"]] <- system.time(
      rho[["clime"]] <- {
        # The f cv.clime and clime likely comes from fastclime but better confirm
        re.cv <- cv.clime(clime(tX))
        clime(tX, standardize = FALSE, re.cv$lambdaopt)
      }
    )
  }

  if ("BigQuic" %in% method) {
    timings[["BigQuic"]] <- system.time(
      rho[["BigQuic"]] <- {
        lambda <- seq(from = 0.1, to = 1, by = 0.1) - 0.01
        # The f BigQuic and clime likely comes from BigQuic but better confirm
        res <- BigQuic(as.matrix(tX),
          lambda = lambda,
          numthreads = 10, memory_size = 512 * 8,
          use_ram = TRUE
        )
        # The f BigQuic.select and clime likely comes from BigQuic but better confirm
        res <- BigQuic.select(res)
        which.opt <- max(which(res$lambda <= res$opt.lambda))
        res$precision_matrices[[which.opt]]
      }
    )
  }

  for (i in 1:length(rho)) {
    rownames(rho[[i]]) <- colnames(rho[[i]]) <- colnames(tX)
    if (sum(diag(rho[[i]])) == 0) diag(rho[[i]]) <- 1
    rho[[i]] <- Matrix::cov2cor(as.matrix(rho[[i]]))
  }

  ## compute average partical correlation (over all methods)
  nother <- length(setdiff(names(rho), c("cor", "rho")))
  meta.pcor <- NULL
  if (nother > 0) {
    flat.rho <- as.vector(rho[["cor"]])
    flat.prho <- sapply(rho[which(names(rho) != "cor")], as.vector)
    jj <- which(abs(flat.prho) > abs(flat.rho)) ## not reliable
    flat.prho[jj] <- NA
    nx <- ncol(tX)
    meta.pcor <- matrix(apply(flat.prho, 1, mean, na.rm = TRUE), nx, nx)
    rownames(meta.pcor) <- colnames(meta.pcor) <- colnames(tX)
  }

  timings <- do.call(rbind, timings)
  res <- list(rho = rho, meta.pcor = meta.pcor, timings = timings)
  if ("cor" %in% names(rho)) res$cor <- rho[["cor"]]
  names(res)
  return(res)
}






#' Test correlation of phenotype with expression data
#'
#' @param df Data frame containing phenotype data
#' @param plot Logical, whether to generate a correlation heatmap plot
#' @param cex Text size factor for plot labels
#'
#' @return A matrix of correlation p-values
#'
#' @description Test correlation between phenotype variables and gene expression, using
#' appropriate methods based on data types.
#'
#' @details This function takes a data frame \code{df} containing phenotype data, including both
#' continuous and discrete variables. It calculates correlation p-values between the phenotype
#' variables and gene expression data in a PGX object, using either Pearson or ANOVA. Discrete
#' variables are compared between groups using Fisher's exact test. The results are returned as
#' a matrix of p-values, with rows corresponding to expression profiles. If \code{plot=TRUE},
#' a correlation plot of the significance (-log10 p-values) is also generated.
#'
#' @export
pgx.testPhenoCorrelation <- function(df, plot = TRUE, cex = 1, compute.pv = TRUE) {
  cl <- sapply(df, class)
  nlev <- apply(df, 2, function(x) length(unique(x[!is.na(x)])))
  cvar <- which(cl %in% c("numeric", "integer") & nlev >= 2)
  dvar <- which(cl %in% c("factor", "character") & nlev >= 2)
  dc <- df[, cvar, drop = FALSE]
  dd <- df[, dvar, drop = FALSE]

  ## generalized correlation matrix
  ddx <- expandPhenoMatrix(dd, drop.ref = FALSE)
  Rx <- cor(cbind(dc, ddx), use = "pairwise")
  rvar <- sub("=.*", "", colnames(Rx))
  Rx[is.nan(Rx)] <- 0
  Rx[is.na(Rx)] <- 0
  R <- tapply(1:nrow(Rx), rvar, function(i) apply(Rx[c(i, i), , drop = FALSE], 2, max, na.rm = TRUE))
  R <- do.call(rbind, R)
  R <- tapply(1:ncol(R), rvar, function(i) apply(R[, c(i, i), drop = FALSE], 1, max, na.rm = TRUE))
  if (length(R) == 1) {
    R <- matrix(R)
  } else {
    R <- do.call(cbind, R)
  }
  R <- t(R / sqrt(diag(R))) / sqrt(diag(R))
  R[is.nan(R)] <- NA

  P <- Q <- NULL
  if (compute.pv) {
    ## discrete vs discrete -> Fisher test
    fisher.P <- NULL
    if (ncol(dd)) {
      fisher.P <- matrix(NA, ncol(dd), ncol(dd))
      if(nrow(fisher.P) == 1 && ncol(fisher.P) == 1) {
          tb <- table(dd[, 1], dd[, 1])
          fisher.P[1, 1] <- stats::fisher.test(tb, simulate.p.value = TRUE)$p.value
      } else {
          i <- 1
          j <- 2
          for (i in 1:(ncol(dd) - 1)) {
              kk <- which(!is.na(dd[, i]) & !is.na(dd[, j]))
              if (length(unique(dd[kk, i])) < 2 || length(unique(dd[kk, j])) < 2) next
              for (j in (i + 1):ncol(dd)) {
                  tb <- table(dd[, i], dd[, j])
                  fisher.P[i, j] <- stats::fisher.test(tb, simulate.p.value = TRUE)$p.value
              }
          }
      }
      rownames(fisher.P) <- colnames(dd)
      colnames(fisher.P) <- colnames(dd)
    }

    ## discrete vs continuous -> ANOVA or Kruskal-Wallace
    kruskal.P <- NULL
    if (ncol(dc) > 0) {
      kruskal.P <- matrix(NA, ncol(dd), ncol(dc))
      for (i in 1:ncol(dd)) {
        for (j in 1:ncol(dc)) {
          kk <- which(!is.na(dc[, j]) & !is.na(dd[, i]))
          if (length(unique(dd[kk, i])) < 2) next
          kruskal.P[i, j] <- stats::kruskal.test(dc[kk, j], dd[kk, i])$p.value
        }
      }
      rownames(kruskal.P) <- colnames(dd)
      colnames(kruskal.P) <- colnames(dc)
    }

    ## continuous vs continuous -> correlation test
    cor.P <- NULL
    if (ncol(dc) > 1) {
      cor.P <- matrix(NA, ncol(dc), ncol(dc))
      i <- 1
      j <- 2
      for (i in 1:(ncol(dc) - 1)) {
        for (j in (i + 1):ncol(dc)) {
          cor.P[i, j] <- stats::cor.test(dc[, i], dc[, j])$p.value
        }
      }
      rownames(cor.P) <- colnames(dc)
      colnames(cor.P) <- colnames(dc)
    }

    P <- matrix(NA, ncol(df), ncol(df))
    rownames(P) <- colnames(P) <- colnames(df)

    if (!is.null(fisher.P)) {
      ii <- match(rownames(fisher.P), rownames(P))
      jj <- match(colnames(fisher.P), colnames(P))
      P[ii, jj] <- fisher.P
    }

    if (!is.null(kruskal.P)) {
      ii <- match(rownames(kruskal.P), rownames(P))
      jj <- match(colnames(kruskal.P), colnames(P))
      P[ii, jj] <- kruskal.P
    }

    if (!is.null(cor.P)) {
      ii <- match(rownames(cor.P), rownames(P))
      jj <- match(colnames(cor.P), colnames(P))
      P[ii, jj] <- cor.P
    }

    ij <- which(!is.na(P), arr.ind = TRUE)
    qv <- stats::p.adjust(P[ij], method = "BH")
    Q <- P
    Q[ij] <- qv

    P[is.na(P)] <- 0
    P <- (P + t(P)) / 2
    Q[is.na(Q)] <- 0
    Q <- (Q + t(Q)) / 2
  }

  BLUERED <- grDevices::colorRampPalette(c(omics_colors("brand_blue"), "white", omics_colors("red")))

  if (plot == TRUE) {
    if (compute.pv) {
      X <- -log10(Q + 1e-8)
    } else {
      X <- R
    }
    diag(X) <- 0
    X[is.na(X)] <- 0
    corrplot::corrplot(
      corr = X,
      type = "upper",
      col = BLUERED(25),
      ## is.corr = (!compute.pv),
      is.corr = FALSE,
      mar = c(0, 0, 0, 2),
      p.mat = Q,
      sig.level = 0.05,
      tl.cex = cex,
      tl.col = "#3b4252",
      tl.offset = 1,
      cl.align.text = "l",
      cl.offset = 0.25,
      cl.cex = 0.7,
      pch.col = "#3b4252",
      order = "hclust",
      number.digits = 2
    )
  }

  return(list(R = R, Rx = Rx, P = P, Q = Q))
}

#' @export
pgx.testPhenoCorrelation2 <- function(df, plot = TRUE, cex = 1, compute.pv = TRUE) {
  cl <- sapply(df, class)
  nlev <- apply(df, 2, function(x) length(unique(x[!is.na(x)])))
  cvar <- which(cl %in% c("numeric", "integer") & nlev >= 2)
  dvar <- which(cl %in% c("factor", "character") & nlev >= 2)
  dc <- df[, cvar, drop = FALSE]
  dd <- df[, dvar, drop = FALSE]

  ## generalized correlation matrix
  dY <- expandPhenoMatrix(dd, drop.ref = FALSE)
  R <- cor(dY, use = "pairwise")
  R[is.na(R)] <- 0

  P <- NULL
  if (compute.pv) {
    testRes <- corrplot::cor.mtest(dY, conf.level = 0.95)
    P <- testRes$p
  }

  if (plot) {
    BLUERED <- grDevices::colorRampPalette(c("blue3", "white", "red3"))
    corrplot::corrplot(
      corr = R, order = "hclust", type = "lower",
      p.mat = P, insig = "blank", ## sig.level = 0.05,
      tl.col = "black", tl.srt = 90, tl.cex = cex,
      col = rev(corrplot::COL2("RdBu"))
    ) -> corrRes
    if (0) {
      p1 <- corrRes$corrPos
      jj <- which(p1$p.value < 0.05)
      text(p1$x[jj], p1$y[jj], round(p1$corr[jj], 2))
    }
  }

  return(list(R = R, Rx = NULL, P = P, Q = NULL))
}


#' @title Get correlation of a gene with other genes
#'
#' @param gene Character string specifying the target gene name.
#' @param xref List of gene expression matrices to compute correlation against.
#'
#' @return Matrix of correlation values between the target gene and other genes.
#'
#' @description Computes correlation between a given gene and all other genes using
#' external gene expression data.
#'
#' @details This function takes a gene name and a list of gene expression matrices.
#' It calculates the correlation between the target gene and all other genes in
#' each expression data matrix.
#'
#' The function returns a matrix with genes in rows and datasets in columns.
#' Each column contains the correlation values between the target gene and
#' other genes in that dataset.
#'
#' @export
pgx.getGeneCorrelation <- function(gene, xref) {
  rho.genes <- unlist(lapply(xref, rownames))
  rho.genes <- sort(unique(rho.genes))

  R <- NULL
  ## correlation using external datasets
  k <- 1
  for (k in 1:length(xref)) {
    has.gene <- (toupper(gene) %in% toupper(rownames(xref[[k]])))
    if (has.gene) {
      ## take geometric mean with correlation in TISSUE
      xx <- log2(1 + xref[[k]])
      jj <- match(toupper(gene), toupper(rownames(xx)))
      tx <- xx[jj, ]
      if (inherits(xx, "dgCMatrix")) {
        suppressWarnings(rho1 <- qlcMatrix::corSparse(xx, cbind(tx))[, 1])
      } else {
        suppressWarnings(rho1 <- stats::cor(t(xx), tx)[, 1])
      }
      if (sum(is.na(rho1)) > 0) rho1[which(is.na(rho1))] <- 0

      rho1 <- rho1[match(rho.genes, names(rho1))]
      R <- cbind(R, rho1)
      if (NCOL(R) == 1) R <- matrix(R, ncol = 1)
      rownames(R) <- rho.genes
      colnames(R)[ncol(R)] <- names(xref)[k]
    }
  }


  if (!is.null(R) && NCOL(R) > 0) {
    R[is.na(R)] <- 0
    R[is.nan(R)] <- 0
    R <- R[which(rowSums(R != 0, na.rm = TRUE) > 0), , drop = FALSE]
    ## geneset rho has no sign (=cosine correlation) so we use the sign of others
    if (ncol(R) > 1 && "gene sets" %in% colnames(R)) {
      k <- which(colnames(R) == "gene sets")
      R[, "gene sets"] <- R[, "gene sets"] * sign(rowMeans(R[, -k, drop = FALSE], na.rm = TRUE))
    }
  }

  return(R)
}
