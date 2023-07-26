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
  jj <- Matrix::head(order(-rowMeans(rho**2)), nmax)
  tX <- t(X[jj, ])

  vX <- var(tX)
  res <- glasso::glasso(vX, 0.1)
  res$cor <- Matrix::cov2cor(vX)
  res$sampleSize <- nrow(tX)

  res$pcor <- -cov2cor(res$wi)
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
  R <- cov2cor(res$w)
  diag(R) <- 0
  tail(sort(abs(R)), 20)
  P <- cov2cor(res$wi)
  diag(P) <- 0
  tail(sort(abs(P)), 20)
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
    par(mar = c(1, 1, 1, 1) * 0)
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
  jj <- Matrix::head(order(-rowMeans(rho**2)), nmax)
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
        r <- -QUIC(cov(tX), 1e-1)$X
        diag(r) <- -diag(r)
        r
      }
    )
  }

  if ("glasso" %in% method) {
    timings[["glasso"]] <- system.time(
      res <- glasso::glasso(cov(tX), 1e-2)
    )
    r <- -cov2cor(res$wi)
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


#' @title Plot Partial Correlation Network Around a Gene
#'
#' @description 
#' Visualize a partial correlation network for a given gene using the correlation results from pgx.pcor().
#'
#' @param res The correlation analysis result list returned by \code{\link{pgx.pcor}}.
#' @param gene The gene name to extract partial correlations around.  
#' @param rho.min Minimum absolute partial correlation threshold.
#' @param pcor.min Minimum absolute meta-partial correlation threshold.
#' @param nsize Node size parameter.
#' @param main Plot title.
#' @param what Type of plot, "cor", "pcor" or "graph".
#' @param edge.width Edge width for graph.
#' @param layout Layout algorithm for graph plot.
#'
#' @details
#' This function extracts the partial correlation matrix around the specified 
#' \code{gene} from the \code{res} object.
#' It then optionally thresholds the correlations, and visualizes them as a 
#' correlogram or network graph, with the central gene in the middle.
#'
#' For the graph, the \code{rho.min} and \code{pcor.min} arguments filter edges to 
#' plot. Node size is proportional to the number of connections. The graph can be 
#' drawn using different layout algorithms like \code{"fr"} or \code{"circle"}.
#'
#' @return 
#' Either a correlogram or graph plot centered on the given gene is produced.
#'
#' @export
pgx.plotPartialCorrelationAroundGene <- function(res, gene, rho.min = 0.8, pcor.min = 0,
                                                 nsize = -1, main = "",
                                                 what = c("cor", "pcor", "graph"),
                                                 edge.width = 10, layout = "fr") {
  rho <- res$rho
  j <- which(colnames(res$rho[[1]]) %in% gene)
  rho2 <- lapply(res$rho, function(x) x[j, ])
  M <- t(as.matrix(do.call(rbind, rho2)))
  M[is.infinite(M) | is.nan(M)] <- NA
  M <- M[order(-rowMeans(M**2, na.rm = TRUE)), , drop = FALSE]


  ## ------------------------------------------------------------
  ## Correlation barplots
  ## ------------------------------------------------------------

  if ("cor" %in% what) {
    par(mar = c(10, 4, 4, 2))
    r1 <- M[, "cor"] ## / ncol(M)
    if (nsize > 0) r1 <- Matrix::head(r1[order(-abs(r1))], nsize)
    barplot(sort(r1, decreasing = TRUE),
      beside = FALSE, las = 3,
      ylab = "marginal correlation"
    )
    title(main, cex.main = 1.2)
  }

  if ("pcor" %in% what) {
    r2 <- M[, which(colnames(M) != "cor")] ## / ncol(M)
    if (nsize > 0) r2 <- Matrix::head(r2[order(-rowMeans(r2**2)), ], nsize)
    r2 <- r2[order(-rowMeans(r2)), ]
    r2 <- r2 / ncol(r2)
    par(mar = c(10, 4, 4, 2))
    barplot(t(r2),
      beside = FALSE, las = 3,
      ylim = c(-1, 1) * 0.3,
      ylab = "partial correlation (average)"
    )
    title(main, cex.main = 1.2)
    if (length(gene) == 1) {
      legend("topright", rev(colnames(r2)),
        fill = rev(grey.colors(ncol(r2))),
        cex = 0.85, y.intersp = 0.85
      )
    }
  }

  if ("graph" %in% what) {
    ## ------------------------------------------------------------
    ## extract core graph
    ## ------------------------------------------------------------
    R <- res$rho[["cor"]]
    if (nsize > 0) {
      top20 <- Matrix::head(unique(c(gene, rownames(M))), nsize)
      R <- res$rho[["cor"]][top20, top20] ## marginal correlation
    }
    gr1 <- igraph::graph_from_adjacency_matrix(
      abs(R),
      mode = "undirected", diag = FALSE, weighted = TRUE
    )
    ee <- igraph::get.edges(gr1, igraph::E(gr1))
    igraph::E(gr1)$rho <- R[ee]
    igraph::V(gr1)$name

    has.pcor <- !is.null(res$meta.pcor)
    P <- NULL
    if (has.pcor) {
      ii <- rownames(R)
      P <- res$meta.pcor[ii, ii]
      P[is.na(P)] <- 0
      P[is.nan(P)] <- 0
      igraph::E(gr1)$pcor <- P[ee]
    }
    if (is.null(rho.min)) rho.min <- 0

    ## threshold edges
    sel.ee <- (abs(igraph::E(gr1)$rho) >= rho.min)
    if (has.pcor) sel.ee <- sel.ee & (abs(igraph::E(gr1)$pcor) >= pcor.min)
    gr2 <- igraph::subgraph.edges(
      gr1, which(sel.ee),
      delete.vertices = FALSE
    )
    subgraphs <- igraph::decompose.graph(gr2)
    k <- which(sapply(sapply(subgraphs, function(g) igraph::V(g)$name), function(vv) gene %in% vv))
    gr2 <- subgraphs[[k]]

    ## set edge width
    igraph::V(gr2)$name
    igraph::E(gr2)$width <- edge.width * 1.0
    if (has.pcor) {
      pw <- abs(igraph::E(gr2)$pcor) / mean(abs(igraph::E(gr2)$pcor), na.rm = TRUE)
      igraph::E(gr2)$width <- edge.width * pw
    }

    ## calculate layout
    ly <- switch(layout,
      "fr" = igraph::layout_with_fr(gr2),
      "kk" = igraph::layout_with_kk(gr2, weights = 1 / igraph::E(gr2)$weight),
      "graphopt" = igraph::layout_with_graphopt(gr2),
      "tree" = igraph::layout_as_tree(gr2),
      igraph::layout_nicely(gr2)
    )
    rownames(ly) <- igraph::V(gr2)$name
    ly <- ly[igraph::V(gr2)$name, , drop = FALSE]

    add.alpha <- function(col, alpha) {
      apply(cbind(t(col2rgb(klr)), alpha), 1, function(x) {
        rgb(x[1], x[2], x[3], alpha = x[4], maxColorValue = 255)
      })
    }

    p1 <- 0.8 * igraph::E(gr2)$rho
    klrpal <- colorRampPalette(c("red3", "grey90", "grey30"))(64)

    klr <- klrpal[32 + 31 * p1]
    igraph::E(gr2)$color <- add.alpha(klr, 64 * abs(p1))

    j <- which(igraph::V(gr2)$name == gene)
    igraph::V(gr2)$label.cex <- 0.8
    igraph::V(gr2)$size <- 10
    igraph::V(gr2)$size[j] <- 16
    igraph::V(gr2)$label.cex[j] <- 1.1
    igraph::V(gr2)$color <- "skyblue"
    igraph::V(gr2)$frame.color <- "grey50"
    igraph::V(gr2)$frame.color[j] <- "black"

    par(mar = c(0, 0, 0, 0))
    plot(gr2, layout = ly)
  }
}


#' Test relationships between gene expression and traits
#'
#' @param me A gene expression matrix 
#' @param df A data frame containing trait data
#' @param plot Logical indicating whether to plot results. Default is TRUE.
#' @param cex Size of text in plots. Default is 1.
#'
#' @return A list with correlation and ANOVA results
#' 
#' @description
#' Tests for relationships between gene expression (me) and trait data (df),
#' using correlation for continuous variables and ANOVA for discrete variables.
#'
#' @details
#' This function takes a gene expression matrix (me) and a data frame (df) containing 
#' trait data. It separates df into continuous (numeric) and discrete (factor) variables.
#' 
#' For continuous variables, it calculates the correlation between each trait and each gene.
#' For discrete variables, it performs an ANOVA between each trait and each gene.
#' 
#' The results are returned as matrices of p-values for the correlations and ANOVAs.
#' If plot = TRUE, scatterplots and boxplots are generated to visualize the results.
#' 
#' @export
pgx.testTraitRelationship <- function(me, df, plot = TRUE, cex = 1) {
  df <- type.convert(df, as.is = TRUE)
  cl <- sapply(df, class)
  cvar <- which(cl %in% c("numeric", "integer"))
  dvar <- which(cl %in% c("factor", "character"))
  dc <- df[, cvar, drop = FALSE]
  dd <- df[, dvar, drop = FALSE]

  ## contious vs continous -> correlation
  rho.P <- NULL
  if (ncol(dc)) {
    rho.P <- matrix(NA, ncol(me), ncol(dc))
    i <- 1
    j <- 2
    rho <- stats::cor(me, dc, use = "pairwise")
    rho.P <- cor.pvalue(P, nrow(me))
  }

  ## continous vs discrete -> ANOVA
  anova.P <- NULL
  if (ncol(dd)) {
    anova.P <- matrix(NA, ncol(me), ncol(dd))
    colnames(anova.P) <- colnames(dd)
    rownames(anova.P) <- colnames(me)
    for (i in 1:ncol(dd)) {
      y <- dd[, i]
      res <- gx.limmaF(t(me), y, fdr = 1, lfc = 0)
      anova.P[, i] <- res[colnames(me), "P.Value"]
    }
    anova.P
  }

  P <- cbind(rho.P, anova.P)
  P <- P[, colnames(df)]

  if (plot == TRUE) {
    sigP <- -log10(P + 1e-8)
    sigP <- (1 - P)**1

    corrplot::corrplot(
      sigP,
      is.corr = FALSE, #
      mar = c(0, 0, 0, 2),
      tl.cex = cex, tl.col = "black", tl.offset = 1,
      cl.align.text = "l", cl.offset = 0.25, cl.cex = 0.7,
      pch.col = "grey50"
    )
  }
  return(P)
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
pgx.testPhenoCorrelation <- function(df, plot = TRUE, cex = 1) {
  cl <- sapply(df, class)
  nlev <- apply(df, 2, function(x) length(unique(x[!is.na(x)])))
  cvar <- which(cl %in% c("numeric", "integer") & nlev >= 2)
  dvar <- which(cl %in% c("factor", "character") & nlev >= 2)
  dc <- df[, cvar, drop = FALSE]
  dd <- df[, dvar, drop = FALSE]

  ## discrete vs discreate -> Fisher test
  fisher.P <- NULL
  if (ncol(dd)) {
    fisher.P <- matrix(NA, ncol(dd), ncol(dd))
    i <- 1
    j <- 2
    for (i in 1:(ncol(dd) - 1)) {
      kk <- which(!is.na(dd[, i]) & !is.na(dd[, j]))
      if (length(unique(dd[kk, i])) < 2 || length(unique(dd[kk, j])) < 2) next
      for (j in (i + 1):ncol(dd)) {
        tb <- table(dd[, i], dd[, j])
        fisher.P[i, j] <- fisher.test(tb, simulate.p.value = TRUE)$p.value
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
        kruskal.P[i, j] <- kruskal.test(dc[kk, j], dd[kk, i])$p.value
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
        cor.P[i, j] <- cor.test(dc[, i], dc[, j])$p.value
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
  qv <- p.adjust(P[ij], method = "BH")
  Q <- P
  Q[ij] <- qv

  P[is.na(P)] <- 0
  P <- (P + t(P)) / 2
  Q[is.na(Q)] <- 0
  Q <- (Q + t(Q)) / 2

  BLUERED <- colorRampPalette(c("blue3", "white", "red3"))

  if (plot == TRUE) {
    logP <- -log10(P + 1e-8)
    logQ <- -log10(Q + 1e-8)
    diag(logQ) <- 0
    corrplot::corrplot(
      corr = logQ,
      type = "upper",
      col = BLUERED(25),
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

  return(list(P = P, Q = Q))
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
