#' @describeIn knnImputeMissing Randomly impute missing values
#' @export
randomImputeMissing <- function(x) {
  i <- 1
  for (i in 1:ncol(x)) {
    jj <- which(is.na(x[, i]) | x[, i] == "NA")
    if (length(jj)) {
      rr <- sample(x[-jj, i], length(jj), replace = TRUE)
      x[jj, i] <- rr
    }
  }
  return(x)
}
#' @export
pgx.scQualityControlPlots.NOTFINISHED <- function(counts) {
  mt.genes <- grep("^MT-", rownames(counts), ignore.case = TRUE, value = TRUE)
  rb.genes <- grep("^RP[SL]", rownames(counts), ignore.case = TRUE, value = TRUE)
  percent.mito <- Matrix::colSums(counts[mt.genes, ]) / Matrix::colSums(counts) * 100
  percent.ribo <- Matrix::colSums(counts[rb.genes, ]) / Matrix::colSums(counts) * 100
  nfeature <- Matrix::colSums(counts > 0)
  ncounts <- Matrix::colSums(counts)

  T1 <- table(sample.id, qc.ok)
  barplot(t(T1), ylab = "number of cells")
  legend("topright",
    legend = c("qc=TRUE", "qc=FALSE"),
    fill = c("grey80", "grey30"), cex = 0.9
  )


  P <- data.frame(percent.mito, percent.ribo)
  P1 <- apply(P, 2, function(x) tapply(x, sample.id, mean))
  barplot(t(P1), ylab = "percentage (%)")
  legend("topright", legend = c("mito", "ribo"), fill = c("grey80", "grey30"), cex = 0.9)

  P2 <- apply(P, 2, function(x) tapply(x, cell.type, mean))
  barplot(t(P2), ylab = "percentage (%)")
  legend("topright", legend = c("mito", "ribo"), fill = c("grey80", "grey30"), cex = 0.9)
}
#' Activation matrix heatmap for PGX results
#'
#' @param pgx PGX object containing results
#' @param features Features to include (default NULL for all)
#' @param contrasts Contrasts to include (default NULL for all)
#' @param n Number of top features to display
#' @param qsize Color by adjusted p-value? (default TRUE)
#' @param cex Overall text size
#' @param cex.row Row text size
#' @param cex.col Column text size
#' @param srt Rotation angle for text labels
#' @param flip Logical to flip sign of fold changes
#' @param level Data level to use ("gene", "geneset")
#' @param clust.x Cluster columns? (default TRUE)
#' @param clust.y Cluster rows? (default TRUE)
#' @param plotlib Plotting library to use (default "base")
#'
#' @return A heatmap grob object
#'
#' @description
#' Generates an activation matrix heatmap from PGX results
#'
#' @details
#' This function takes a PGX object and generates a clustered heatmap
#' visualization of the top results. Fold changes are shown with color scale.
#' The top n features can be displayed, optionally clustering rows and columns.
#' Sign of changes can be flipped and results can be extracted at the gene or
#' gene set level. Colors can be based on adjusted p-values instead.
#'
#' Useful for visualizing patterns of activation across comparisons.
#'
#' @export
pgx.ActivationMatrix <- function(pgx, features = NULL, contrasts = NULL,
                                 n = 50, qsize = TRUE,
                                 cex = 1, cex.row = 1, cex.col = 1, srt = 90,
                                 flip = FALSE, level = "geneset",
                                 clust.x = TRUE, clust.y = TRUE, plotlib = "base") {
  if (level == "geneset") {
    out <- pgx.getMetaFoldChangeMatrix(pgx, what = "meta", level = "geneset")
  } else {
    out <- pgx.getMetaFoldChangeMatrix(pgx, what = "meta", level = "gene")
  }
  F <- out$fc
  Q <- out$qv
  if (!is.null(contrasts)) {
    ct <- intersect(contrasts, colnames(F))
    F <- F[, ct, drop = FALSE]
    Q <- Q[, ct, drop = FALSE]
  }
  if (!is.null(features)) {
    features <- intersect(features, rownames(F))
    F <- F[features, , drop = FALSE]
    Q <- Q[features, , drop = FALSE]
  }

  F1 <- Matrix::head(F[order(-apply(F, 1, sd)), ], n = n)
  F1 <- Matrix::head(F[order(-rowMeans(F**2)), ], n = n)

  ## cluster
  ii <- 1:nrow(F1)
  jj <- 1:ncol(F1)
  if (clust.y) ii <- fastcluster::hclust(dist(F1))$order
  if (clust.x) jj <- fastcluster::hclust(dist(t(F1)))$order
  F1 <- F1[ii, jj]
  Q1 <- Q[rownames(F1), colnames(F1)]

  cpal <- rev(RColorBrewer::brewer.pal(11, "RdYlBu"))
  cpal <- colorRampPalette(c("blue4", "lightskyblue1", "lightyellow", "rosybrown1", "red4"))(64)
  cpal <- colorspace::diverge_hcl(64, c = 60, l = c(30, 100), power = 1)

  p <- NULL
  if (plotlib == "base") {
    F2 <- F1
    if (qsize) F2 <- F2 * Q1
    if (flip) {
      F2 <- t(F2)
    }
    colnames.f2 <- colnames(F2)
    colnames(F2) <- rep("", ncol(F2))
    corrplot::corrplot(
      F2, #
      is.corr = FALSE, col = cpal,
      tl.col = "grey0", tl.cex = 0.8 * cex.row,
      tl.pos = "lt", ## tl.srt=45, tl.offset=1,
      cl.cex = 0.6,
      mar = c(10, 0, 1, 0)
    )
    text(1:ncol(F2), 0, colnames.f2,
      srt = srt, xpd = TRUE, adj = 1, cex = 0.8 * cex.col
    )
  }

  df <- NULL
  if (plotlib %in% c("ggplot", "plotly")) {
    df <- reshape2::melt(F1)
    df$size <- df$value
    tt.size <- "value"
    if (qsize) {
      df$size <- -log10(0.0001 + Q[cbind(df$Var1, df$Var2)])
      tt.size <- "-log10q"
    }
    p <- ggplot2::ggplot(ggplot2::aes(x = Var2, y = Var1, color = value, size = size), data = df) +
      ggplot2::theme_minimal() +
      ggplot2::geom_point() +
      ggplot2::scale_size(range = c(0.1, 5 * cex)) +
      ggplot2::scale_color_gradientn(colors = cpal) +
      ggplot2::xlab(NULL) +
      ggplot2::ylab(NULL) +
      ggplot2::labs(size = tt.size, color = "value") +
      ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = srt)) +
      ggplot2::theme(
        plot.margin = ggplot2::margin(5, 0, 0, 10),
        legend.text = ggplot2::element_text(size = 9),
        legend.key.size = grid::unit(10, "pt"),
        legend.key.height = grid::unit(12, "pt")
      )


    if (!qsize) {
      p <- p + ggplot2::guides(size = FALSE)
    }
    if (flip) {
      p <- p + ggplot2::coord_flip() +
        ggplot2::scale_x_discrete(guide = ggplot2::guide_axis(angle = 0)) +
        ggplot2::scale_y_discrete(guide = ggplot2::guide_axis(angle = srt))
    }
  }
  if (plotlib == "plotly") {
    p <- p + ggplot2::scale_size(range = c(0.1, 3 * cex))
    p <- plotly::ggplotly(p) %>%
      plotly::layout(xaxis = list(tickangle = -srt, side = "bottom"))
  }
  p
}
#' @export
pgx.plotForwardProjection <- function(gr, gene, cex = 1, fx = NULL,
                                      features = NULL, main = NULL, plot = TRUE) {
  if (!is.null(features)) {
    gr <- pgx.createVipGeneLayer(gr, genes)
  }

  vtype <- gsub("\\}.*|^\\{", "", rownames(gr$layout))
  tsne_genes <- gr$layout[which(vtype == "gene"), ]
  tsne_gsets <- gr$layout[which(vtype == "geneset"), ]

  uscale <- function(x) (x - min(x)) / (max(x) - min(x)) - 0.5
  pos1 <- apply(tsne_gsets[, 1:2], 2, uscale)
  pos2 <- apply(tsne_genes[, 1:2], 2, uscale)
  pos1 <- t(t(pos1) + c(+0.6, 0))
  pos2 <- t(t(pos2) + c(-0.6, 0))

  ## ---------------- get all edges/paths ---------------


  tt <- ""
  to <- from <- NULL
  if (!is.null(geneset)) {
    gs <- paste0("{geneset}", geneset)

    if (gs %in% igraph::V(gr)$name) {
      nb <- igraph::V(gr)[igraph::neighbors(gr, gs)]$name
      gg <- intersect(nb, rownames(pos2))
      from <- pos1[gs, ]
      to <- pos2[gg, , drop = FALSE]
    }
    tt <- geneset
  }

  if (!is.null(gene)) {
    gg <- paste0("{gene}", gene)

    if (gg %in% igraph::V(gr)$name) {
      nb <- igraph::V(gr)[igraph::neighbors(gr, gg)]$name
      gs <- intersect(nb, rownames(pos1))
      from <- pos2[gg, ]
      to <- pos1[gs, , drop = FALSE]
    } else {
      cat("WARNING::", gg, " not in omicsgraph\n")
    }
    tt <- gene
  }

  if (plot == TRUE) {
    cex1 <- 1 + 1 * (rownames(pos1) %in% igraph::V(gr)$name)
    klr1 <- c("grey70", "grey10")[cex1]
    cex2 <- 1 + 1 * (rownames(pos2) %in% igraph::V(gr)$name)
    klr2 <- c("grey70", "grey10")[cex2]
    if (!is.null(fx)) {
      fx1 <- fx[match(rownames(pos1), names(fx))]
      fx2 <- fx[match(rownames(pos2), names(fx))]


      fx1 <- fx1 / max(abs(fx1), na.rm = TRUE)
      fx2 <- fx2 / max(abs(fx2), na.rm = TRUE)
      klr1 <- gplots::bluered(32)[16 + round(15 * fx1)]
      klr2 <- gplots::bluered(32)[16 + round(15 * fx2)]
      ix1 <- cut(fx1, breaks = c(-99, -0.1, 0.1, 99))
      ix2 <- cut(fx2, breaks = c(-99, -0.1, 0.1, 99))
      klr1 <- c("blue", "gray40", "red")[as.integer(ix1)]
      klr2 <- c("blue", "gray40", "red")[ix2]
      klr1[which(is.na(klr1))] <- "grey80"
      klr2[which(is.na(klr2))] <- "grey80"
    }

    cex <- 0.02
    if (nrow(pos1) < 1000) cex <- 0.4
    pch <- "."
    pch <- 20
    j1 <- 1:length(klr1)
    if (length(klr1) > 5000) {
      j1 <- c(sample(grep("grey", klr1), 5000), which(klr1 != "grey"))
    }
    plot(pos1[j1, ],
      pch = pch, xlim = c(-1.1, 1.1), ylim = c(-0.5, 0.5),
      xaxt = "n", yaxt = "n", xlab = "", ylab = "", bty = "n",
      col = klr1[j1], cex = cex * cex1[j1]
    )

    j2 <- 1:nrow(pos2)
    if (length(klr2) > 5000) {
      j2 <- c(sample(grep("grey", klr2), 5000), which(klr2 != "grey"))
    }
    points(pos2[j2, ], pch = pch, cex = cex * cex2[j2], col = klr2[j2], main = "genes")

    legend("bottomleft", "GENES", bty = "n", col = "grey60")
    legend("bottomright", "GENE SETS", bty = "n", col = "grey60")
    if (!is.null(from)) {
      points(from[1], from[2], pch = 20, cex = 1, col = "green3")
      points(to[, 1], to[, 2], pch = 20, cex = 0.4, col = "green3")
      arrows(from[1], from[2], to[, 1], to[, 2],
        length = 0.05,
        lwd = 0.5, col = paste0(gplots::col2hex("green3"), "33")
      )
    }
    if (!is.null(main)) tt <- main
    mtext(tt, line = 0.5, at = 0, font = 2, cex = 1.1)
  }
  invisible(rownames(to))
}
#' @title Infer cell types using LM22 signature
#'
#' @description This function infers the proportion of cell types in a bulk
#' RNA-seq sample using the LM22 signature.
#'
#' @param counts A counts matrix with genes as rows and samples as columns.
#' @param min.count The minimum count threshold. Default is 3.
#' @param markers A data frame with positive and negative cell type markers.
#'
#' @details This function checks if the canonical cell type markers from the
#' LM22 signature are expressed above a minimum threshold in the input counts
#' matrix. It returns a list with logical vectors indicating which cell types
#' are detected based on their marker expression.
#'
#' @return A list with logical vectors for each cell type indicating marker presence.
#'
#' @export
pgx.inferCellTypeLM22 <- function(counts, low.th = 0.01, add.unknown = FALSE,
                                  normalize.mat = TRUE, method = "NNLM",
                                  min.count = 3, celltype0 = NULL, min.prob = 0.2) {
  ## Two-pass (2-level) cell type identification using LM22
  M <- playdata::LM22
  M <- as.matrix(M)
  colnames(M) <- gsub(" ", "_", sub(" ", ".", sub(" ", "_", colnames(M))))
  colnames(M) <- sub("Macrophages_", "Macrophages.", colnames(M))

  P0 <- NULL
  if (is.null(celltype0)) {
    res <- pgx.inferCellType(counts,
      low.th = 0, add.unknown = add.unknown,
      normalize.mat = TRUE, method = "NNLM",
      min.prob = min.prob
    )
    celltype0 <- res$celltype
    P0 <- res$probs
  }

  ## Filter count matrix
  X <- counts
  X <- X[(rowMeans(X >= min.count) > low.th), ] ## OK???

  ## Match matrices
  rownames(X) <- toupper(rownames(X))
  rownames(M) <- toupper(rownames(M))
  gg <- intersect(rownames(X), rownames(M))
  X <- X[gg, ]
  M <- M[gg, ]

  ## 2nd stage

  celltype <- rep(NA, ncol(X))
  ct <- "Macrophages"
  ct <- "Monocytes"
  ct <- "T_cells"
  P1 <- matrix(NA, nrow = ncol(X), ncol = 0)
  rownames(P1) <- colnames(X)
  for (ct in unique(celltype0)) {
    msel <- grep(ct, colnames(M))
    M1 <- M[, msel, drop = FALSE]
    jj <- which(celltype0 == ct)
    jj <- 1:length(celltype0)
    ct3 <- celltype0[jj]

    if (ncol(M1) > 0) {
      X1 <- X[, jj, drop = FALSE]
      X1 <- X1 / (1e-3 + rowMeans(X1)) ## center feature means??
      M1 <- M1 / (1e-3 + rowMeans(M1)) ## center feature means??
      res1 <- pgx.deconvolution(
        X1,
        ref = M1, methods = "NNLM", normalize.mat = TRUE,
        add.unknown = FALSE
      )
      C1 <- res1$results[["NNLM"]]
      C1 <- C1 / rowSums(C1)
      ct3 <- colnames(C1)[max.col(C1)]
      C0 <- P0[colnames(X1), ct]
      C1 <- C1 * C0
      C1 <- C1[match(rownames(P1), rownames(C1)), , drop = FALSE]
      C1[is.na(C1)] <- 0
      P1 <- cbind(P1, C1)
    } else {
      C1 <- P0[, ct]
      P1 <- cbind(P1, C1)
      colnames(P1)[ncol(P1)] <- ct
    }
    celltype[jj] <- ct3
  }

  celltype <- colnames(P1)[max.col(P1)]

  ## collapse small groups to 'other_cells'
  low.ct <- names(which(table(celltype) < low.th * length(celltype)))
  celltype[celltype %in% low.ct] <- "other_cells"
  names(celltype) <- colnames(counts)

  res <- list(celltype = celltype, probs = P1)
  return(res)
}
#' @export
fastcorSAVE <- function(x) {
  sx <- scale(x)
  (t(sx) %*% sx) / (nrow(x) - 1)
}
#' @title Purify tumor expression profiles
#'
#' @param expr Numeric gene expression matrix, rows are genes, columns are samples.
#' @param normal Numeric gene expression matrix from normal samples.
#' @param method Character vector of deconvolution methods to apply.
#'
#' @return List containing purified expression profiles and estimated tumor purity.
#'
#' @description Estimates and removes normal cell contamination from tumor expression profiles.
#'
#' @details This function takes a tumor expression matrix and an expression matrix from normal samples.
#' It applies computational deconvolution methods like ISOpure, DeMixT, etc. to estimate normal contamination.
#'
#' The normal profile is removed from each tumor sample to generate a purified tumor expression profile.
#' Results from all methods are returned as a list, along with numeric vectors of estimated tumor purity.
#'
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
    normal.frac <- (normalX %*% cf)
    alpha0 <- (1 - Matrix::colSums(normal.frac) / Matrix::colSums(tumorX))
    alpha0

    xhat[["nnlm"]] <- pmax(tumorX - normal.frac, 0)
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

    x.hat <- ISOpureS2model$cc_cancerprofiles
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
#' Remove principal components from gene expression data
#'
#' @title Remove principal components from gene expression data
#'
#' @description
#' Removes the top principal components from a gene expression matrix.
#'
#' @param X Gene expression matrix, genes in rows, samples in columns.
#' @param nv Number of principal components to remove.
#'
#' @details
#' This function calculates the top \code{nv} principal components of the gene expression matrix \code{X} using \code{irlba::irlba}.
#' It then removes these principal components from \code{X} using \code{limma::removeBatchEffect}.
#'
#' The goal is to remove technical noise or batch effects captured by the top principal components.
#'
#' @return The gene expression matrix \code{X} with the top \code{nv} principal components removed.
#'
#' @export
pgx.removePC <- function(X, nv) {
  suppressWarnings(suppressMessages(
    pc <- irlba::irlba(X, nv = nv)$v
  ))
  cX <- limma::removeBatchEffect(X, covariates = pc)
  return(cX)
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
#' Check Cell Type Markers
#'
#' @description
#' Checks if known cell type marker genes are detected (positive markers)
#' or not detected (negative markers) in a gene expression matrix.
#'
#' @param counts Gene expression matrix with genes in rows and samples in
#' columns.
#' @param min.count Minimum expression count threshold. Default 3.
#' @param markers Named list of cell type marker genes.
#'
#' @details
#' This function takes a gene expression matrix and checks if known cell type
#' marker genes are detected above a minimum expression threshold (positive markers)
#' or not detected (negative markers).
#'
#' It compares the matrix rownames against lists of known canonical markers for major
#' immune and cancer cell types. Markers can be specified, otherwise built-in
#' canonical markers are used.
#'
#' For each cell type, it identifies the positive and negative markers based on
#' expression above or below min.count. This provides an overview of which cell types
#' are potentially present or absent based on their expression markers in the data.
#'
#' @return Named list with logical vectors indicating detected (TRUE) vs undetected
#' (FALSE) markers for each cell type.
#'
#' @export
pgx.checkCellTypeMarkers <- function(counts, min.count = 3, markers = NULL) {
  ## Checks if markers are expressed (positive markers) or
  ## not-expressed (negative markers).
  ##
  ##
  CANONICAL.MARKERS.SAVE <- list(
    "B cells" = c("Ms4a1", "CD79a", "CD79b", "Fcmr", "Ebr1"),
    "CD4 T cells" = c("Il17r", "Ccr7", "Cd3e", "Cd3d", "Cd3g"),
    "CD8 T cells" = c("Cd8a", "Cd8b"),
    "NK cells" = c("Gnly", "Nkg7", "Gzma", "Klrb1c", "Klrk1", "Klra4"),
    "Dendritic cells" = c("Fcer1a", "Cst3", "Siglech", "Fscn1", "Ccl22"),
    "Macrophages" = c("C1qa", "C1qb", "C1qc", "Lyz2"),
    "Monocytes" = c("Ly6c2", "Lyz", "Cd14", "Fcgr3a", "Ms4a7"),
    "B16.melanoma" = c("Mlana", "Dct", "Tyrp1", "Mt1", "Mt2", "Pmel", "Pgk1")
  )
  CANONICAL.MARKERS2 <- list(
    "B_cells" = c("Ms4a1+", "Cd79a+", "Cd79b+"),
    "T_cells" = c("Cd3e+", "Cd3d+", "Cd3g+"),
    "NK_cells" = c("Gzma+", "Klrb1c+", "Klra4+"),
    "Dendritic_cells" = c("Cst3+", "H2-Ab1+", "Adgre1-"),
    "Macrophages" = c("C1qa+", "C1qb+", "C1qc+"),
    "Monocytes" = c("Ly6c2+", "Lyz2+", "Ccr2+"),
    "Immune_cell" = c("Ptprc+")
  )
  if (is.null(markers)) {
    markers <- CANONICAL.MARKERS2
  }

  marker.genes <- sort(unique(unlist(markers)))
  marker.genes <- gsub("[-+]$", "", marker.genes)

  M <- matrix(0, nrow = length(marker.genes), ncol = length(markers))
  dimnames(M) <- list(marker.genes, names(markers))
  k <- 7
  k <- 1
  for (k in 1:ncol(M)) {
    mm <- markers[[k]]

    m.neg <- sub("[-]$", "", grep("[-]$", mm, value = TRUE))
    m.pos <- setdiff(gsub("[+-]$", "", mm), m.neg)
    if (length(m.pos)) M[match(m.pos, rownames(M)), k] <- +1
    if (length(m.neg)) M[match(m.neg, rownames(M)), k] <- -1
  }

  X <- counts
  rownames(M) <- toupper(rownames(M))
  rownames(X) <- toupper(rownames(X))
  gg <- intersect(rownames(M), rownames(X))
  counts0 <- counts[match(gg, toupper(rownames(counts))), , drop = FALSE]

  X1 <- 1 * (X[gg, , drop = FALSE] >= min.count)
  M1 <- M[gg, , drop = FALSE]
  check.pos <- t(pmax(M1, 0)) %*% X1
  check.neg <- t(pmin(M1, 0)) %*% X1
  check.pos <- t(M1 == +1) %*% X1 / Matrix::colSums(M1 == 1)
  check.neg <- t(M1 == -1) %*% X1 / Matrix::colSums(M1 == -1)
  check <- 1 * t(check.pos & check.neg)

  return(list(check = check, marker.expr = counts0, markers = markers))
}
#' @export
pgx.computeShortestPath <- function(graph, contrast, niter = 1, r = 0.01,
                                    output = "vpath") {
  ## add source/sink
  graph <- pgx._addSourceSink(graph)

  ## calculate weights
  fc <- graph$foldchange[, contrast]
  ee <- igraph::get.edges(graph, igraph::E(graph))
  f1 <- pmax(fc[ee[, 1]], 0)
  f2 <- pmax(fc[ee[, 2]], 0)
  node.values <- sqrt(f1 * f2) ## strictly positive
  edge.rho <- pmax(igraph::E(graph)$weight, 0)
  score <- node.values * edge.rho
  weights0 <- -log(pmax(score / max(score), 1e-8))

  ## ----------------------------------------------------------
  ## solve single SP (SOURCE to SINK)
  ## ----------------------------------------------------------
  vpath <- list()
  epath <- list()
  if (niter > 0) {
    for (i in 1:niter) {
      sd0 <- r * sd(weights0)
      weights1 <- 0.0 + pmax(weights0 + sd0 * rnorm(length(weights0)), 0)

      system.time(
        sp.out <- igraph::shortest_paths(graph,
          from = "SOURCE", to = "SINK", mode = "all",
          weights = weights1, output = output
        )
      )
      nv <- length(sp.out$vpath[[1]])
      ne <- length(sp.out$epath[[1]])
      vpath[[i]] <- sp.out$vpath[[1]]$name[2:(nv - 1)]
      epath[[i]] <- sp.out$epath[[1]][2:(ne - 1)]
    }
  }

  vfreq <- sort(table(unlist(vpath)), decreasing = TRUE)



  res <- list(vpath = vpath, epath = epath, vfreq = vfreq)
  return(res)
}
#' @export
silac.plotDistribution <- function(obj, samples, minq = 3, main = NULL) {
  if (!is.null(samples) && inherits(samples[1], "character") &&
    samples[1] %in% names(obj$groups)) {
    if (is.null(main)) main <- samples
    samples <- obj$groups[[samples]]
  }
  ## Ranking according to the average at 24h
  Q <- as.matrix(obj$LFQ.ratio[, samples])
  valid <- (rowSums(Q > 0) >= minq)
  Q <- Q[valid, ]
  meanQ <- rowMeans(Q, na.rm = TRUE)
  plot(sort(meanQ, decreasing = TRUE), main = main)
  hist(meanQ, breaks = 20, col = "darkgrey", main = main)
}
#' @export
fastcor <- function(a) {
  a <- Matrix::Matrix(a)
  n <- nrow(a)
  muX <- Matrix::colMeans(a)
  system.time(covmat <- (as.matrix(Matrix::crossprod(a)) - n * Matrix::tcrossprod(muX)) / (n - 1))
  sdvec <- sqrt(diag(covmat))
  r <- covmat / tcrossprod(sdvec)
  r
}
#' @export
pgx.makeTriSystemGraph <- function(data, Y, nfeat = 25, numedge = 100, posonly = FALSE) {
  if (is.null(names(data))) stop("X.list must be named")
  if (!all(sapply(data, ncol) == length(Y))) {
    stop("data columns must match Y")
  }
  if (is.null(posonly)) posonly <- FALSE

  .getLoadingImportance <- function(k, res.pls, nfeat) {
    U <- res.pls$loadings[[k]]
    V <- res.pls$variates[[k]]
    U <- U[which(rowSums(abs(U)) > 0), ]
    U <- U[order(-rowSums(U * U)), ]
    y <- res.pls$Y
    vm <- apply(V, 2, function(x) tapply(x, y, mean))
    R0 <- U %*% t(V)
    R <- U %*% t(vm)

    R <- Matrix::head(R[order(-rowSums(R**2)), ], nfeat)

    W <- exp(1.5 * R / sd(R))

    W <- W / (1e-5 + rowSums(W))
    W <- W * rowSums(R * R, na.rm = TRUE)**0.5
    return(W)
  }

  .detectSelfLoops <- function(gr, posonly) {
    v1 <- grep("1:", igraph::V(gr)$name, value = TRUE)
    v3 <- grep("3:", igraph::V(gr)$name, value = TRUE)

    wt <- (igraph::E(gr)$rho * igraph::E(gr)$importance)
    wt <- wt / max(abs(wt))
    if (posonly) wt <- pmax(wt, 0)
    wt <- -log(pmin(1e-8 + abs(wt), 1))

    self.epath <- vector("list", length(v1))
    self.vpath <- vector("list", length(v1))
    self.dist <- rep(NA, length(v1))
    names(self.vpath) <- names(self.epath) <- names(self.dist) <- v1
    i <- 1
    for (i in 1:length(v1)) {
      d1 <- igraph::distances(gr, v1[i], to = v3, mode = "out", weights = wt)
      d3 <- igraph::distances(gr, v1[i], to = v3, mode = "in", weights = wt)
      suppressWarnings(sp1 <- igraph::shortest_paths(
        gr, v1[i],
        to = v3, weights = wt,
        mode = "out", output = "both"
      ))
      suppressWarnings(sp3 <- igraph::shortest_paths(
        gr, v1[i],
        to = v3, weights = wt,
        mode = "in", output = "both"
      ))
      dd <- d1 + d3
      j <- which.min(dd)
      self.dist[i] <- dd[j]
      self.vpath[[i]] <- NA
      self.epath[[i]] <- NA
      if (!is.infinite(dd[j])) {
        self.vpath[[i]] <- c(sp1$vpath[[j]]$name, v1[i])
        self.epath[[i]] <- c(sp1$epath[[j]], sp3$epath[[j]])
      }
    }
    jj <- order(self.dist)
    res <- list(
      dist = self.dist[jj], vpath = self.vpath[jj],
      epath = self.epath[jj]
    )
    return(res)
  }


  .makeTriGraph <- function(res.pls, data, nfeat, numedge, posonly) {
    cat("<makeTriSystemGraph:.makeTriGraph> called\n")
    W.list <- list()
    X.list <- list()
    k <- 1
    for (k in 1:3) {
      W.list[[k]] <- .getLoadingImportance(k, res = res.pls, nfeat = nfeat)
      gg <- rownames(W.list[[k]])
      X.list[[k]] <- data[[k]][gg, ]
      rownames(X.list[[k]]) <- rownames(W.list[[k]])
    }
    names(W.list) <- names(X.list) <- names(res.pls$X)

    cat("<makeTriSystemGraph:.makeTriGraph> 1\n")
    cat("<makeTriSystemGraph:.makeTriGraph> posonly=", posonly, "\n")


    edge.list <- c()
    k <- 1
    for (k in 1:3) {
      ## add correlation lines
      p1 <- rownames(W.list[[k]])
      m <- ifelse(k < 3, k + 1, 1)
      p2 <- rownames(W.list[[m]])
      rho <- stats::cor(res.pls$X[[k]][, p1], res.pls$X[[m]][, p2])
      if (posonly) rho <- pmax(rho, 0)
      q0 <- -1
      if (numedge > 0) q0 <- Matrix::tail(sort(abs(rho)), numedge)[1]
      idx <- which(abs(rho) > q0, arr.ind = TRUE)
      pp <- NULL
      pp <- data.frame(from = p1[idx[, 1]], to = p2[idx[, 2]], rho = rho[idx])

      jj <- which(sub("[0-9]:", "", pp$from) == sub("[0-9]:", "", pp$to))
      if (length(jj)) {
        pp$rho[jj] <- 0.01 ## no self loops across levels
      }
      self.edges <- data.frame(from = p1, to = p1, rho = 0)
      edge.list[[k]] <- rbind(pp, self.edges)
    }

    cat("<makeTriSystemGraph:makeGraph> making graph object\n")

    ee <- do.call(rbind, edge.list)
    gr <- igraph::graph_from_edgelist(as.matrix(ee[, 1:2]), directed = TRUE)
    ww <- lapply(W.list, function(w) rowMeans(w * w)**0.5)
    names(ww) <- NULL
    ww <- unlist(ww)
    igraph::V(gr)$importance <- ww[igraph::V(gr)$name] ## node importance

    ## set importance as edge weights
    igraph::E(gr)$rho <- ee$rho
    igraph::E(gr)$importance <- abs(ee$rho) * sqrt(ww[ee$from] * ww[ee$to])


    gr <- igraph::delete_edges(gr, igraph::E(gr)[which_loop(gr)])
    edges <- data.frame(igraph::get.edgelist(gr),
      rho = igraph::E(gr)$rho,
      importance = igraph::E(gr)$importance
    )
    colnames(edges)[1:2] <- c("from", "to")

    gr
    loops <- .detectSelfLoops(gr, posonly)
    top.loops <- Matrix::head(loops$epath, 10)
    top.loops
    edges$looping <- (igraph::E(gr) %in% unlist(top.loops))
    igraph::E(gr)$looping <- (igraph::E(gr) %in% unlist(top.loops))
    cat("<makeTriSystemGraph:makeGraph> done!\n")

    names(W.list) <- names(res.pls$X)
    out <- list(graph = gr, edges = edges, W = W.list)
    return(out)
  }

  ## ----------------------------------------------------------------------
  ## ----------------------------------------------------------------------
  ## ----------------------------------------------------------------------

  ## prepend index
  for (i in 1:length(data)) {
    rownames(data[[i]]) <- paste0(i, ":", rownames(data[[i]]))
  }

  ## set number of components and features
  NCOMP <- 3
  nfeat1 <- min(nfeat, nrow(data[[1]]))
  nfeat2 <- min(nfeat, nrow(data[[2]]))
  nfeat3 <- min(nfeat, nrow(data[[3]]))

  list.keepX <- list(rep(nfeat1, NCOMP), rep(nfeat2, NCOMP), rep(nfeat3, NCOMP))
  names(list.keepX) <- names(data)

  ## set up a full design where every block is connected
  design <- matrix(1,
    ncol = length(data), nrow = length(data),
    dimnames = list(names(data), names(data))
  )
  diag(design) <- 0
  design[1, 3] <- 0.5
  design[3, 1] <- 0.5
  design


  cat("<makeTriSystemGraph> calling block.splsda!\n")
  res.pls <- mixOmics::block.splsda(
    X = lapply(data, t), Y = Y, ncomp = NCOMP,
    keepX = list.keepX, design = design
  )

  cat("<makeTriSystemGraph> calling .makeTriGraph\n")

  res.gr <- .makeTriGraph(
    res.pls, data,
    nfeat = nfeat, numedge = numedge, posonly = posonly
  )

  res <- c(res.pls, res.gr)
  attr(res, "class") <- attr(res.pls, "class")
  return(res)
}
#' @title Prepare count data for DESeq2 differential analysis
#'
#' @description This function takes raw count data and prepares it for differential
#' expression analysis with DESeq2. It handles filtering, normalization, and batch correction.
#'
#' @param counts Data frame with raw counts, or named list with counts, samples, and genes components.
#' @param samples Data frame with sample metadata, required if counts is just a matrix.
#' @param genes Data frame with gene metadata, required if counts is just a matrix.
#' @param remove.batch Logical indicating whether to estimate and remove batch effects.
#' @param test Statistical test to use, either "Wald" or "LRT".
#' @param prior.cpm Filtering threshold for counts per million. Default is 0.
#' @param filter Logical indicating whether to filter low counts.
#'
#' @return A DESeqDataSet object ready for analysis with DESeq2.
#'
#' @details This function handles conversion of raw count data into a properly formatted
#' DESeqDataSet object for differential analysis with DESeq2. It checks that the sample and gene
#' metadata matches the counts. Filtering, normalization, and batch correction can be applied.
#'
#' The returned DESeqDataSet contains the normalized, filtered count matrix, along with sample and
#' gene data frames, ready for analysis with DESeq2 functions like DESeq() and results().
#'
#' @seealso
#' \code{\link[DESeq2]{DESeqDataSet}}, \code{\link[DESeq2]{DESeq}}, \code{\link[DESeq2]{results}}
#'
#' @export
ngs.cookForDESEQ2 <- function(counts, samples, genes, remove.batch = TRUE,
                              test = "Wald", prior.cpm = 0, filter = TRUE) {
  if (all(c("counts", "samples", "genes") %in% names(counts))) {
    samples <- counts$samples
    genes <- counts$genes
    counts <- counts$counts
  }

  #
  #
  if (!all(colnames(counts) == rownames(samples))) stop("samples do not match")
  if (!all(rownames(counts) == rownames(genes))) stop("genes do not match")
  if (is.null(samples)) stop("need samples specified")
  if (is.null(genes)) stop("need genes specified")

  ## ------------------------------------------------------------------
  ## remove genes with zero values (deseq does this like this?)
  ## ------------------------------------------------------------------
  if (filter) {
    keep <- rowSums(counts) > 0
    keep <- (rowSums(counts) > 100) > 3
    keep <- rowSums(counts) > 0 ## check original count because of prior.cpm
    table(keep)
    counts <- counts[keep, ]
    genes <- genes[keep, ]
  }

  ## ------------------------------------------------------------------
  ## prior count regularization
  ## ------------------------------------------------------------------

  prior.cpm
  if (!is.null(prior.cpm) && prior.cpm > 0) {
    cat("adding prior counts at prior.cpm=", prior.cpm, "\n")
    CPM.FACTORS <- Matrix::colSums(counts) / 1e6
    prior.counts <- (prior.cpm * CPM.FACTORS)
    prior.counts <- pmax(prior.counts, 1)
    summary(prior.counts)
    #
    counts <- t(t(counts) + prior.counts)
  }

  ## ------------------------------------------------------------
  ## Now create an DEseq2 object (see also tximport Vignette)
  ## ------------------------------------------------------------
  Matrix::head(samples)

  if (!("group" %in% colnames(samples))) {
    stop("samples must have 'group' column for DESeq")
  } else {
    cat("found 'group' column in sample table\n")
  }
  counts <- round(counts)
  mode(counts) <- "integer"
  colnames(samples)
  nbatch <- sum(grepl("batch", colnames(samples)))
  nbatch
  if (remove.batch && nbatch > 1) {
    ## It is good practice to add the batch parameters in the
    ## final model, even if they have been remove before. But with
    ## multiple batch effects it often creates '(in)dependency'
    ## problems...
    ##
    cat("found multiple 'batch' columns in sample table\n")
    batch <- colnames(samples)[grep("^batch", colnames(samples))]
    if (length(batch) > 1) {
      ## ONLY ONE BATCH FOR NOW!!!!!
      batch <- batch[1]
      cat("WARNING: only correcting one batch implemented\n")
    }
    design.formula <- formula(paste(c("~ 0 + group", batch), collapse = " + ")) ## multiple batch
    ## } else if(is.null(design)) {
  } else {
    design.formula <- formula(" ~ 0 + group")
  }
  cat("using model design: ", as.character(design.formula), "\n")

  rownames.counts <- rownames(counts)
  rownames(counts) <- NULL
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts, design = design.formula, colData = data.frame(samples)
  )
  rownames(counts) <- rownames.counts
  ## dds <- DESeq2::DESeqDataSetFromMatrix(
  ## counts, design = ~ 0 + group + batch, colData = samples[,c("group","batch")])

  ## ------------------------------------------------------------------
  ## to collapse TECHNICAL replicates???
  ## ------------------------------------------------------------------
  if (FALSE && "replicate" %in% colnames(samples)) {
    repl.id <- samples$replicate
    dds <- DESeq2::collapseReplicates(dds, repl.id, dds$sample)
  }

  ## The following code estimates size factors to account for
  ## differences in sequencing depth. So the value are typically
  ## centered around 1. If all the samples have exactly the same
  ## sequencing depth, you expect these numbers to be near 1.
  ## sf <- DESeq2::estimateSizeFactorsForMatrix(DESeq2::counts(dds))

  ## Run DESeq : Modeling counts with generic 'group'
  fitType <- "parametric"
  if (prior.cpm != 0) fitType <- "local"
  dds <- DESeq2::DESeq(dds, fitType = fitType, test = test)
  DESeq2::resultsNames(dds) # lists the coefficients

  ## we add the gene annotation here (not standard...)
  SummarizedExperiment::rowData(dds)$genes <- genes ## does this work??

  return(dds)
}
#' Differential expression testing between groups
#'
#' @title Differential expression testing between groups
#'
#' @description Performs differential expression testing between two groups.
#' Supports wilcoxon, t-test, limma, and fisher's exact test.
#'
#' @param sig Gene expression matrix with genes in rows and samples in columns.
#' @param class.label Vector of group labels for each sample. Must have two unique values.
#' @param fdr FDR threshold for identifying differentially expressed genes.
#' @param test.method Test method to use. Options are "wilcox", "limma", "ttest", or "fisher".
#' @param running.name Optional name prefix for output files.
#' @param output.to.file Logical for writing results to file.
#'
#' @details This function performs differential expression testing between two groups defined by \code{class.label}.
#' It handles filtering out NA values and outputs p-values and FDR for each gene.
#' The test method can be specified, with "limma" as the default.
#'
#' @return Data frame with gene p-values and FDR.
#'
#' @export
gx.test.groups <- function(sig, class.label, fdr = 0.20,
                           test.method = c("wilcox", "limma", "ttest", "fisher"),
                           running.name = NULL,
                           output.to.file = TRUE) {
  bb <- sort(unique(as.character(class.label)))
  if (length(bb) != 2) {
    stop("meti.limma:: currently only 2 classes")
  }
  if (length(test.method) > 1) test.method <- test.method[1]

  ## cleanup NA
  kk <- which(!is.na(class.label))
  sig <- sig[, kk]
  class.label <- class.label[kk]

  ## do LIMMA
  pv <- rep(NA, nrow(sig))
  names(pv) <- rownames(sig)
  if (test.method == "limma") {
    cat("performing Limma test\n")

    design <- cbind(1, class.label == bb[2])
    if (is.null(running.name)) running.name <- paste(bb[2], bb[1], sep = "vs")
    colnames(design) <- c("WT", running.name)
    d1 <- colnames(design)[1]
    d2 <- colnames(design)[2]
    fit <- limma::lmFit(sig, design)
    fit <- limma::eBayes(fit)
    tt <- limma::topTable(fit, coef = d2, number = nrow(sig))
    pv <- tt$P.Value[match(names(pv), tt$ID)]
  } else if (test.method == "wilcox") {
    cat("performing Wilcox rank test\n")
    for (i in 1:nrow(sig)) {
      pv[i] <- wilcox.test(sig[i, ] ~ class.label)$p.value
    }
  } else if (test.method == "ttest") {
    cat("performing T-test\n")
    for (i in 1:nrow(sig)) {
      pv[i] <- t.test(sig[i, ] ~ class.label)$p.value
    }
  } else if (test.method == "fisher") {
    cat("performing Fisher test\n")
    cbeta <- matrix(cut(sig, breaks = c(-1, median(sig), 99), label = c(0, 1)),
      nrow = nrow(sig)
    )
    ii <- which(apply(cbeta, 1, function(x) length(setdiff(unique(x), c("NA", NA))) > 1))
    for (i in ii) {
      pv[i] <- fisher.test(cbeta[i, ], class.label)$p.value
    }
  } else {
    stop("unknown test method")
  }

  ## qvalue

  qv <- rep(NA, nrow(sig))
  kk <- which(!is.na(pv))
  qv[kk] <- qvalue::qvalue(pv[kk])$qvalue

  ## return object
  rr <- data.frame(ID = rownames(sig))
  rr$Sig0 <- rowMeans(sig[, which(class.label == bb[1])], na.rm = TRUE)
  rr$Sig1 <- rowMeans(sig[, which(class.label == bb[2])], na.rm = TRUE)
  rr$DiffSig <- (rr$Sig1 - rr$Sig0)
  rr$P.Value <- pv
  rr$Q.Value <- qv
  colnames(rr) <- sub("Sig0", paste("AveSig.", bb[1], sep = ""), colnames(rr))
  colnames(rr) <- sub("Sig1", paste("AveSig.", bb[2], sep = ""), colnames(rr))

  ## order on absolute difference
  rr <- rr[which(rr$Q.Value < fdr), ]
  rr <- rr[order(-abs(rr$DiffSig)), ]


  return(rr)
}
#' @export
pgx.reduceCells <- function(counts, method, ncells, pheno = NULL, group.id = NULL) {
  if (ncol(counts) > ncells) {
    if (method == "pool") {
      message(">> Pooling cells...")


      pc <- pgx.poolCells(counts, ncells, groups = group.id)
      counts <- pc$counts
      if (!is.null(pheno)) {
        pheno1 <- apply(pheno, 2, function(x) {
          tapply(x, pc$cluster.id, function(aa) names(which.max(table(aa))))
        })
        pheno <- pheno1[colnames(counts), ]
        rownames(pheno) <- colnames(counts)
      }
    } else if (method == "subsample" && is.null(group.id)) {
      message(">> Subsampling cells...")
      sel <- sample(ncol(counts), ncells)
      counts <- counts[, sel]
      if (!is.null(pheno)) {
        pheno <- pheno[sel, ]
        rownames(pheno) <- colnames(counts)
      }
    } else if (method == "subsample" && !is.null(group.id)) {
      message(">> Subsampling cells in groups...")
      n1 <- round(ncells / length(unique(group.id)))
      sel <- tapply(1:ncol(counts), group.id, function(i) Matrix::head(sample(i), n1))
      sel <- unlist(sel)
      counts <- counts[, sel]
      if (!is.null(pheno)) {
        pheno <- pheno[sel, ]
        rownames(pheno) <- colnames(counts)
      }
    } else {
      stop("FATAL ERROR:: unknown method", method)
    }
  }
  list(counts = data.matrix(counts), pheno = data.frame(pheno))
}
#' GSEA Enrichment Plot (Up/Down Genes)
#'
#' Generates an enrichment plot for GSEA analysis with separate traces for up-regulated and down-regulated genes.
#'
#' @param rnk A numeric vector representing the ranked list.
#' @param gset.up A character vector specifying the gene set for up-regulated genes.
#' @param gset.dn A character vector specifying the gene set for down-regulated genes.
#' @param names A character vector of names for the ranked list.
#' @param main The main title of the plot.
#' @param decreasing Logical indicating whether the ranked list should be sorted in decreasing order.
#' @param cex.main The size of the main title text.
#' @param len.main The maximum number of characters for the main title before line breaks are applied.
#' @param sum.trace Logical indicating whether to plot the sum of the traces for up-regulated and down-regulated genes.
#' @param res The resolution of the plot.
#' @param xlab The label for the x-axis.
#' @param ylab The label for the y-axis.
#' @export
#' @return a plot
gsea.enplot.UPDN <- function(rnk, gset.up, gset.dn, names = NULL, main = NULL,
                             decreasing = TRUE, cex.main = 0.9, len.main = 40,
                             sum.trace = TRUE, res = 1200,
                             xlab = "Rank in ordered dataset",
                             ylab = "Ranked list metric") {
  if (!is.null(names)) names(rnk) <- names
  rnk <- rnk[!is.na(rnk)]
  rnk <- rnk[order(rnk + 1e-8 * rnorm(length(rnk)), decreasing = decreasing)]

  ## ranked list metric
  ii <- (1:length(rnk))
  if (length(ii) > res) ii <- ii[seq(1, length(ii), length(ii) / res)]
  qq <- quantile(rnk[ii], probs = c(0.01, 0.99), na.rm = TRUE)
  y1 <- qq[2]
  y0 <- qq[1]
  dy <- 0.25 * (y1 - y0)
  plot(ii, rnk[ii],
    type = "h", col = "grey", ylim = c(y0 - 1.25 * dy, y1),
    xlab = NA, ylab = ylab, xaxt = "n"
  )
  mtext(xlab, 1, line = 1.2, cex = 0.7)
  abline(h = 0, lty = 2, lwd = 0.5)

  ## gene set barcode
  jj <- match(gset.dn, names(rnk))
  arrows(jj, (y0 - dy), jj, y0 - 0.5 * dy, col = "grey30", lwd = 1, length = 0)
  jj <- match(gset.up, names(rnk))
  arrows(jj, (y0 - 0.5 * dy), jj, y0, col = "grey30", lwd = 1, length = 0)

  ## color legend
  kk <- c(seq(1, length(rnk) * 0.95, floor(length(rnk) / 10)), length(rnk))
  i <- 1
  for (i in 1:(length(kk) - 1)) {
    r <- mean(rnk[kk[c(i, i + 1)]])
    r1 <- (r / max(abs(rnk), na.rm = TRUE))
    r1 <- abs(r1)**0.66 * sign(r1)
    irnk <- round(10 + 10 * r1)
    cc <- gplots::bluered(20)[irnk]
    rect(kk[i], y0 - 1.30 * dy, kk[i + 1], y0 - 1.05 * dy, col = cc, border = NA)
  }

  ## weighted cumulative random walk (down genes)
  x0 <- 1 * (names(rnk) %in% gset.dn)
  x0 <- x0 * abs(rnk)
  n0 <- sum(!(names(rnk) %in% gset.dn))
  n1 <- sum(names(rnk) %in% gset.dn)
  r0 <- cumsum(x0 == 0) / sum(x0 == 0)
  r1 <- cumsum(x0) / (1e-4 + sum(x0))
  trace.dn <- (r1 - r0)
  trace.dn <- trace.dn / max(abs(trace.dn)) * 0.8
  if (max(trace.dn) >= abs(min(trace.dn))) trace.dn <- trace.dn * abs(y1)
  if (max(trace.dn) < abs(min(trace.dn))) trace.dn <- trace.dn * abs(y0)
  if (!decreasing) trace.dn <- -1 * trace.dn
  lines(ii, trace.dn[ii], col = "orange", type = "l", lwd = 2.4)

  ## weighted cumulative random walk
  x0 <- 1 * (names(rnk) %in% gset.up)
  x0 <- x0 * abs(rnk)
  n0 <- sum(!(names(rnk) %in% gset.up))
  n1 <- sum(names(rnk) %in% gset.up)
  r0 <- cumsum(x0 == 0) / sum(x0 == 0)
  r1 <- cumsum(x0) / (1e-4 + sum(x0))
  trace.up <- (r1 - r0)
  trace.up <- trace.up / max(abs(trace.up)) * 0.8
  if (max(trace.up) >= abs(min(trace.up))) trace.up <- trace.up * abs(y1)
  if (max(trace.up) < abs(min(trace.up))) trace.up <- trace.up * abs(y0)
  if (!decreasing) trace.up <- -1 * trace.up
  lines(ii, trace.up[ii], col = "green", type = "l", lwd = 2.4)

  if (sum.trace) lines(ii, trace.up[ii] + trace.dn[ii], col = "blue", type = "l", lwd = 2.4)

  if (is.null(main)) main <- "Enrichment plot"
  tt.main <- as.character(main)
  if (nchar(tt.main) > len.main) {
    tt.main <- paste(substr(tt.main, 1, len.main),
      substr(tt.main, len.main + 1, 2 * len.main),
      substr(tt.main, 2 * len.main + 1, 3 * len.main),
      sep = "\n"
    )
  }
  title(main = tt.main, cex.main = cex.main, line = 0.6)
}
#' @export
matrix.prod.SAVE <- function(matlist) {
  matx <- lapply(matlist, function(x) as.vector(x))
  matx <- as.matrix(do.call(cbind, matx))
  nna <- rowSums(!is.na(matx))
  matx <- apply(matx, 1, prod, na.rm = TRUE)**(1 / nna)
  xp <- matrix(matx, nrow = nrow(matlist[[1]]), ncol = ncol(matlist[[1]]))
  return(xp)
}
#' Make contrasts from group pairs
#'
#' This function creates a contrasts matrix from pairs of groups.
#'
#' @param main.group Character vector of the main groups in the comparisons.
#' @param ref.group Character vector of the reference groups in the comparisons.
#' @param groups Character vector of all groups. If NULL, groups are taken from main.group and ref.group.
#' @param comparisons Character vector of comparison names. If NULL, comparisons are generated from main.group and ref.group.
#'
#' @details The function first splits the main.group and ref.group into individual groups.
#' It then creates a contrasts matrix with rows for each group and columns for each comparison.
#' The matrix entries are 1 for groups in main.group, -1 for groups in ref.group, and 0 otherwise.
#'
#' @return A contrasts matrix with rows for groups and columns for comparisons.
#'
#' @export
makeContrastsFromPairs <- function(main.group, ref.group, groups = NULL, comparisons = NULL) {
  if (is.null(comparisons)) comparisons <- paste0(main.group, "_vs_", ref.group)
  main.group <- as.character(main.group)
  ref.group <- as.character(ref.group)
  groups1 <- unlist(strsplit(main.group, split = "[+]"))
  groups2 <- unlist(strsplit(ref.group, split = "[+]"))
  if (is.null(groups)) groups <- sort(unique(c(groups1, groups2)))
  contr.matrix <- matrix(0, nrow = length(groups), ncol = length(comparisons))
  i <- 1
  for (i in 1:length(comparisons)) {
    main <- strsplit(main.group[i], split = "[+]")[[1]]
    ref <- strsplit(ref.group[i], split = "[+]")[[1]]
    ct <- 1 * (groups %in% main) - 1 * (groups %in% ref)
    ct
    contr.matrix[, i] <- ct
  }
  colnames(contr.matrix) <- comparisons
  rownames(contr.matrix) <- groups

  ## normalize to zero mean and symmetric sum-to-one. Any NA to zero.
  for (i in 1:ncol(contr.matrix)) {
    m <- contr.matrix[, i]
    m[is.na(m)] <- 0
    contr.matrix[, i] <- 1 * (m > 0) / sum(m > 0) - 1 * (m < 0) / sum(m < 0)
  }

  return(contr.matrix)
}
#' Write data to a CLS file
#'
#' This function writes data to a CLS (CLustering Solutions) file format.
#' The CLS format is commonly used to store class labels or clustering assignments.
#'
#' @param y The vector of class labels or clustering assignments.
#' @param file The file path to save the CLS file.
#' @param name (Optional) The name of the CLS file.
#'
#' @export
#'
#' @return None
#'
write.cls <- function(y, file, name = "") {
  ## prepare class file
  is.numeric <- (length(setdiff(unique(y), NA)) > 3)
  if (is.numeric) {
    write("#numeric", file = file, append = FALSE)
    suppressWarnings(write(paste0("#", name), file = file, append = TRUE))
    suppressWarnings(write(paste(y, collapse = " "), file = file, append = TRUE))
  } else {
    write(paste(length(y), "2 1", sep = " "), file = file, append = FALSE)
    cls.types <- paste("#", paste(unique(y), collapse = " "), sep = "")
    suppressWarnings(write(cls.types, file = file, append = TRUE))
    suppressWarnings(write(paste(y, collapse = " "), file = file, append = TRUE))
  }
}
#' Collapse duplicate genes in NGS dataset
#'
#' @param ngs An NGS dataset object
#'
#' @return The NGS dataset object with collapsed genes
#'
#' @details
#' This function collapses duplicate gene names in an NGS dataset by summing their counts.
#'
#' It first identifies genes that occur only once or more than once in the dataset.
#' For uniquely occurring genes, it keeps them as is.
#' For duplicate genes, it sums their counts by gene name.
#'
#' The gene and count matrices in the NGS dataset object are then updated accordingly.
#'
#' @export
ngs.collapseByGene <- function(ngs) {
  gene <- as.character(ngs$genes$gene_name)
  p1 <- names(which(table(gene) == 1))
  p2 <- names(which(table(gene) > 1))
  if (length(p2) == 0) {
    gene <- as.character(ngs$genes$gene_name)
    rownames(ngs$genes) <- gene
    rownames(ngs$counts) <- gene
    return(ngs)
  }
  j1 <- which(gene %in% p1)
  j2 <- which(gene %in% p2)
  x1 <- ngs$counts[j1, , drop = FALSE]
  rownames(x1) <- ngs$genes$gene_name[j1]
  x2 <- apply(ngs$counts[j2, , drop = FALSE], 2, function(x) tapply(x, gene[j2], sum))
  if (length(p2) == 1) {
    x2 <- matrix(x2, nrow = 1)
    rownames(x2) <- p2
    colnames(x2) <- colnames(x1)
  }
  x1 <- rbind(x1, x2)
  x1 <- x1[!(rownames(x1) %in% c(NA, "", "NA")), , drop = FALSE]
  x1 <- x1[order(rownames(x1)), , drop = FALSE]
  ngs$genes <- ngs$genes[match(rownames(x1), ngs$genes$gene_name), ]
  ngs$counts <- x1
  rownames(ngs$genes) <- rownames(ngs$counts) <- rownames(x1)
  return(ngs)
}
#' @title Extract gene-specific subgraph from PubMed network
#'
#' @param gr An igraph network object generated from PubMed data
#' @param gene Character vector of gene symbols to extract network for
#' @param nmin Minimum number of genes in connected components to keep
#'
#' @return An igraph object containing the gene-specific subgraph
#'
#' @description Extracts a subgraph related to the specified genes from a PubMed network graph.
#'
#' @details This function takes a PubMed network graph \code{gr} and a vector of \code{gene} symbols.
#' It extracts the nodes containing those genes and the edges between them.
#' Connected components with fewer than \code{nmin} genes are removed.
#'
#' Edge weights and gene annotations are added to the subgraph.
#'
#' @export
pmid.extractGene <- function(gr, gene, nmin = 3) {
  jj <- c()
  for (g in gene) {
    j1 <- which(sapply(igraph::V(gr)$genes, function(s) (gene %in% s)))
    jj <- c(jj, j1)
  }
  jj <- unique(jj)
  ngene1 <- sapply(igraph::V(gr)$genes[jj], length)
  gr1 <- igraph::induced_subgraph(gr, jj)
  if (verbose > 0) cat("annotating edges...\n")
  gr1 <- pmid.annotateEdges(gr1)
  nshared <- unlist(parallel::mclapply(igraph::E(gr1)$genes, length, mc.cores = NCORE()))
  ee <- which(!(nshared == 1 & sapply(igraph::E(gr1)$genes, "[", 1) %in% gene))
  gr1 <- igraph::subgraph.edges(gr1, ee, delete.vertices = TRUE)
  cmp <- igraph::components(gr1)
  jj <- which(cmp$membership %in% which(cmp$csize >= nmin))
  gr1 <- igraph::induced_subgraph(gr1, jj)
  return(gr1)
}
#' Scatterplot matrix using ggplot2
#'
#' @param F A data frame or matrix.
#' @param F2 Optional second data frame or matrix to plot on y-axis. Default is NULL.
#' @param title_cex Title text size. Default is 2.
#' @param no.axes Hide axis labels and tick marks. Default is FALSE.
#' @param ... Other arguments passed to ggplot2 functions.
#'
#' @return A ggplot scatterplot matrix object
#'
#' @title Scatterplot matrix using ggplot2
#'
#' @description Generate a scatterplot matrix from a data frame or matrix using ggplot2.
#'
#' @details This function takes a data frame or matrix \code{F} and generates a scatterplot matrix.
#' Each variable in \code{F} is plotted against each other variable. An optional second data frame
#' or matrix \code{F2} can be provided for the y-axis variables. Titles and axis labels are added
#' automatically. \code{title_cex} controls the title size. \code{no.axes} hides axes.
#'
#' @export
plot_ggsplom <- function(F, F2 = NULL, title_cex = 2, no.axes = FALSE, ...) {
  if (is.null(F2)) {
    F2 <- F
  }
  lim0 <- range(F)
  lim1 <- range(F2)
  symm <- all(colnames(F) == colnames(F2))

  gg <- intersect(rownames(F), rownames(F2))
  F <- F[gg, ]
  F2 <- F2[gg, ]

  plt <- list()
  k <- 1
  for (i in 1:ncol(F)) {
    for (j in 1:ncol(F2)) {
      if (i == j && symm) {
        p1 <- plot_ggscatter(0, 0, cex = 0) +
          ggplot2::theme_void() +
          ggplot2::geom_text(x = 0, y = 0, label = colnames(F)[i], size = title_cex)
      } else {
        p1 <- plot_ggscatter(F[, i], F[, j]) +
          ggplot2::xlim(lim0) + ggplot2::ylim(lim1) +
          ggplot2::xlab(colnames(F)[i]) +
          ggplot2::ylab(colnames(F2)[j])
      }
      plt[[k]] <- p1
      k <- k + 1
    }
  }

  blankx <- ggplot2::theme(
    plot.margin = ggplot2::margin(0, 0, 0, 0),
    axis.text.x = ggplot2::element_blank(),
    axis.title.x = ggplot2::element_blank(),
    axis.ticks.x = ggplot2::element_blank()
  )
  blanky <- ggplot2::theme(
    plot.margin = ggplot2::margin(0, 0, 0, 0),
    axis.text.y = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_blank(),
    axis.ticks.y = ggplot2::element_blank(),
  )

  np <- length(plt)
  nc <- ncol(F)
  nr <- ncol(F2)
  border.y <- ((0:(np - 1) %% nc) == 0)
  border.x <- ((0:(np - 1) %/% nc) == (nr - 1))
  for (i in 1:np) {
    if (no.axes || !border.y[i]) plt[[i]] <- plt[[i]] + blanky
    if (no.axes || !border.x[i]) plt[[i]] <- plt[[i]] + blankx
  }
  cowplot::plot_grid(
    plotlist = plt, labels = NA,
    nrow = nr, ncol = nc,
    rel_widths = c(1.2, rep(1, nc - 1)),
    rel_heights = c(rep(1, nr - 1), 1.2)
  )
}
#' @export
prot.calcCopyNumber <- function(data, mol, y) {
  ## data should be a numeric matrix
  ## mol is the molecular weight
  ## y is the mass of the cell in PICOGRAM (pg)
  TotalIntenisty <- apply(data, 2, sum)
  ## mass in gram
  zz <- vector(length = 0)
  Mass <- vector(length = 0)
  MassProtein <- for (i in 1:length(TotalIntenisty)) {
    zz <- (data[, i] * y) / TotalIntenisty[i] * 10^-12
    Mass <- cbind(Mass, zz)
  }
  colnames(Mass) <- colnames(data)
  # calculate moles
  Mol <- Mass / (mol * 1000)
  Copy <- Mol * 6.022140857 * 10^23
  return(Copy)
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
#' @describeIn read.csv3 Read delimited table automatically determine separator
#' @export
read.csv3.BAK <- function(file, ...) {
  ## read delimited table automatically determine separator
  line1 <- as.character(read.csv(file, comment.char = "#", sep = "\n", nrow = 1)[1, ])
  sep <- names(which.max(sapply(c("\t", ",", ";"), function(s) length(strsplit(line1, split = s)[[1]]))))
  message("[read.csv3] sep = ", sep)
  read.csv(file, comment.char = "#", sep = sep, ...)
}
#' @title Make specific contrasts
#'
#' @param df Data frame with sample metadata
#' @param contrasts Character vector of contrast names to extract
#' @param mingrp Minimum group size for inclusion in contrasts
#'
#' @return Data frame with specified contrast columns
#'
#' @description This function extracts specific named contrasts from a sample metadata data frame.
#' It searches for contrast names matching the input patterns and returns those contrast columns.
#'
#' @details Contrast names are searched for patterns like "Group1:Group2_vs_Group3:Group4".
#' The groups before and after the "_vs_" are extracted.
#' Dummy variable columns are generated based on these groups.
#' Groups smaller than mingrp are excluded.
#'
#' The final output is a data frame containing only the specified contrast columns,
#' with dummy variable encodings of 0/1/-1 for the sample groups.
#'
#' @export
pgx.makeSpecificContrasts <- function(df, contrasts, mingrp = 3) {
  K <- c()
  i <- 1
  for (i in 1:length(contrasts)) {
    ct <- contrasts[i]
    if (!(grepl("[:]", ct) && grepl("_vs_", ct))) next()
    ph <- sub("[:].*", "", ct)
    if (!ph %in% colnames(df)) next()

    groups <- strsplit(sub(".*[:]", "", ct), split = "_vs_")[[1]]
    y <- -1 * (df[, ph] == groups[2]) + 1 * (df[, ph] == groups[1])
    K <- cbind(K, y)
    colnames(K)[ncol(K)] <- ct
  }
  rownames(K) <- rownames(df)

  K0 <- contrastAsLabels(K)
  group <- pgx.getConditions(K0)
  if (length(levels(group)) > 0.5 * nrow(K)) {
    cat("WARNING:: contrast matrix looks degenerate. consider removing a contrast.\n")
  }

  contr.matrix <- K[which(!duplicated(group)), , drop = FALSE]
  rownames(contr.matrix) <- group[which(!duplicated(group))]
  res <- list(contr.matrix = contr.matrix, group = group)

  return(res)
}
#' @export
pgx.fastq2counts <- function(fastq_dir, destdir, indexdir, nthread = 4, do.qc = FALSE,
                             species = c("human", "mouse", "rat"),
                             quant.method = "salmon",
                             trimming = TRUE, trimmethod = "trimmomatic",
                             instrument = "HiSeq", library_layout = "SINGLE") {
  species <- species[1]

  run_fastqc1 <- function(destdir, fastq_dir, nthread) {
    cat(paste("Running FastQC... ", Sys.time(), "\n", sep = ""))
    dir.exists(paste0(destdir, "/fastqc"))
    if (!dir.exists(paste0(destdir, "/fastqc"))) {
      system(paste0("mkdir -p ", destdir, "/fastqc"))
    } else {
      system(paste0("rm -fr ", destdir, "/fastqc"))
      system(paste0("mkdir -p ", destdir, "/fastqc"))
    }
    fastq_files <- list.files(fastq_dir, pattern = "[.]fastq$", full.names = TRUE)

    cmd <- paste0(
      "fastqc -o ", destdir, "/fastqc/ --threads ",
      nthread, " ", paste(fastq_files, collapse = " ")
    )
    system(cmd)
  }

  run_fastqc2 <- function(destdir, fastq_dir, nthread) {
    cat(paste("Running FastQC... ", Sys.time(), "\n", sep = ""))
    dir.exists(paste0(destdir, "/fastqc_trimmed"))
    if (!dir.exists(paste0(destdir, "/fastqc_trimmed"))) {
      system(paste0("mkdir -p ", destdir, "/fastqc_trimmed"))
    } else {
      system(paste0("rm -fr ", destdir, "/fastqc_trimmed"))
      system(paste0("mkdir -p ", destdir, "/fastqc_trimmed"))
    }
    fastq_files <- list.files(fastq_dir, pattern = "_trimmed.fastq$", full.names = TRUE)
    cmd <- paste0(
      "fastqc -o ", destdir, "/fastqc_trimmed/ --threads ",
      nthread, " ", paste(fastq_files, collapse = " ")
    )
    system(cmd)
  }

  ## ----------- Unzip FASTQ files (if necessary)
  dir(fastq_dir)
  dir(fastq_dir, "[.]gz$", include.dirs = FALSE)
  if (length(dir(fastq_dir, "[.]gz$", include.dirs = FALSE))) {
    cat(">>> Unzipping FASTQ files... \n")

    cmd <- paste0("(cd ", fastq_dir, " && gunzip *gz)")
    cmd
    system(cmd)
  }

  ## remove any previously trimmed sequences
  system(paste0("rm -f ", fastq_dir, "/*_trimmed.fastq"))
  system(paste0("rm -f ", fastq_dir, "/*trimmed*"))
  system(paste0("rm -f ", fastq_dir, "/*trimming*"))

  ## ----------- Run FastQC on each fastq file to generate quality control (QC) reports.
  if (do.qc) {
    cat(">>> Running FastQC ... \n")
    run_fastqc1(destdir = destdir, fastq_dir = fastq_dir, nthread = nthread)
  }

  ## ----------- Run Trimmomatic (some say it is not necessary for salmon...)
  file_id <- sub(".fastq$", "", dir(fastq_dir, pattern = ".fastq$"))
  file_id <- grep("trimmed$", file_id, value = TRUE, invert = TRUE)
  file_id

  if (trimming && trimmethod == "trimmomatic") {
    cat(">>> Running Trimmomatic... \n")
    i <- 1
    for (i in 1:length(file_id)) {
      GREP2::trim_fastq(
        srr_id = file_id[i], fastq_dir = fastq_dir,
        instrument = instrument, library_layout = library_layout,
        destdir = fastq_dir, n_thread = nthread
      )
    }
  }
  if (trimming && trimmethod == "trimgalore") {
    cat(">>> Running TrimGalore... \n")
    fastq1 <- fastq2 <- NULL
    if (library_layout == "SINGLE") {
      fastq1 <- file.path(fastq_dir, paste0(file_id, ".fastq"))
      fastq2 <- NULL
    } else {
      fastq1 <- file.path(fastq_dir, paste0(file_id, "_1.fastq"))
      fastq2 <- file.path(fastq_dir, paste0(file_id, "_2.fastq"))
    }

    trimgalore_fastq(fastq1,
      fastq2 = fastq2,
      adapter1 = NULL, adapter2 = NULL,
      dest.dir = fastq_dir,
      threads = nthread,
      trimgalore = "trim_galore"
    )


    fq <- dir(fastq_dir, pattern = "_trimmed.fq$", full.names = TRUE)
    fq
    for (f in fq) file.rename(from = f, to = sub("fq$", "fastq", f))
  }

  if (quant.method == "salmon") {
    ## ----------- Before running Salmon, you will have to build index first.
    dir.exists(indexdir)
    index2 <- file.path(indexdir, paste0(species, "_transcripts_release99_index"))
    index2
    dir.exists(index2)
    if (!dir.exists(index2)) {
      cat(">>> Building index for species", species, "... \n")
      system(paste("mkdir -p", indexdir))
      # The f build_index likely comes from package GREP2, but better confirm
      build_index(species = species, kmer = 31, ens_release = 99, destdir = indexdir)
    } else {
      cat(">>> Found index folder at ", index2, "\n")
    }
    dir(indexdir)

    ## ----------- Run Salmon
    cat(">>> Running Salmon... \n")
    i <- 1
    for (i in 1:length(file_id)) {
      pgx.run_salmon(
        srr_id = file_id[i], library_layout = library_layout,
        index_dir = index2, destdir = destdir,
        fastq_dir = fastq_dir, use_trimmed_fastq = trimming,
        other_opts = "--validateMappings", nthread = nthread
      )
    }
  }
  if (quant.method == "kallisto") {
    ## ----------- Before running Kallisto, you will have to build index first.
    indexdir
    dir.exists(indexdir)
    kallisto_index <- file.path(indexdir, paste0(species, "_transcripts.idx"))
    if (!file.exists(kallisto_index)) {
      cat(">>> Building index for species", species, "... \n")
      system(paste("mkdir -p", indexdir))


      if (species == "human") ref.genome <- file.path(indexdir, "Homo_sapiens.GRCh38.cdna.all.fa.gz")
      if (species == "mouse") ref.genome <- file.path(indexdir, "Mus_musculus.GRCm38.cdna.all.fa.gz")
      kallisto_index
      ref.genome
      system(paste("kallisto index -i", kallisto_index, ref.genome))
    } else {
      cat(">>> Found index at ", kallisto_index, "\n")
    }

    ## ----------- Run Kallisto
    cat(">>> Running Kallisto... \n")
    i <- 1
    fq <- dir(fastq_dir, pattern = ".fastq$", full.names = TRUE)
    fq <- grep("trimmed", fq, invert = TRUE, value = TRUE)
    if (trimming) fq <- dir(fastq_dir, pattern = "_trimmed.fastq$", full.names = TRUE)
    kallisto_outdir <- file.path(destdir, "kallisto")
    system(paste("mkdir -p", kallisto_outdir))

    cmd <- "kallisto quant"
    if (library_layout == "SINGLE") cmd <- "kallisto quant --single -l 250 -s 50" ## NEED RETHINK!!!
    cmd <- paste(cmd, "-i", kallisto_index, "-b 100 -t 32")
    cmd
    i <- 1
    for (i in 1:length(fq)) {
      f1 <- gsub(".*/|.fastq$|_trimmed.fastq$", "", fq[i])
      out1 <- file.path(kallisto_outdir, f1)
      cmd1 <- paste(cmd, "-o", out1, fq[i])
      cmd1
      system(cmd1)
    }
  }

  ## ----------- Run MultiQC
  if (do.qc) {
    cat(">>> Running MultiQC... \n")
    # The f run_multiqc likely comes from package GREP2, but better confrim
    run_multiqc(
      fastqc_dir = file.path(destdir, "fastqc"),
      salmon_dir = file.path(destdir, quant.method),
      destdir = destdir
    )
  }

  ## ----------- Run tximport
  if (quant.method == "salmon") {
    cat(">>> Running TxImport on Salmon files... \n")
    # The f run_tximport likely comes from package GREP2, but better confrim
    txi <- run_tximport(
      srr_id = file_id,
      species = species,
      salmon_dir = paste0(destdir, "/salmon"),
      countsFromAbundance = "lengthScaledTPM"
    )
  }

  if (quant.method == "kallisto") {
    cat(">>> Running TxImport on Kallisto files... \n")
    txi <- run_tximport_kallisto(
      srr_id = file_id,
      species = species,
      kallisto_dir = paste0(destdir, "/kallisto")
    )
  }

  ## ----------- Extract counts
  genes <- txi$gene_counts[, 2:3]
  counts <- as.matrix(txi$gene_counts[, 4:ncol(txi$gene_counts), drop = FALSE])
  rownames(genes) <- rownames(counts) <- txi$gene_counts$ENSEMBL

  ## ----------- Collapse multiple probes to single gene by summing up counts
  counts <- apply(counts, 2, function(x) tapply(x, genes$SYMBOL, sum))
  genes <- genes[match(rownames(counts), genes$SYMBOL), ]
  rownames(genes) <- genes$SYMBOL
  counts <- round(counts, digits = 3)

  ## ----------- Save counts as text files for Omics Playground
  colnames(genes) <- c("gene_name", "gene_title")
  system(paste0("mkdir -p ", destdir, "/files_csv"))
  cat(">>> Writing CSV files to", paste0(destdir, "/files_csv/ \n"))
  write.csv(counts, file = paste0(destdir, "/files_csv/counts.csv"))
  write.csv(genes, file = paste0(destdir, "/files_csv/genes.csv"))
  samples <- data.frame(sample = colnames(counts), group = NA, phenotype1 = NA) ## empty template
  write.csv(samples, file = paste0(destdir, "/files_csv/samples.csv"), row.names = FALSE)

  cat(">>> done!\n")
}
#' @export
kcomp <- function(g, k = 3, minsize = 0) {
  ## Return  K-largest components
  cmp <- igraph::components(g)
  sort(cmp$csize)
  cmp.size <- cmp$csize
  names(cmp.size) <- 1:length(cmp.size)
  cmp.size <- cmp.size[which(cmp.size >= minsize)]
  cmp.top <- names(cmp.size)
  if (!is.null(k)) cmp.top <- Matrix::head(names(sort(-cmp.size)), k)
  vtx <- which(cmp$membership %in% cmp.top)
  igraph::induced_subgraph(g, igraph::V(g)[vtx])
}
#' Find Louvain clusters using shared nearest neighbor graph
#'
#' @title Find Louvain clusters using shared nearest neighbor graph
#'
#' @param expr Expression matrix with genes in rows and samples in columns
#' @param k Number of nearest neighbors for graph construction
#' @param resolution Resolution parameter for Louvain clustering
#'
#' @return A vector of Louvain cluster assignments
#'
#' @description
#' Clusters samples using a shared nearest neighbor graph and Louvain clustering.
#'
#' @details
#' This function constructs a shared nearest neighbor graph from the expression matrix,
#' where samples are connected if they appear in each other's k-nearest neighbor lists.
#'
#' Louvain clustering is then applied to detect communities in the graph. The resolution
#' parameter controls the number and size of clusters.
#'
#' The output is a vector of cluster assignments for each sample.
#'
#' @export
pgx.findLouvainClusters.SNN <- function(X, prefix = "c", level = 1, gamma = 1, small.zero = 0.01) {
  ## find clusters using graph clustering method
  message("perform Louvain clustering...")


  if (level == 1) {
    suppressMessages(suppressWarnings(
      gr <- scran::buildSNNGraph(t(X), d = 50)
    ))
  } else {
    ## finer clusters
    suppressMessages(suppressWarnings(
      gr <- scran::buildSNNGraph(t(X), d = 50, k = 2)
    ))
  }

  idx <- igraph::cluster_louvain(gr)$membership
  idx <- paste0(prefix, idx)

  if (!is.null(idx) && small.zero > 0) {
    ## ------------ zap small clusters to "0"
    min.size <- pmax(3, small.zero * length(idx))
    small.clusters <- names(which(table(idx) < min.size))
    idx[which(idx %in% small.clusters)] <- "0"
  }

  ## rename levels with largest cluster first
  idx <- factor(idx, levels = names(sort(-table(idx))))
  levels(idx) <- paste0(prefix, 1:length(levels(idx)))
  idx <- as.character(idx)
  message("Found ", length(unique(idx)), " clusters...")
  return(idx)
}
#' @title Expand Annotation Matrix
#'
#' @description This function expands an annotation matrix by
#' converting categorical variables into dummy variables and ranking numerical variables.
#'
#' @param A A data frame representing the annotation matrix to be expanded.
#'
#' @details The function takes an annotation matrix `A` as input,
#' where each column represents a variable and each row represents an
#' observation. For each column, if the variable is numerical, the function
#' ranks the values and stores them in a new column. If the variable is
#' categorical, the function creates dummy variables for each level of the c
#' ategorical variable and stores them in new columns.
#'
#' The resulting expanded annotation matrix is returned as a data frame,
#' where each row represents an observation and each column represents a
#' variable or a level of a categorical variable.
#'
#' @return A data frame representing the expanded annotation matrix.
#'
#' @export
#' @export
expandAnnotationMatrixSAVE <- function(A) {
  ## get expanded annotation matrix
  nlevel <- apply(A, 2, function(x) length(unique(x)))
  y.isnum <- apply(A, 2, is.num)

  i <- 1
  m1 <- list()
  for (i in 1:ncol(A)) {
    if (is.num(A[, i])) {
      m0 <- matrix(rank(A[, i], na.last = "keep"), ncol = 1)
      colnames(m0) <- colnames(A)[i]
    } else {
      x <- as.character(A[, i])
      x[is.na(x)] <- "_"
      m0 <- model.matrix(~ 0 + x)
      colnames(m0) <- sub("^x", "", colnames(m0))
    }
    if (NCOL(m0) > 1) {
      colnames(m0) <- paste0(colnames(A)[i], "=", colnames(m0))
    }
    m1[[i]] <- m0
  }
  names(m1) <- colnames(A)

  M <- do.call(cbind, m1)
  rownames(M) <- rownames(A)
  return(M)
}
#' Normalize a contrasts matrix
#'
#' @param contr.matrix A contrasts matrix with groups as rows and comparisons as columns
#'
#' @return The normalized contrasts matrix
#'
#' @description
#' Normalizes a contrasts matrix to have zero row means and symmetric column sums.
#'
#' @details
#' This function normalizes a contrasts matrix in the following way:
#'
#' - Subtracts the row means, centering each contrast around 0
#' - Divides each column by the absolute column sum, making columns symmetric sum-to-1
#' - Replaces any NA values with 0
#'
#' The result is a normalized contrasts matrix suitable for statistical testing.
#'
#' @export
normalizeContrasts <- function(contr.matrix) {
  ## normalize to zero mean and symmetric sum-to-one. Any NA to zero.
  for (i in 1:ncol(contr.matrix)) {
    m <- contr.matrix[, i]
    m[is.na(m)] <- 0
    contr.matrix[, i] <- 1 * (m > 0) / sum(m > 0) - 1 * (m < 0) / sum(m < 0)
  }
  contr.matrix
}
#' @export
pgx.getTCGAproteomics <- function() {
  GenomicDataCommons::status()

  qfiles <- GenomicDataCommons::files() %>% filter(~ cases.project.project_id == "TCGA-BRCA" &
    type == "gene_expression" &
    analysis.workflow_type == "HTSeq - Counts")
  manifest_df <- qfiles %>% GenomicDataCommons::manifest()

  fnames <- GenomicDataCommons::gdcdata(manifest_df$id[1:2], progress = FALSE)

  resp <- GenomicDataCommons::cases() %>%
    filter(~ project.project_id == "TCGA-BRCA" &
      project.project_id == "TCGA-BRCA") %>%
    GenomicDataCommons::facet("samples.sample_type") %>%
    GenomicDataCommons::aggregations()
  resp$samples.sample_type

  resp <- GenomicDataCommons::cases() %>%
    filter(~ project.project_id == "CPTAC-3" &
      project.project_id == "CPTAC-3") %>%
    GenomicDataCommons::facet("samples.sample_type") %>%
    GenomicDataCommons::aggregations()
  resp$samples.sample_type
}
#' Convert GMT to Matrix
#'
#' This function converts a GMT (Gene Matrix Transposed) file into a sparse matrix representation.
#'
#' @param gmt A list containing gene sets in GMT format.
#' @param bg A character vector specifying the background genes to consider. If NULL, the most frequently occurring genes in the GMT will be used as the background.
#' @param use.multicore Logical indicating whether to use parallel processing for faster computation.
#'
#' @return A sparse matrix representation of the GMT.
#' @export
gmt2mat.nocheck <- function(gmt, bg = NULL, use.multicore = TRUE) {
  if (is.null(bg)) {
    bg <- names(sort(table(unlist(gmt)), decreasing = TRUE))
  }
  gmt <- lapply(gmt, function(s) intersect(bg, s))
  j <- 1
  idx <- c()
  if (use.multicore) {
    idx <- parallel::mclapply(gmt, function(s) match(s, bg))
    idx[sapply(idx, length) == 0] <- 0
    idx <- sapply(1:length(idx), function(i) rbind(idx[[i]], i))
    idx <- matrix(unlist(idx[]), byrow = TRUE, ncol = 2)
    idx <- idx[!is.na(idx[, 1]), ]
    idx <- idx[idx[, 1] > 0, ]
  } else {
    idx <- c()
    for (j in 1:length(gmt)) {
      ii0 <- which(bg %in% gmt[[j]])
      if (length(ii0) > 0) {
        idx <- rbind(idx, cbind(ii0, j))
      }
    }
  }
  D <- Matrix::sparseMatrix(idx[, 1], idx[, 2],
    x = rep(1, nrow(idx)),
    dims = c(length(bg), ncol = length(gmt))
  )
  dim(D)
  rownames(D) <- bg
  colnames(D) <- names(gmt)
  D
}
#' Update PGX-info file with new pgx object (WIP)
#'
#' @export
pgxinfo.addPgx <- function(pgx.dir, pgx, file = "datasets-info.csv",
                           update.fc = TRUE) {
  info.file <- file.path(pgx.dir, file)
  pgxinfo <- read.csv(info.file, row.names = 1)
  pgxinfo <- pgx.updateInfoPGX(pgxinfo, pgx, remove.old = TRUE)
  write.csv(pgxinfo, file = info.file)

  allfc.file <- file.path(pgx.dir, "datasets-allFC.csv")
  tsne.file <- file.path(pgx.dir, "datasets-tsne.csv")
  h5.file <- file.path(pgx.dir, "datasets-sigdb.h5")
  pgxname <- sub("[.]pgx$", "", pgx$name)
  pgxname

  ## add FC columns to allFC file
  if (update.fc && !file.exists(allfc.file)) {
    ## complete update of all files
    pgx.initDatasetFolder(pgx.dir)
    return()
  }

  ## add FC columns to allFC file
  allfC <- NULL
  if (update.fc && file.exists(allfc.file)) {
    allFC <- data.table::fread(allfc.file, check.names = FALSE) ## HEADER!!!
    allFC <- as.matrix(allFC, rownames = 1)
    del <- grep(paste0("\\[", pgxname, "\\]"), colnames(allFC))
    if (length(del) > 0) {
      allFC <- allFC[, -del, drop = FALSE]
    }
    F1 <- pgx.getMetaMatrix(pgx)$fc
    F1 <- F1[match(rownames(allFC), rownames(F1)), ]
    colnames(F1) <- paste0("[", pgxname, "] ", colnames(F1))
    allFC <- cbind(allFC, F1)
    allFC <- round(allFC, digits = 4)
    allFC1 <- data.frame(gene = rownames(allFC), allFC, check.names = FALSE)
    data.table::fwrite(allFC1, allfc.file) ## HEADER!!!

    ## NEED RETHINK!!!! This could be done perhaps more efficient
    ## when updating with one extra dataset.
    pgx.createSignatureDatabaseH5.fromMatrix(h5.file,
      X = allFC,
      update.only = FALSE
    )

    ## update t-SNE file (in the future we do not need this file)
    ## Just get the tSNE from sigdb
    cn <- rhdf5::h5read(h5.file, "data/colnames")
    tsne <- rhdf5::h5read(h5.file, "clustering/tsne2d")
    rownames(tsne) <- cn
    colnames(tsne) <- paste0("tsne.", 1:ncol(tsne))
    write.csv(tsne, file = tsne.file)
  }
}
#' @export
pgx.readMatrixH5 <- function(h5.file, select = NULL, rows = NULL) {
  if (is.null(select) && is.null(rows)) {
    X <- rhdf5::h5read(h5.file, "data/matrix")
    rn <- rhdf5::h5read(h5.file, "data/rownames")
    cn <- rhdf5::h5read(h5.file, "data/colnames")
  }
  if (!is.null(select) || !is.null(rows)) {
    X <- rhdf5::h5read(h5.file, "data/matrix", index = list(rows, select))
    rn <- rhdf5::h5read(h5.file, "data/rownames", index = list(rows))
    cn <- rhdf5::h5read(h5.file, "data/colnames", index = list(select))
  }
  rownames(X) <- rn
  colnames(X) <- cn
  X[which(X < -999999)] <- NA
  as.matrix(X)
}
#' @export
prot.normalizeCounts <- function(counts, scale = 1e6, scaling = "by.column",
                                 qnormalize = TRUE, zero.offset = 0.01,
                                 zero.thr = 0.25) {
  ## ------------------------------------------------------------
  ## start with original counts from MaxQuant
  ## ------------------------------------------------------------
  X <- counts

  X[X <= zero.thr] <- NA ## treat zero as NA for the moment
  which.zeros <- Matrix::which(X == 0 | is.na(X), arr.ind = TRUE)

  ## ------------------------------------------------------------
  ## normalize to CPM (or otherwise)
  ## ------------------------------------------------------------

  if (is.null(scale)) {
    scale <- mean(Matrix::colSums(X, na.rm = TRUE))
    message("set scale parameter to = ", scale)
  }

  if (scaling == "by.column") {
    message("scaling by columns to ", scale, " total counts")
    X <- t(t(X) / Matrix::colSums(X, na.rm = TRUE)) * scale ## counts-per-million
  } else if (scaling == "global") {
    message("global average scaling to ", scale, " counts")
    X <- X / mean(Matrix::colSums(X, na.rm = TRUE)) * scale ## counts-per-million
  } else {
    stop("FATAL:: unknown scaling method. scaling = ", scaling)
  }
  Matrix::colSums(X, na.rm = TRUE)

  ## ------------------------------------------------------------
  ## set zero offset
  ## ------------------------------------------------------------

  if (zero.offset > 0 && zero.offset < 1) {
    q0 <- quantile(X, probs = zero.offset, na.rm = TRUE)
    log2(q0)
    message("set log2(x)=0 to ", log2(q0))
    X <- X / q0 ## scale to unit
  }

  ## ------------------------------------------------------------
  ## choose quantile or median normalization (on not-missing values)
  ## ------------------------------------------------------------
  if (qnormalize) {
    jj <- sample(ncol(X), 10)
    X <- limma::normalizeQuantiles(X)
    X <- pmax(X, 0)
  }

  ## ------------------------------------------------------------
  ## put back 'missing' values to zero
  ## ------------------------------------------------------------
  X[which.zeros] <- 0

  return(X)
}
#' @export
silac.ttest <- function(X, group1, group2, method = "limma") {
  out1 <- NULL ## important to avoid old results
  p1 <- p2 <- p3 <- NULL ## just for testing
  if (method == "geiger") {
    ## Your code
    out1 <- volcano(X[, group1], X[, group2], rownames(X))
    out1$Gene <- NULL
    colnames(out1) <- c("logFC", "P.Value")
    p1 <- out1
  } else if (method == "genefilter") {
    ## faster t-test

    jj <- which(colnames(X) %in% c(group1, group2))
    y <- 1 * (colnames(X)[jj] %in% group1) + 2 * (colnames(X)[jj] %in% group2)
    out1 <- genefilter::rowttests(X[, jj], factor(y))
    out1$statistic <- NULL
    colnames(out1) <- c("logFC", "P.Value")
    p2 <- out1
  } else if (method == "limma") {
    ## See e.g. https://bioconductor.org/help/course-materials/2010/BioC2010/limma2.pdf

    jj <- which(colnames(X) %in% c(group1, group2))
    y <- 1 * (colnames(X)[jj] %in% group1)
    design <- model.matrix(~y)
    colnames(design) <- c("Intercept", "group1_vs_group2")
    fit <- limma::eBayes(limma::lmFit(X[, jj], design))
    out1 <- limma::topTable(fit, coef = 2, sort.by = "none", number = Inf)
    out1 <- out1[, c("logFC", "P.Value", "adj.P.Val")]
    p3 <- out1
  } else {
    stop("ERROR:: unknown method")
  }


  ## set NaN p.values to 1??
  out1$P.Value[which(is.nan(out1$P.Value))] <- 1

  ## add FDR q.value
  if (!("adj.P.Val" %in% names(out1))) {
    out1$adj.P.Val <- p.adjust(out1$P.Value, method = "fdr")
  }

  ## good practice to add group means
  avg <- cbind(rowMeans(X[, group1], na.rm = TRUE), rowMeans(X[, group2], na.rm = TRUE))

  colnames(avg) <- paste0("avg.", c("group1", "group2"))
  out1 <- cbind(out1, avg)

  return(out1)
}
#' @export
matrix.mean <- function(...) {
  matlist <- list(...)
  if (inherits(matlist[[1]], "list")) matlist <- matlist[[1]]
  matnna <- lapply(matlist, function(x) 1 * !is.na(x))
  p <- Reduce("+", matnna)
  matlist <- lapply(matlist, function(x) {
    x[is.na(x)] <- 0
    x
  })
  xsum <- Reduce("+", matlist) / p
  return(xsum)
}
#' @param x A numeric vector of values to plot.
#'
#' @title Plot Histogram with ggplot2
#'
#' @description Generate a histogram using ggplot2.
#'
#' @details This function takes a numeric vector \code{x} and generates a histogram
#' plot using ggplot2. The data is passed to \code{ggplot(data=x)} to generate a histogram
#' geom automatically. Titles, axis labels and other customizations can be added via
#' additional ggplot commands.
#'
#' @return A ggplot histogram plot object.
#'
#' @export
plot_gghist <- function(x) {
  dplyr::as_tibble(x) %>%
    tidyr::pivot_longer(
      cols = everything(),
      names_to = "cell",
      values_to = "expression"
    ) %>%
    ggplot2::ggplot(ggplot2::aes(x = expression, group = cell)) +
    ggplot2::geom_density() +
    ggplot2::coord_cartesian(ylim = c(0, 0.6), xlim = c(0, 3))
}
#' @export
matrix.prodSPARSE <- function(..., na.fill = 1) {
  matlist <- list(...)
  if (inherits(matlist[[1]], "list")) matlist <- matlist[[1]]
  if (is.null(sparse)) sparse <- is(matlist[[1]], "sparseMatrix")
  p <- Reduce("+", lapply(matlist, function(x) 1 * !is.na(x)))
  i <- 1
  for (i in 1:length(matlist)) {
    if (any(is.na(matlist[[i]]))) {
      mx <- 1
      if (na.fill == "mean") {
        mx <- mean(matlist[[i]], na.rm = TRUE)
      }
      jj <- which(is.na(matlist[[i]]), arr.ind = TRUE)
      matlist[[i]][jj] <- mx
    }
  }
  xprod <- Reduce("*", matlist)
  if (na.fill != "mean") {
    xprod <- xprod**(1 / p)
  }
  jj <- which(p == 0, arr.ind = TRUE)
  if (sparse) {
    xprod[is.na(xprod)] <- 0
    xprod <- Matrix::Matrix(xprod, sparse = TRUE)
  }
  if (length(jj) > 0) {
    if (sparse) jj <- jj[which(jj[, 2] == 1 | jj[, 1] == 1), ]
    xprod[jj] <- NA
  }
  return(xprod)
}
#' Perform differential expression analysis using limma
#'
#' Performs differential expression analysis on a gene expression matrix using limma.
#' Handles single sample case by duplicating sample. Auto-detects reference level.
#'
#' @param X A gene expression matrix with genes in rows and samples in columns
#' @param pheno A vector or factor of phenotype labels for the samples
#' @param fdr False discovery rate threshold for identifying differentially expressed genes
#' @param compute.means Logical indicating whether to compute group means
#' @param lfc Log fold change cutoff for differential expression
#' @param max.na Maximum proportion of missing values allowed in a gene
#' @param ref Character vector of reference phenotype labels
#' @param trend Logical indicating whether to fit a trend line
#' @param verbose Verbosity level for status messages
#'
#' @return A list containing the limma results, top table, and plot objects
#'
#' @export
gx.limma.SAVE <- function(X, pheno, fdr = 0.05, compute.means = TRUE, lfc = 0.20,
                          max.na = 0.20, ref = c(
                            "ctrl", "ctr", "control", "dmso", "nt", "0", "0h", "0hr",
                            "non", "no", "not", "neg", "negative", "ref", "veh", "vehicle",
                            "wt", "wildtype", "untreated", "normal", "false", "healthy"
                          ),
                          trend = FALSE, verbose = 1) {
  if (sum(duplicated(rownames(X))) > 0) {
    cat("WARNING:: matrix has duplicated rownames\n")
  }
  ## detect single sample case
  is.single <- (max(table(pheno)) == 1)
  if (is.single) {
    cat("WARNING:: no replicates, no stats...\n")
    X <- cbind(X, X)
    pheno <- c(pheno, pheno)
  }

  ## filter probes and samples
  ii <- which(rowMeans(is.na(X)) <= max.na)
  jj <- which(!is.na(pheno))
  if (verbose > 0) cat(sum(is.na(pheno) > 0), "with missing phenotype\n")
  X0 <- X[ii, jj]
  pheno0 <- as.character(pheno[jj])
  X0 <- X0[!(rownames(X0) %in% c(NA, "", "NA")), ]
  if (verbose > 0) {
    cat("analyzing", ncol(X0), "samples\n")
    cat("testing", nrow(X0), "features\n")
  }

  ## auto-detect reference
  pheno.ref <- c()
  ref.detected <- FALSE
  ref <- toupper(ref)

  is.ref <- (toupper(pheno0) %in% toupper(ref))
  ref.detected <- (sum(is.ref) > 0 && sum(!is.ref) > 0)

  if (ref.detected) {
    pheno.ref <- unique(pheno0[which(toupper(pheno0) %in% toupper(ref))])
    if (verbose > 0) cat("setting reference to y=", pheno.ref, "\n")
    bb <- c(pheno.ref, sort(setdiff(unique(pheno0), pheno.ref)))
  } else {
    if (verbose > 0) cat("WARNING: could not auto-detect reference\n")
    bb <- as.character(sort(unique(pheno0)))
    if (verbose > 0) cat("setting reference to first class", bb[1], "\n")
  }
  if (length(bb) != 2) {
    stop("gx.limma::fatal error:only two class comparisons")
    return
  }
  design <- cbind(1, pheno0 == bb[2])
  colnames(design) <- c("WT", "2vs1")
  d1 <- colnames(design)[1]
  d2 <- colnames(design)[2]
  fit <- limma::lmFit(X0, design)
  fit <- limma::eBayes(fit, trend = trend)
  top <- limma::topTable(fit, coef = d2, number = nrow(X0))
  if ("ID" %in% colnames(top)) {
    rownames(top) <- top$ID
    top$ID <- NULL
  }
  top <- top[rownames(X0), ]
  Matrix::head(top)

  ## only significant
  top <- top[which(top$adj.P.Val <= fdr & abs(top$logFC) >= lfc), ]
  if (verbose > 0) cat("found", nrow(top), "significant at fdr=", fdr, "and minimal FC=", lfc, "\n")

  if (compute.means && nrow(top) > 0) {
    avg <- t(apply(
      X0[rownames(top), ], 1,
      function(x) tapply(x, pheno0, mean, na.rm = TRUE)
    ))
    avg <- avg[, as.character(bb), drop = FALSE]
    colnames(avg) <- paste0("AveExpr.", colnames(avg))
    top <- cbind(top, avg)
  }
  top$B <- NULL

  if (is.single) {
    top$P.Value <- NA
    top$adj.P.Val <- NA
    top$t <- NA
  }

  ## reorder on fold change
  top <- top[order(abs(top$logFC), decreasing = TRUE), ]


  ## unlist???

  return(top)
}
#' @export
mat.downsample <- function(mat, m, n = -1, FUN = c("mean", "max"),
                           clust = FALSE) {
  downsample.rows <- function(x, n = 1000, FUN, clust = FALSE) {
    if (NCOL(x) <= n) {
      return(x)
    }
    rowSymmMax <- function(a) {
      mx <- cbind(
        apply(a, 1, max, na.rm = TRUE),
        -apply(a, 1, min, na.rm = TRUE)
      )
      mx[is.infinite(mx)] <- 0 ## OK???
      apply(mx, 1, max) * c(1, -1)[max.col(mx)]
    }
    if (clust) {
      x <- x[, nclust(x, distance = "euclidean", link = "complete")$order]
    }
    nn <- ncol(x)
    nd <- round(nn / n)
    jj <- seq(1, nn, nd) ## start index of block
    col.mx <- colMeans(x, na.rm = TRUE)
    col.sd <- apply(x, 2, sd, na.rm = TRUE)
    col.mx0 <- sapply(jj, function(i) mean(col.mx[i:min(nn, i + nd - 1)], na.rm = TRUE))
    #
    dx <- scale(x, scale = FALSE) ## just mean center
    ii <- lapply(jj, function(i) {
      i:min(nn, i + nd - 1)
    })
    sdx <- apply(dx, 2, sd, na.rm = TRUE)
    ii <- lapply(ii, function(a) a[order(-sdx[a])]) ## order on SD
    if (FUN == "max") {
      sx <- sapply(ii, function(i) rowSymmMax(dx[, i, drop = FALSE]))
    } else if (FUN == "mean") {
      sx <- sapply(ii, function(i) rowMeans(dx[, i, drop = FALSE], na.rm = TRUE))
    } else {
      stop("error. unknown function")
    }
    if (inherits(sx, "numeric")) sx <- matrix(sx, nrow = 1)
    sx <- t(t(sx) + col.mx0) ## add back pooled col.mean

    ## set rownames
    aa <- sapply(ii, function(i) paste(colnames(x)[i], collapse = " "))
    colnames(sx) <- aa
    rownames(sx) <- rownames(x)
    sx
  }
  if (NCOL(mat) == 0 && n > 0) {
    return(mat)
  }
  #
  FUN <- FUN[1]
  if (NCOL(mat) == 1 && !inherits(mat, "matrix")) {
    dx <- matrix(mat, ncol = 1)
    rownames(dx) <- names(mat)
    n <- -1
  } else {
    dx <- mat
  }
  if (nrow(dx) == 0 && m > 0) {
    return(mat)
  }
  if (n < ncol(dx) && n > 0) dx <- downsample.rows(dx, n, FUN = FUN, clust = clust)
  if (m < nrow(dx) && m > 0) dx <- t(downsample.rows(t(dx), m, FUN = FUN, clust = clust))
  if (NCOL(mat) == 1) dx <- dx[, 1]
  return(dx)
}
#' Calculate mean F statistics from limma differential expression analysis
#'
#' @title Calculate mean F statistics
#'
#' @description Calculates the mean F statistic across all genes from running limma differential expression analysis.
#' Useful for comparing overall separation between groups.
#'
#' @param X A gene expression matrix with genes in rows and samples in columns.
#' @param pheno A phenotype factor or vector indicating the phenotype of each sample.
#'
#' @details Runs limma differential expression between the phenotype groups for each gene.
#' Extracts the F statistics from the limma results and calculates the mean F statistic across all genes.
#' Higher mean F indicates better separation between phenotype groups.
#'
#' @return The mean F statistic across all genes.
#'
#' @export
gx.meanFstats <- function(X, pheno) {
  getF <- function(x, y) {
    ii <- which(!is.na(y))
    y1 <- y[ii]
    if (inherits(y1, "factor")) y1 <- factor(as.character(y1))
    design <- model.matrix(~y1)
    fit <- limma::lmFit(x[, ii], design)
    fit <- limma::eBayes(fit, trend = TRUE)
    top <- limma::topTableF(fit, number = nrow(x))
    mean(top$F)
  }
  fstat <- c()
  px <- playbase::tidy.dataframe(pheno) ## get variable types correct
  for (p in c("random", colnames(px))) {
    if (p == "random") {
      y <- sample(c("a", "b"), ncol(X), replace = TRUE)
    } else {
      y <- px[, p]
      p
    }
    fstat[p] <- getF(X, y)
  }
  fstat
}
#' @export
human2mouse.SLLOWWW <- function(x) {
  human <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse <- biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 <- biomaRt::getLDS(
    attributes = c("hgnc_symbol"),
    filters = "hgnc_symbol",
    values = x, mart = human,
    attributesL = c("mgi_symbol"),
    martL = mouse,
    uniqueRows = T
  )
  genesx <- unique(genesV2[, 2])


  return(genesx)
}
#' Plot GSEA Enrichment Plots
#'
#' Plot the GSEA enrichment plots from the GSEA analysis results.
#'
#' @param gsea The GSEA analysis results obtained from \code{\link{getGseaOutput}}.
#' @param gsets Optional character vector specifying the gene sets to plot.
#'              If NULL (default), all available gene sets will be plotted.
#' @param ncol The number of columns in the grid arrangement of plots.
#'             Defaults to 5.
#'
#' @export
#'
#' @return None
gseaPlotEnplots <- function(gsea, gsets = NULL, ncol = 5) {
  enplots <- gsea$enplots
  kk <- 1:length(enplots)
  if (!is.null(gsets)) kk <- match(gsets, names(enplots))
  imgs <- lapply(enplots[kk], function(p) {
    grid::rasterGrob(as.raster(png::readPNG(p)),
      interpolate = TRUE
    )
  })
  gridExtra::grid.arrange(grobs = imgs, ncol = ncol)
}
#' @export
pgx.createSeurateFigures <- function(obj) {
  caption1 <- paste("Project:", obj@project.name, "   Date:", Sys.Date())
  caption1
  fig <- list()

  ## ----------------------------------------------------------------------
  ## QC
  ## ----------------------------------------------------------------------
  Seurat::Idents(obj) <- "library_id" ## important!!
  Seurat::Idents(obj) <- "phenotype" ## important!!
  Seurat::Idents(obj) <- "batch" ## important!!
  Seurat::Idents(obj) <- "orig.indent" ## important!!
  vplot <- Seurat::VlnPlot(
    obj,
    features = c("nCount_RNA", "nFeature_RNA", "percent.mito", "percent.ribo"),
    pt.size = 0.15, ncol = 4
  )
  vplot <- vplot & ggplot2::xlab(NULL)

  q1 <- Seurat::FeaturePlot(obj, features = "nCount_RNA")
  q2 <- Seurat::FeaturePlot(obj, features = "nFeature_RNA")
  q3 <- Seurat::FeaturePlot(obj, features = "percent.mito")
  q4 <- Seurat::FeaturePlot(obj, features = "percent.ribo")

  plot1 <- Seurat::FeatureScatter(obj,
    feature1 = "nCount_RNA", feature2 = "nFeature_RNA",
    pt.size = 0.3
  )
  plot2 <- Seurat::FeatureScatter(obj,
    feature1 = "nCount_RNA", feature2 = "percent.mito",
    pt.size = 0.3
  )
  plot3 <- Seurat::FeatureScatter(obj,
    feature1 = "nCount_RNA", feature2 = "percent.ribo",
    pt.size = 0.3
  )

  theme0 <- ggplot2::theme(
    plot.title = ggplot2::element_text(size = 12, face = "bold"),
    legend.title = ggplot2::element_text(size = 11),
    legend.text = ggplot2::element_text(size = 9),
    legend.key.size = grid::unit(0.55, "lines"),
    axis.text.x = ggplot2::element_text(size = 10),
    axis.text.y = ggplot2::element_text(size = 10),
    axis.title.x = ggplot2::element_text(size = 11),
    axis.title.y = ggplot2::element_text(size = 11)
  )


  qq <- (q1 | q2 | q3 | q4)
  qq <- qq & ggplot2::guides(color = ggplot2::guide_colourbar(barwidth = 0.5, barheight = 3))

  pp <- (plot1 | plot2 | plot3) & ggplot2::ggtitle(NULL)
  pp <- pp & ggplot2::theme(legend.key.size = grid::unit(0.45, "lines"))
  ## pp

  v1 <- cowplot::plot_grid(vplot & theme0, ncol = 1)
  v2 <- cowplot::plot_grid(qq & theme0, ncol = 1)
  v3 <- cowplot::plot_grid(pp & theme0, ncol = 1)

  fig1 <- (v1 / v2 / v3)
  fig1 <- fig1 +
    patchwork::plot_annotation(
      title = "Seurat QC plots",
      subtitle = "These plots show the QC of your experiment.",
      caption = caption1,
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, face = "bold")),
      tag_levels = "a"
    )
  ## fig1
  fig[["QC"]] <- fig1

  ## ----------------------------------------------------------------------
  ## Variable features & PCA
  ## ----------------------------------------------------------------------

  vg <- obj[["RNA"]]@var.features
  hvg <- Matrix::head(Seurat::VariableFeatures(obj), 20)
  v1 <- Seurat::VariableFeaturePlot(obj) %>% Seurat::LabelPoints(points = hvg, repel = TRUE)
  v1 <- v1 & ggplot2::theme(legend.position = c(0.02, 0.95))


  r1 <- Seurat::ElbowPlot(obj) + ggplot2::ggtitle("PCA elbow plot")
  r2 <- base2grob::base2grob(~ {
    par(mfrow = c(1, 1), mar = c(0, 0, 0, 0), oma = c(0, 0, 1, 0))
    Seurat::DimHeatmap(obj,
      dims = 1:9, ncol = 3, cells = 500,
      nfeatures = 16, balanced = TRUE
    )
  })

  v1 <- v1 + ggplot2::ggtitle("Variable features") + theme0 +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 1, 0.5, 0, "cm"))

  r1 <- r1 + ggplot2::ggtitle("PCA elbow plot") + theme0 +
    ggplot2::theme(plot.margin = ggplot2::margin(0, 0.5, 0.5, 1, "cm"))

  fig3 <- (v1 | r1) / r2 + patchwork::plot_layout(design = "A\nB\nB") &
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0))

  fig3 <- fig3 +
    patchwork::plot_annotation(
      title = "Seurat PCA plots",
      subtitle = "These plots show the PCA of your experiment.",
      caption = caption1,
      theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, face = "bold")),
      tag_levels = "a"
    )

  ## fig3
  fig[["PCA"]] <- fig3

  ## ----------------------------------------------------------------------
  ## Cluster markers
  ## ----------------------------------------------------------------------
  Seurat::Idents(obj) <- "seurat_clusters"

  markers <- Seurat::FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)

  markers$score <- with(markers, -log(p_val) * avg_logFC * pct.1 * (1 - pct.2))


  top1 <- markers %>%
    plotly::group_by(cluster) %>%
    dplyr::top_n(1, -p_val)
  top1 <- markers %>%
    plotly::group_by(cluster) %>%
    dplyr::top_n(1, score)
  top1$gene
  vplot2 <- Seurat::VlnPlot(obj, features = top1$gene, pt.size = 0.15, ncol = 5)
  vplot2 <- vplot2 & ggplot2::xlab(NULL)
  vplot2 <- vplot2 & ggplot2::ylab("expression")
  vplot2 <- vplot2 & ggplot2::theme(plot.margin = ggplot2::margin(1, 3, 1, 3, "mm"))
  vplot2 <- vplot2 & ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0))

  ntop <- floor(65 / length(unique(obj$seurat_clusters)))

  top <- markers %>%
    plotly::group_by(cluster) %>%
    dplyr::top_n(ntop, score)
  top
  h1 <- Seurat::DoHeatmap(obj, features = top$gene, angle = 0, size = 4) + Seurat::NoLegend()

  vplot2 <- cowplot::plot_grid(vplot2 & theme0, ncol = 1)
  h1x <- cowplot::plot_grid(h1, ncol = 1)
  fig4 <- (vplot2) / h1x + patchwork::plot_layout(design = "A\nB\nB")

  fig4 <- fig4 + patchwork::plot_annotation(
    title = "Seurat cluster markers",
    subtitle = "These plots show the clusters of your experiment.",
    caption = caption1,
    theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, face = "bold")),
    tag_levels = "a"
  )
  ## fig4

  fig[["cluster.markers"]] <- fig4

  ## ----------------------------------------------------------------------
  ## Celltype assignment
  ## ----------------------------------------------------------------------

  CANONICAL.MARKERS <- list(
    "B_cells" = c("Ms4a1", "Cd79a", "Cd79b", "Fcmr", "Ebr1"),
    "CD4 T cells" = c("Il17r", "Ccr7", "Cd3e", "Cd3d", "Cd3g", "Cd4"),
    "CD8 T cells" = c("Cd8a", "Cd8b"),
    "NK_cells" = c("Gnly", "Nkg7", "Gzma", "Klrb1c", "Klrk1", "Klra4"),
    "Dendritic_cells" = c("Fcer1a", "Cst3", "Siglech", "Fscn1", "Ccl22"),
    "Macrophages" = c("C1qa", "C1qb", "C1qc", "Lyz2"),
    "Monocytes" = c("Cd14", "Lyz", "Ly6c2", "Fcgr3a", "Ms4a7"),
    "Platelet" = "Ppbp"
  )

  ct1 <- pgx.inferCellType(obj[["RNA"]]@counts, add.unknown = FALSE, low.th = 0.01)

  tapply(ct1, obj$seurat_clusters, function(x) table(x))
  ct1x <- tapply(ct1, obj$seurat_clusters, function(x) names(which.max(table(x))))
  ct1x

  obj$cell.type <- ct1x[as.character(obj$seurat_clusters)]

  d1 <- Seurat::DimPlot(obj, group.by = "seurat_clusters", label = TRUE) +
    ggplot2::ggtitle("Seurat clusters") + Seurat::NoLegend()
  d4 <- Seurat::DimPlot(obj, group.by = "cell.type", label = TRUE) +
    ggplot2::ggtitle("Cell type") + Seurat::NoLegend()

  sel <- unlist(tapply(1:ncol(obj), obj$cell.type, head, 300))
  markers$cell.type <- ct1x[as.character(markers$cluster)]

  ntop <- floor(45 / length(unique(ct1x)))
  top <- markers %>%
    plotly::group_by(cell.type) %>%
    dplyr::top_n(ntop, score)
  h2 <- Seurat::DoHeatmap(obj[, sel],
    features = top$gene, group.by = "cell.type",
    hjust = 0.5, angle = 0, size = 4
  ) + Seurat::NoLegend()
  h2 <- h2 + ggplot2::theme(plot.margin = ggplot2::margin(0, 0, 0, 0, "mm"))


  fig5 <- h2 / (d1 + theme0 | d4 + theme0) + patchwork::plot_layout(design = "A\nA\nB")


  fig5 <- fig5 + patchwork::plot_annotation(
    title = "Seurat cell type identification",
    subtitle = "Assignment of cell type identity to clusters.",
    caption = caption1,
    theme = ggplot2::theme(
      plot.title = ggplot2::element_text(size = 16, face = "bold")
    ),
    tag_levels = "a"
  )
  ## fig5
  fig[["cell.type"]] <- fig5

  ## ----------------------------------------------------------------------
  ## Phenotype dimension plots
  ## ----------------------------------------------------------------------

  not.ph <- grep("ident|nCount|nFeat|mito|[.]mt|ribo|_snn|^percent", colnames(obj@meta.data), value = TRUE)
  ph <- setdiff(colnames(obj@meta.data), not.ph)
  ph
  dd <- list()
  for (p in ph) {
    dd[[p]] <- Seurat::DimPlot(obj, group.by = p, label = TRUE) +
      ggplot2::ggtitle(p) + Seurat::NoLegend() + theme0
  }


  fig2 <- patchwork::wrap_plots(dd)

  fig2 <- fig2 + patchwork::plot_annotation(
    title = "Seurat Phenotype plots",
    subtitle = "Distribution of phenotypes on the t-SNE/UMAP.",
    caption = caption1,
    theme = ggplot2::theme(plot.title = ggplot2::element_text(size = 16, face = "bold")),
    tag_levels = "a"
  )
  ## fig2

  fig[["phenotypes"]] <- fig2

  ## ----------------------------------------------------------------------
  ## Return all figures
  ## ----------------------------------------------------------------------
  return(fig)
}
#' @export
silac.volcanoPlot <- function(obj, group1, group2, psig = 0.05, lfc = 0.2) {
  if (inherits(group1[1], "character") && group1[1] %in% names(obj$groups)) {
    group1 <- obj$groups[[group1]]
  }
  if (inherits(group2[1], "character") && group2[1] %in% names(obj$groups)) {
    group2 <- obj$groups[[group2]]
  }

  ## Test differences naive memory
  P <- as.matrix(obj$LFQ.ratio)
  group1 <- intersect(group1, colnames(P))
  group2 <- intersect(group2, colnames(P))
  genes <- as.character(obj$proteins$Gene.names)
  pdata <- volcano(P[, group1], P[, group2], genes)
  jj <- which(!is.na(pdata$P) & !is.nan(pdata$P))
  sig <- which(pdata$P < psig & abs(pdata$EFFECTSIZE) > lfc)
  highlight <- pdata$Gene[sig]

  volcanoly(pdata[jj, ],
    snp = "Gene", highlight = highlight,
    effect_size_line = c(-1, 1) * lfc,
    genomewideline = -log10(0.05)
  )
  invisible(pdata)
}
#' @title Surrogate variable analysis batch correction
#'
#' @description
#' Performs batch correction using surrogate variable analysis (SVA).
#' Estimates surrogate variables and removes them from the expression data.
#'
#' @param X Gene expression matrix, genes in rows and samples in columns.
#' @param batch Factor or matrix indicating batch groups for each sample.
#' @param model.par Vector of column names in \code{pheno} to use as model covariates.
#' @param pheno Dataframe containing sample metadata/covariates. Must match colnames of \code{X}.
#' @param n.sv Number of surrogate variables to estimate. If NULL uses a data-driven approach.
#'
#' @return Batch corrected gene expression matrix.
#'
#' @details
#' This function performs batch correction using the sva R package.
#'
#' It constructs a design matrix containing the specified batch groups and covariates.
#' Surrogate variables are estimated to capture remaining batch effects.
#'
#' The top surrogate variables are included in the design matrix and removed from the expression data.
#' This accounts for batch effects in a completely unsupervised way.
#'
#' @export
pgx.svaCorrect <- function(X, pheno, nmax = -1) {
  ##
  ## IK: not sure about this SVA correction stuff...
  if (NCOL(pheno) == 1) {
    pheno <- data.frame(pheno = pheno)
  }
  X <- as.matrix(X)

  ## add some random... sometimes necessary

  X <- X + 1e-8 ## ??

  ## setup model matrix
  mod1 <- c()
  pheno1 <- pheno
  colnames(pheno1) <- paste0(colnames(pheno1), "_IS_")
  for (v in colnames(pheno1)) {
    expr <- paste0("model.matrix(~", v, ",data=pheno1)")
    m1 <- eval(parse(text = expr))[, -1, drop = FALSE]
    mod1 <- cbind(mod1, m1)
  }
  colnames(mod1) <- sub("_IS_", "=", colnames(mod1))

  mod1x <- cbind(1, mod1)
  mod0x <- mod1x[, 1, drop = FALSE]


  message("Estimating number of surrogate variables...")

  ## fast method using SmartSVA

  pp <- paste0(colnames(pheno), collapse = "+")
  lm.expr <- paste0("lm(t(X) ~ ", pp, ", data=pheno)")
  X.r <- t(stats::resid(eval(parse(text = lm.expr))))
  n.sv <- isva::EstDimRMT(X.r, FALSE)$dim + 1
  suppressWarnings(min1 <- min(sapply(apply(pheno, 2, table), min)))
  n.sv <- min(n.sv, min1)

  message("Calculating SVA...")
  vX <- X
  if (nmax > 0) {
    vX <- Matrix::head(X[order(-apply(X, 1, sd)), ], nmax)
  }
  sv <- sva::sva(vX, mod1x, mod0x, n.sv = n.sv)$sv


  message("Perform batch correction...")

  cX <- limma::removeBatchEffect(X, covariates = sv, design = mod1x)


  ## recenter on old feature means
  cX <- cX - Matrix::rowMeans(cX, na.rm = TRUE) + Matrix::rowMeans(X, na.rm = TRUE)

  return(cX)
}
#' @export
silac.readDataFile <- function(datafile, remove.outliers = TRUE) {
  cat(">>> reading datafile from", datafile, "\n")
  df <- read.delim(datafile)
  rownames(df) <- paste0("tag", 1:nrow(df)) ## give proper rownames

  ## ------------ Filter contaminants (important)
  keep <- (df$Reverse != "+" & df$Only.identified.by.site != "+" & df$Potential.contaminant != "+")
  df1 <- df[keep, ]

  ## use LFQ values that are normalized across the samples

  LFQ.L <- as.matrix(df1[, grep("LFQ.intensity.L.", colnames(df1))])
  LFQ.H <- as.matrix(df1[, grep("LFQ.intensity.H.", colnames(df1))])

  ## --------- rename samples

  LFQ.names <- gsub(".*MB[0-9]*_|.*MB_RG_|\"", "", colnames(LFQ.L))

  ## Use new sample nomenclature
  names.new.df <- read.delim("samples_renamed.txt", stringsAsFactors = FALSE)

  sum(duplicated(names.new.df$Renamed))

  ## match names
  org.names <- gsub(".*MB[0-9]*_|.*MB_RG_|\"", "", names.new.df$Original)
  org.names <- sub("Rest", "Resting", org.names)
  names.new.df2 <- names.new.df[match(LFQ.names, org.names), ]

  ## clean names.... (OMG)
  clean.names <- names.new.df2$Renamed
  clean.names <- gsub("_SILAC", "-SILAC", clean.names)
  clean.names <- gsub("_Naive", "-Naive", clean.names)
  clean.names <- gsub("_Rest", "-Rest", clean.names)
  clean.names <- gsub("_Tcell", "-Tcell", clean.names)
  clean.names <- gsub("3-Methyladenosine", "3Methyladenosine", clean.names)

  clean.names <- gsub("^DorE", "DoE", clean.names)
  clean.names <- gsub("_", "", clean.names)

  ## write.csv( cbind(LFQ.names=LFQ.names,


  ##                  clean.names=clean.names ), file="checknames.csv")

  ## check duplicated
  ndup <- sum(duplicated(clean.names))
  ndup
  if (ndup > 0) {
    stop("FATAL. Duplicated clean names.")
  }

  ## update names
  colnames(LFQ.H) <- clean.names
  colnames(LFQ.L) <- clean.names

  ## ------------ annotation table (keep separate) ---------
  ## Include proteinID, gene names and MW in dataframe A
  proteins <- df1[, c("Protein.IDs", "Gene.names", "Mol..weight..kDa.")]

  ## ------------ filter samples  ---------------------------------------
  ## Remove wrong sample and bad donors
  if (remove.outliers) {
    bad.samples <- grep("Do109|Do4221-Tcell-Mem", colnames(LFQ.H), value = TRUE)
    keep.samples <- setdiff(colnames(LFQ.H), bad.samples)
    LFQ.L <- LFQ.L[, keep.samples]
    LFQ.H <- LFQ.H[, keep.samples]
  }

  ## ------------ filter probes ---------------------------------------
  ## Delete Factors from the list with severe miss identifications in the control group
  LFQ.ratio <- LFQ.H / (LFQ.H + LFQ.L + 1e-8) ## only temporary, later we calculate again
  Ctrls <- LFQ.ratio[, grep("Treat=NO-SILAC=0h", colnames(LFQ.ratio))]

  keep <- (rowSums(Ctrls) < 1.9)

  LFQ.L <- LFQ.L[keep, ]
  LFQ.H <- LFQ.H[keep, ]
  proteins <- proteins[keep, ]

  ## ------------ collapse by gene??
  gene1 <- as.character(proteins$Gene.names)

  sum(duplicated(gene1))
  LFQ.L <- apply(LFQ.L, 2, function(x) tapply(x, gene1, sum, na.rm = TRUE))
  LFQ.H <- apply(LFQ.H, 2, function(x) tapply(x, gene1, sum, na.rm = TRUE))
  jj <- match(rownames(LFQ.L), gene1)
  proteins <- proteins[jj, ]
  rownames(proteins) <- rownames(LFQ.L)

  ## ------------ create sample annotation

  samples <- sapply(colnames(LFQ.H), strsplit, split = "-")
  samples <- data.frame(do.call(rbind, samples))
  colnames(samples) <- c("donor", "cell.type", "subtype", "state", "treatment", "SILAC")
  rownames(samples) <- colnames(LFQ.H) ## should be already

  samples$treatment <- gsub("Treat=", "", samples$treatment)
  samples$SILAC <- gsub("SILAC=", "", samples$SILAC)

  ## add mass in picograms (PLEASE CHECK!!!)
  samples$mass.pg <- 25
  samples$mass.pg[which(samples$state %in% c("Act23h", "Act48h"))] <- 75

  ## ------------ define groups
  groups <- list()
  nav <- (samples$subtype == "Naive" & samples$state == "Rest" & samples$treatment == "NO")
  groups[["NaiveRest_0h"]] <- which(nav & samples$SILAC == "0h")
  groups[["NaiveRest_6h"]] <- which(nav & samples$SILAC == "6h")
  groups[["NaiveRest_12h"]] <- which(nav & samples$SILAC == "12h")
  groups[["NaiveRest_24h"]] <- which(nav & samples$SILAC == "24h")
  groups[["NaiveRest_48h"]] <- which(nav & samples$SILAC == "48h")

  mem <- (samples$subtype == "Mem" & samples$state == "Rest" & samples$treatment == "NO")
  groups[["MemRest_0h"]] <- which(mem & samples$SILAC == "0h")
  groups[["MemRest_6h"]] <- which(mem & samples$SILAC == "6h")
  groups[["MemRest_12h"]] <- which(mem & samples$SILAC == "12h")
  groups[["MemRest_24h"]] <- which(mem & samples$SILAC == "24h")
  groups[["MemRest_48h"]] <- which(mem & samples$SILAC == "48h")

  groups[["NaiveRest"]] <- which(nav)
  groups[["MemRest"]] <- which(mem)

  ## special groups
  groups[["NaiveRest_24h_Baf"]] <- which(samples$subtype == "Naive" & samples$state == "Rest" &
    samples$treatment == "Bafilomycin24h")
  groups[["NaiveRestWashout"]] <- which(samples$subtype == "Naive" & samples$state == "Rest" &
    samples$treatment == "CHX24hwashout24h")
  groups[["NaiveRestCHX"]] <- which(samples$subtype == "Naive" & samples$state == "Rest" &
    samples$treatment == "CHX24h")
  groups[["MemoryRest_24h_IL2"]] <- which(samples$subtype == "Mem" & samples$state == "Rest" &
    samples$treatment == "IL2")
  groups[["MemRestWashout"]] <- which(samples$subtype == "Mem" & samples$state == "Rest" &
    samples$treatment == "CHX24hwashout24h")
  groups[["MemRestCHX"]] <- which(samples$subtype == "Mem" & samples$state == "Rest" &
    samples$treatment == "CHX24h")

  ## IMPORTANT:: immediately convert indices to names. much safer!!!
  groups <- lapply(groups, function(ii) rownames(samples)[ii])
  sapply(groups, length) ## ALWAYS CHECK!!!


  ## ------------ add derived quantities: ratio
  LFQ.total <- (LFQ.L + LFQ.H) ## total expression
  LFQ.ratio <- LFQ.H / (LFQ.L + LFQ.H + 1e-8) ## notice adding small  number to avoid NaN

  molwt <- proteins$Mol..weight..kDa.
  copy.number <- silac.calcCopyNumber(data = LFQ.total, mol.weight = molwt, y = samples$mass.pg)

  ## prepare output object
  output <- list(
    samples = samples, groups = groups, proteins = proteins,
    LFQ.ratio = LFQ.ratio, LFQ.H = LFQ.H, LFQ.L = LFQ.L,
    copy.number = copy.number
  )

  return(output)
}
#' @export
hamming.match <- function(a, b, d = 0) {
  any(rowSums(tagged.hamming(a, b)) < d)
}
#' @export
length_similarityEXACT <- function(x, p = 1) {
  D <- 1.0 * abs(outer(x, x, FUN = "-"))
  L <- 1.0 * abs(outer(x, x, FUN = "pmax"))

  R <- exp(-p * (D / L))
  return(R)
}
#' Parse GEO series matrix file
#'
#' @param SERIES_FILE Character string. Path to GEO series matrix file.
#' @param PLATFORM_FILE Character string. Path to GEO platform file.
#' @param GENE_COLUMN Character string. Name of column containing gene symbols. Default "GENE_SYMBOL".
#' @param EXPRESSION_OUTPUT_FILE Character string. Name of expression matrix output file. Default "expression.csv".
#' @param ANNOTATION_OUTPUT_FILE Character string. Name of annotation output file. Default "annotation.csv".
#' @param write.file Logical. Whether to write output files. Default TRUE.
#'
#' @return List containing expression matrix and annotation data frame if write.file is FALSE.
#'
#' @details Parses a GEO series matrix file and corresponding platform file to extract expression data and annotations.
#' The gene symbols, expression values, and annotations are extracted into R objects.
#' The expression matrix and annotation data frame can be optionally written to CSV files.
#'
#' @export
parse_geo_series_matrix <- function(SERIES_FILE,
                                    PLATFORM_FILE,
                                    GENE_COLUMN = "GENE_SYMBOL",
                                    EXPRESSION_OUTPUT_FILE = "expression.csv",
                                    ANNOTATION_OUTPUT_FILE = "annotation.csv",
                                    write.file = TRUE) {
  ##
  ## https://gist.github.com/SimonLarsen
  ## https://gist.github.com/SimonLarsen/66f27c188039129f3510669b992b2c99

  ## Read characteristics
  con <- file(SERIES_FILE, "r")
  characteristics <- c()
  while (TRUE) {
    line <- readLines(con, n = 1)
    if (length(line) == 0) {
      break
    } else if (startsWith(line, "!Sample_title")) {
      titles <- unlist(strsplit(line, "\t"))[-1]
      titles <- gsub("\\\"", "", titles)
    } else if (startsWith(line, "!Sample_characteristics")) {
      characteristics <- c(characteristics, line)
    } else if (startsWith(line, "!Sample_geo_accession")) {
      accession <- unlist(strsplit(line, "\t"))[-1]
      accession <- gsub("\\\"", "", accession)
    }
  }
  close(con)

  ## Parse characteristics
  anno <- data.frame(lapply(characteristics, function(x) {
    values <- unlist(strsplit(x, "\t"))[-1]
    values <- gsub("\\\"", "", values)
    parts <- strsplit(values, ": ")

    name <- parts[[1]][[1]]
    values <- sapply(parts, function(x) x[2])

    out <- list()
    out[[name]] <- values
    return(out)
  }))

  anno <- data.table(sample = accession, title = titles, anno)

  ## Read probe-level expression data
  D <- fread(SERIES_FILE, header = TRUE, skip = "\"ID_REF\"", fill = TRUE, na.strings = c("", "NA", "null"))
  D <- D[1:(nrow(D) - 1), ] # remove table end marker
  D <- as.matrix(D[, -1], rownames = D$ID_REF)


  ## Read platform data table
  ref <- read.table(PLATFORM_FILE, header = TRUE, sep = "\t", quote = "", comment.char = "#", fill = TRUE)[, c("ID", GENE_COLUMN)]
  colnames(ref) <- c("ID_REF", "gene")
  ref$gene <- as.character(ref$gene)
  ref <- subset(ref, !is.na(gene) & gene != "")

  gene_split <- strsplit(ref$gene, " /// ")
  ref <- data.frame(
    ID_REF = rep(ref$ID_REF, sapply(gene_split, length)),
    gene = unlist(gene_split)
  )

  ## Merge tables to map Gene genes ids

  gene <- ref$gene[match(rownames(D), as.character(ref$ID_REF))]

  ## Aggregate duplicate genes by median
  sum(duplicated(gene))

  ex <- do.call(rbind, tapply(1:nrow(D), gene, function(i) colMeans(D[i, , drop = FALSE])))

  ## Transpose matrix, add sample column and set gene names as column names


  ## drop empty or duplicated columns
  sel1 <- (apply(anno, 2, function(x) mean(x %in% c(NA, "NA", "null", ""))) < 1)
  sel2 <- !duplicated(t(anno))
  sel <- which(sel1 & sel2)

  anno1 <- anno[, sel, with = FALSE]

  ## Write results to separate expression and annotation files
  if (write.file) {
    cat("writing", EXPRESSION_OUTPUT_FILE, "\n")
    fwrite(data.table(gene = rownames(ex), ex), file = EXPRESSION_OUTPUT_FILE, sep = ",")

    cat("writing", ANNOTATION_OUTPUT_FILE, "\n")
    fwrite(anno1, file = ANNOTATION_OUTPUT_FILE, sep = ",")
  }

  res <- list()
  res$values <- ex
  res$anno <- anno1
  res
}
#' Estimate copy number from gene expression
#'
#' @param ngs An NGS object containing gene expression data
#' @param nsmooth Smoothing window size for copy number estimation. Default 40.
#'
#' @return A list with estimated copy number values, chromosome, position,
#' and ideogram image for each gene.
#'
#' @description Estimates copy number variation from gene expression data by smoothing
#' relative expression values within genomic windows.
#'
#' @details This function takes an NGS object containing normalized gene expression data.
#' It calculates the relative expression level of each gene compared to the mean expression.
#' These relative levels are then smoothed within sliding genomic windows of size \code{nsmooth}
#' genes to estimate regional copy number variation. The smoothed values are returned along
#' with the chromosome, genomic position, and an ideogram image for data visualization.
#'
#' @export
pgx.CNAfromExpression <- function(ngs, nsmooth = 40) {
  ## This estimates CNV by local smoothing of relative expression
  ## values.
  symbol <- as.vector(as.list(org.Hs.eg.db::org.Hs.egSYMBOL))
  chrloc <- as.list(org.Hs.eg.db::org.Hs.egCHRLOC)
  chr <- as.vector(sapply(chrloc, function(x) names(x)[1]))
  pos <- abs(as.integer(as.vector(sapply(chrloc, function(x) x[1]))))
  chr[sapply(chr, is.null)] <- NA
  chr <- as.character(unlist(chr))
  chr <- chr[match(ngs$genes$gene_name, symbol)]
  pos <- pos[match(ngs$genes$gene_name, symbol)]
  genes <- data.frame(chr = chr, pos = pos)
  rownames(genes) <- ngs$genes$gene_name

  sel <- which(!is.na(genes$chr) & !is.na(genes$pos))
  genes <- genes[sel, ]

  if (!is.null(ngs$counts)) {
    cna <- log2(100 + edgeR::cpm(ngs$counts)) ## moderated log2
  } else {
    cna <- ngs$X
  }
  gg <- intersect(rownames(genes), rownames(cna))
  cna <- cna[gg, ]
  genes <- genes[gg, ]

  ## ---------------------------------------------------------------------
  ## order genes and matrix according genomic position
  ## ---------------------------------------------------------------------
  jj <- which(genes$chr %in% c(1:22, "X", "Y"))
  genes <- genes[jj, ]
  genes$chr <- factor(genes$chr, levels = c(1:22, "X", "Y"))
  jj <- order(genes$chr, genes$pos)
  genes <- genes[jj, ]
  cna <- cna[rownames(genes), ]
  cna0 <- cna

  ## ---------------------------------------------------------------------
  ## apply 'crude' moving average filter (THIS SHOULD BE IMPROVED!)
  ## ---------------------------------------------------------------------
  mavg <- function(x, n = nsmooth) {
    stats::filter(x, rep(1 / n, n), sides = 2, circular = TRUE)
  }
  cna <- t(scale(t(cna), center = FALSE)) ## z-score
  cna <- apply(cna, 2, mavg)
  cna <- cna - apply(cna, 1, median, na.rm = TRUE)
  rownames(cna) <- rownames(cna0)

  res <- list(cna = cna, chr = genes$chr, pos = genes$pos)
  return(res)
}
#' @export
mat2sp.DEPRECATED <- function(x) {
  idx <- which(x != 0, arr.ind = TRUE)
  dn <- list(rownames(x), colnames(x))
  spX <- Matrix::sparseMatrix(
    i = idx[, 1], j = idx[, 2], x = x[idx],
    dims = c(nrow(x), ncol(x)), dimnames = dn
  )
  spX
}
#' @export
pgx.ReclusterSignatureDatabase <- function(h5.file, reduce.sd = 1000, reduce.pca = 100) {
  h5exists <- function(h5.file, obj) {
    xobjs <- apply(rhdf5::h5ls(h5.file)[, 1:2], 1, paste, collapse = "/")
    obj %in% gsub("^/|^//", "", xobjs)
  }

  X <- rhdf5::h5read(h5.file, "data/matrix")
  rn <- rhdf5::h5read(h5.file, "data/rownames")
  cn <- rhdf5::h5read(h5.file, "data/colnames")
  rownames(X) <- rn
  colnames(X) <- cn
  X[which(X < -999999)] <- NA

  ## --------------------------------------------------
  ## Precalculate t-SNE/UMAP
  ## --------------------------------------------------

  if (!h5exists(h5.file, "clustering")) rhdf5::h5createGroup(h5.file, "clustering")

  pos <- pgx.clusterBigMatrix(
    abs(X), ## on absolute foldchange!!
    methods = c("pca", "tsne", "umap"),
    dims = c(2, 3),
    reduce.sd = reduce.sd,
    reduce.pca = min(reduce.pca, round(ncol(X) / 3))
  )

  rhdf5::h5write(pos[["pca2d"]], h5.file, "clustering/pca2d") ## can write list??
  rhdf5::h5write(pos[["pca3d"]], h5.file, "clustering/pca3d") ## can write list??
  rhdf5::h5write(pos[["tsne2d"]], h5.file, "clustering/tsne2d") ## can write list??
  rhdf5::h5write(pos[["tsne3d"]], h5.file, "clustering/tsne3d") ## can write list??
  rhdf5::h5write(pos[["umap2d"]], h5.file, "clustering/umap2d") ## can write list??
  rhdf5::h5write(pos[["umap3d"]], h5.file, "clustering/umap3d") ## can write list??
  rhdf5::h5closeAll()
}
#' @export
graph_from_knn <- function(pos, k = 10) {
  ## return: edge weight are distances
  if (is.null(rownames(pos))) stop("pos must have rownames")



  if (ncol(pos) > 3 || NCOL(pos) == 1) {
    stop("positions must be 2 or 3 columns\n")
  }
  ## use fast KNN package

  res <- FNN::get.knn(pos, k = k)
  idx <- res$nn.index
  xval <- res$nn.dist
  xval <- as.vector(unlist(xval))
  sp.idx <- do.call(rbind, lapply(1:nrow(idx), function(i) cbind(i, idx[i, ])))
  sp <- Matrix::sparseMatrix(i = sp.idx[, 1], j = sp.idx[, 2], x = xval, dims = c(nrow(pos), nrow(pos)))
  sp <- (sp + t(sp)) / 2
  rownames(sp) <- colnames(sp) <- rownames(pos)
  g <- igraph::graph_from_adjacency_matrix(sp, mode = "undirected", diag = FALSE, weighted = TRUE)
  g$layout <- pos
  return(g)
}
#' @title Import Salmon quantification results
#'
#' @description
#' Imports Salmon transcript-level quantification results into an NGS object.
#'
#' @param sampleTable A data frame with sample info and paths to Salmon sf files.
#' @param gencode Annotations for transcripts (e.g. from Ensembl).
#' @param txi Pre-computed Salmon tximport object. If NULL will be calculated.
#'
#' @details
#' This function takes a data frame \code{sampleTable} containing sample information and paths to Salmon sf files generated by Salmon quantification.
#' It runs \code{\link{tximport}} to import the Salmon results and aggregate to gene-level counts.
#'
#' Gene annotations are added from the provided \code{gencode} data frame based on Ensembl IDs.
#' The output is an \code{\link[edgeR]{DGEList}} object containing the imported counts and sample data ready for differential analysis.
#'
#' If a \code{txi} object is already provided, it will be used instead of re-importing the Salmon results.
#'
#' @return
#' A DGEList containing imported counts, sample data, and gene annotations.
#'
#' @export
ngs.rawSalmon <- function(sampleTable, gencode, txi = NULL) {
  if (is.null(txi)) {
    if (!("sf.file" %in% colnames(sampleTable))) stop("need sf.files in table")
    sf.files <- sampleTable$sf.files
    names(sf.files) <- rownames(sampleTable)
    txi <- ngs.tximportSalmon(sf.files)
  } else {
    if (!all(rownames(sampleTable) == colnames(txi$counts))) {
      stop("provided txi does not match sampleTable")
    }
  }

  keys <- rownames(txi$counts)
  geneTable <- ngs.getGeneAnnot(keys, keytype = "ENSEMBL", gencode)

  ## we like this..
  probegene <- paste0(rownames(geneTable), ":", geneTable$gene_name)
  Matrix::head(probegene)
  rownames(txi$counts) <- probegene
  rownames(txi$abundance) <- probegene
  rownames(geneTable) <- probegene
  all(rownames(txi$counts) == rownames(geneTable))

  raw <- edgeR::DGEList(round(txi$counts), group = NULL) ## we like integer counts...
  raw$group <- NULL
  raw$samples <- sampleTable
  raw$genes <- geneTable

  return(raw)
}
#' @title Build gene association graph from PubMed identifiers
#'
#' @description Builds an undirected gene association graph from PubMed
#' identifiers (PMIDs) mapped to gene symbols.
#'
#' @param P Sparse matrix mapping PMIDs to gene symbols
#'
#' @details This function takes as input a sparse matrix \code{P} that maps
#' PubMed identifiers (PMIDs) to gene symbols.
#' It filters \code{P} to only include PMIDs mapped to 2-10 genes.
#'
#' It then calculates a co-occurrence matrix \code{M} between genes by multiplying
#' \code{P} and its transpose. \code{M[i,j]} indicates the number of PMIDs shared
#' between gene i and gene j.
#'
#' An undirected graph is constructed from \code{M} using the igraph package,
#' with genes as nodes and edge weights proportional to their PMID co-occurrences.
#'
#' @return An igraph undirected graph object representing the gene association network.
#'
#' @export
pmid.buildGraph <- function(P) {
  P <- P[which(Matrix::rowSums(P) <= 10), ]
  P <- P[which(Matrix::rowSums(P) >= 2), ]
  P <- P[, which(Matrix::colSums(P) > 0)]
  P[1:10, 1:10]

  ## create graph from overlap
  M <- P[, ] %*% t(P[, ])
  diag(M) <- 0
  object.size(M)
  gr <- igraph::graph_from_adjacency_matrix(M,
    mode = "undirected",
    diag = FALSE, weighted = TRUE
  )
  igraph::V(gr)$name <- rownames(M)
  gr <- igraph::subgraph.edges(gr, which(igraph::E(gr)$weight > 0))

  P <- P[igraph::V(gr)$name, ]
  pmids <- parallel::mclapply(igraph::V(gr)$name, function(x) gsub("PMID:", "", strsplit(x, split = ",")[[1]]))
  nref <- sapply(pmids, length)
  vgenes <- parallel::mclapply(1:nrow(P), function(i) names(which(P[i, ] != 0)), mc.cores = NCORE())
  vgenes2 <- unlist(sapply(vgenes, paste, collapse = ","))
  igraph::V(gr)$size <- nref
  igraph::V(gr)$genes <- vgenes
  igraph::V(gr)$pmid <- pmids
  return(gr)
}
#' @title Correct gene symbols with September/March suffix
#'
#' @description
#' Corrects gene symbols that use "Sep" or "Mar" suffix by converting them to full month names.
#'
#' @param gg A character vector of gene symbols.
#'
#' @details
#' This function checks if any gene symbols in \code{gg} end with "-Sep", "-SEP", "-Mar", or "-MAR",
#' which represent September or March. If found, it converts these to full month names by replacing
#' "-Sep" with "-September" and "-Mar" with "-March".
#'
#' It first replaces any "-SEP" or "-MAR" with "-Sep" and "-Mar" respectively.
#' Then it uses \code{\link[plyr]{mapvalues}} to match any "-Sep" or "-Mar" suffixes to the
#' corresponding "-September" or "-March" and replaces them.
#'
#' @return
#' The character vector \code{gg} with any September/March suffixes converted to full month names.
#'
#' @export
correctMarchSeptemberGenes <- function(gg) {
  sep.from <- c(paste0("0", 1:9, "-Sep"), paste0(1:19, "-Sep"))
  sep.to <- c(paste0("SEPT", 1:9), paste0("SEPT", 1:19))
  mar.from <- c(paste0("0", 1:9, "-Mar"), paste0(1:19, "-Mar"))
  mar.to <- c(paste0("MARCH", 1:9), paste0("MARCH", 1:19))


  from <- c(sep.from, mar.from)
  to <- c(sep.to, mar.to)
  gg1 <- sub("[.-]Sep$|[.-]SEP$", "-Sep", gg)
  gg1 <- sub("[.-]Mar$|[.-]MAR$", "-Mar", gg1)
  jj <- which(from %in% gg1)
  gg2 <- gg1
  if (length(jj) > 0) {
    cat("Found ", length(jj), "Sept/Mar genes!\n")

    from[jj]
    gg2 <- plyr::mapvalues(gg1, from[jj], to[jj])
  }
  return(gg2)
}
#' @title Match gene names to feature identifiers
#'
#' @description Matches gene names in a character vector to the
#' corresponding feature identifiers in an NGS object.
#'
#' @param ngs An NGS object containing feature annotation data.
#' @param genes Character vector of gene names to match.
#'
#' @details This function takes a character vector of gene names and
#' finds the matching feature identifiers in the feature annotation data
#' contained in the \code{ngs$genes} data frame of an NGS object. It
#' matches the \code{gene_name} column in \code{ngs$genes}.
#'
#' Matching is case-insensitive. The function returns the feature
#' identifiers from the rownames of \code{ngs$genes} that match the
#' provided gene names.
#'
#' @return Character vector of matching feature identifiers
#' (typically Ensembl IDs).
#'
#' @examples
#' \dontrun{
#' # TODO
#' }
#'
#' @export
ngs.matchFeatures <- function(ngs, genes) {
  jj <- match(toupper(genes), toupper(ngs$genes$gene_name))
  rownames(ngs$genes)[jj]
}
#' Detect time variables in expression matrix
#'
#' @param Y Expression matrix with samples in columns
#'
#' @return A character vector of detected time variables
#'
#' @description
#' Detects column names in an expression matrix \code{Y} that represent time variables.
#'
#' @details
#' This function checks the column names of the expression matrix \code{Y} for common time variable names like
#' "time", "hour", "hr", or "day". It returns a character vector of any matched column names that represent time points.
#'
#' This is useful for identifying potential time variables that could be used for time course analyses.
#'
#' @export
pgx.detect_timevar <- function(Y) {
  has.time <- any(grepl("time|hour|hr|day", colnames(Y)))
  has.time
  ttvar <- grep("time|hour|hr|day", colnames(Y), value = TRUE)
  return(ttvar)
}
#' Collapse rows by common symbol
#'
#' @title Collapse rows by common symbol
#'
#' @description
#' Collapses rows of a matrix based on a common symbol, removing duplicates.
#'
#' @param X Numeric matrix with rows to collapse.
#' @param symbol Character vector of symbols corresponding to rows of X.
#'
#' @details
#' This function takes a matrix X and a corresponding symbol vector.
#' It orders the rows by standard deviation.
#' It then keeps the first unique occurrence of each symbol, removing duplicates.
#' The result is a matrix with rows collapsed based on the provided symbols.
#'
#' @return
#' A matrix with rows collapsed based on the provided symbols.
#'
#' @export
gx.collapse2symbol <- function(X, symbol) {
  j1 <- order(-apply(X, 1, sd))
  X <- X[j1, ]
  symbol <- symbol[j1]
  j2 <- which(!duplicated(symbol) & !(symbol %in% c("", " ", "NA", NA)))
  X <- X[j2, ]
  rownames(X) <- symbol[j2]
  return(X)
}
#' @export
length_similaritySPARSE <- function(x, r = 0.1) {
  X <- length_encode(x, r = r)
  S <- cosine_similarity(t(X))
  rownames(S) <- colnames(S) <- NULL
  return(S)
}
#' Detect organism from NGS object
#'
#' @title Detect organism from NGS object
#'
#' @description Detects the organism for an NGS object by
#' looking at the case of gene identifiers.
#'
#' @param ngs An NGS object.
#'
#' @details This function detects the organism for an NGS object
#' by calculating the ratio of lowercase characters in the
#' rownames of the counts matrix. It assumes human data will have a
#' high proportion of uppercase ENSEMBL identifiers, while mouse
#' data will have more lowercase.
#'
#' It calculates the ratio of lowercase characters in the first 100
#' characters of the rownames of ngs$counts.
#' If the ratio is greater than 0.5 it assigns 'mouse', otherwise 'human'.
#'
#' @return A character string indicating the detected organism ('human' or 'mouse').
#'
#' @examples
#' \dontrun{
#' # TODO
#' }
#' @export
ngs.detectOrganism <- function(ngs) {
  lowcase.ratio <- mean(grepl("[a-z]", substring(rownames(ngs$counts), 2, 100)))
  c("human", "mouse")[1 + 1 * (lowcase.ratio > 0.5)]
}
#' @export
downsample_graph.DEPRECATED <- function(g, idx = NULL, grouped = NULL, merge.op = "max") {
  if (is.null(g$layout)) stop("graph must have layout")

  if (is.null(idx)) {
    hc <- hclust_graph(g, k = NULL)
    idx <- hc[, ncol(hc)]
  }

  if (!is.null(grouped)) {
    grouped[is.na(grouped)] <- "_other"
    table(grouped)
    grouped.idx <- idx
    groups <- unique(grouped)
    i <- 1
    for (i in 1:length(groups)) {
      jj <- which(grouped == groups[i])
      grouped.idx[jj] <- paste0(i, ":", idx[jj])
    }
    idx <- grouped.idx
  } else {
    idx <- as.character(idx)
  }
  table(idx)

  ## resort indices
  idx.levels <- names(sort(-table(idx)))
  idx <- factor(idx, levels = idx.levels)
  table(idx)

  ## get adjacency matrix
  adjM <- g[, ]
  diag(adjM) <- 0
  adjM[1:4, 1:4]

  ## downsample by taking mean of merge blocks in adjM

  ngrp <- length(unique(idx))
  ngrp
  if (merge.op == "mean") {
    subM <- tapply(1:ncol(adjM), idx, function(ii) rowMeans(adjM[, ii, drop = FALSE], na.rm = TRUE))
    subM <- parallel::mclapply(subM, function(x) tapply(x, idx, mean, na.rm = TRUE))
  } else if (merge.op == "max") {
    subM <- tapply(1:ncol(adjM), idx, function(ii) apply(adjM[, ii, drop = FALSE], 1, max, na.rm = TRUE))
    subM <- parallel::mclapply(subM, function(x) tapply(x, idx, max, na.rm = TRUE))
  } else {
    stop("unknown merge.op", merge.up, "\n")
  }

  subM.idxnames <- names(subM)
  subM <- matrix(unlist(subM), ncol = ngrp, nrow = ngrp)
  sum(is.na(subM))
  rownames(subM) <- colnames(subM) <- subM.idxnames
  adjM[1:4, 1:4]
  dim(subM)
  subM[1:4, 1:4]

  ## diagonal??
  if (0) {
    for (i in 1:length(subM.idxnames)) {
      ii <- which(idx == subM.idxnames[i])
      subM[i, i] <- min(adjM[ii, ii], na.rm = TRUE)
    }
  }

  ## create downsampled graph
  subG <- igraph::graph_from_adjacency_matrix(
    subM,
    weighted = TRUE, diag = FALSE, mode = "undirected"
  )
  igraph::V(subG)$name <- rownames(subM)
  members <- tapply(igraph::V(g)$name, idx, list)
  names(members) <- rownames(subM)
  igraph::V(subG)$members <- members

  avg.pos <- apply(g$layout[igraph::V(g)$name, ], 2, function(x) tapply(x, idx, mean))
  rownames(avg.pos) <- rownames(subM)
  subG$layout <- avg.pos

  subG <- igraph::simplify(subG)
  names(idx) <- igraph::V(g)$name
  subG$membership <- idx
  return(subG)
}
#' @describeIn eset.getPhenoData Extracts the organism from the phenotype data of an ExpressionSet object
#' @export
eset.getOrganism <- function(eset) unique(as.character(eset.getPhenoData(eset, "organism_ch1")))

#' @describeIn eset.getPhenoData
#' @export
eset.getCH1 <- function(eset) {
  pdata <- pData(phenoData(eset))
  pdata <- pdata[, grep(":ch1$", colnames(pdata)), drop = FALSE]
  clin_info <- NULL
  has.clin <- "Clinical info:ch1" %in% colnames(pdata)
  has.clin
  if (has.clin) {
    clin_ch1 <- as.character(pdata[, "Clinical info:ch1"])
    clin.terms <- sub(":.*$", "", strsplit(clin_ch1[1], split = ";")[[1]])
    clin_info <- t(sapply(clin_ch1, function(s) (strsplit(s, split = ";")[[1]])))
    clin_info <- apply(clin_info, 2, function(s) sub(".*[:] ", "", s))
    clin_info <- data.frame(clin_info, stringsAsFactors = FALSE, check.names = FALSE)
    rownames(clin_info) <- NULL
    colnames(clin_info) <- clin.terms
    pdata <- cbind(
      pdata[, -which(colnames(pdata) == "Clinical info:ch1"), drop = FALSE],
      clin_info
    )
  }
  colnames(pdata) <- sub(":ch1", "", colnames(pdata))

  return(pdata)
}
#' @export
matrix.mean.SAVE <- function(matlist, weights = 1) {
  matx <- lapply(matlist, function(x) as.vector(x))
  matx <- as.matrix(do.call(rbind, matx))
  if (length(weights) == 1) weights <- rep(weights, length(matlist))
  weigths <- weights / sum(weights)
  matx <- matx * weights
  matx <- colMeans(matx * weights, na.rm = TRUE)
  matrix(matx, nrow = nrow(matlist[[1]]), ncol = ncol(matlist[[1]]))
}
#' @export
pgx.selectTCGAstudies <- function(cancertype, variables) {
  ## Scan the available TCGA studies for cancertype and clinical
  ## variables.
  ##
  ##
  ##

  mycgds <- cgdsr::CGDS("http://www.cbioportal.org/")
  all.studies <- sort(cgdsr::getCancerStudies(mycgds)[, 1])
  studies <- grep(cancertype, all.studies, value = TRUE)
  clin <- list()
  samples <- list()
  studies
  mystudy <- studies[1]

  for (mystudy in studies) {
    mystudy
    myprofiles <- cgdsr::getGeneticProfiles(mycgds, mystudy)[, 1]
    myprofiles

    ## mrna datatypes
    mrna.type <- "mrna"
    if (any(grepl("rna_seq_mrna$", myprofiles))) mrna.type <- "rna_seq_mrna"
    if (any(grepl("v2_mrna$", myprofiles))) mrna.type <- "rna_seq_v2_mrna"
    pr.mrna <- grep(paste0(mrna.type, "$"), myprofiles, value = TRUE)
    pr.mrna
    if (length(pr.mrna) == 0) next()

    all.cases <- cgdsr::getCaseLists(mycgds, mystudy)[, 1]
    all.cases


    caselist <- grep(paste0(mrna.type, "$"), all.cases, value = TRUE)
    caselist
    clin0 <- cgdsr::getClinicalData(mycgds, caselist)
    rownames(clin0) <- gsub("[.]", "-", rownames(clin0)) ## correct names...
    clin[[mystudy]] <- clin0
    samples[[mystudy]] <- rownames(clin0)
  }

  sel <- sapply(clin, function(v) any(grepl(variables, colnames(v))))
  sel
  sel.studies <- studies[sel]
  sel.clin <- clin[sel]

  res <- list(
    studies = sel.studies,
    clinicalData = sel.clin
  )
  return(res)
}
#' @export
sparse.corDEPRECATED <- function(a, b = NULL) {
  ## see also corSparse in package 'qlcMatrix'

  if (is.null(b)) b <- a
  a <- Matrix(a, sparse = TRUE)
  b <- Matrix(b, sparse = TRUE) ## make sure it is sparse coded
  n <- nrow(a)
  aMeans <- Matrix::colMeans(a, na.rm = TRUE)
  bMeans <- Matrix::colMeans(b, na.rm = TRUE)
  aSD <- ((Matrix::colSums(a**2, na.rm = TRUE) - n * aMeans**2) / (n - 1))**0.5
  bSD <- ((Matrix::colSums(b**2, na.rm = TRUE) - n * bMeans**2) / (n - 1))**0.5
  mm <- tcrossprod(aMeans, bMeans)
  ss <- tcrossprod(aSD, bSD)
  rr <- as.matrix(crossprod(a, b))
  covmat <- (rr - n * mm) / (n - 1)
  cormat <- covmat / ss
  return(cormat)
}
#' @export
pgx.getTCGA.multiomics.TOBEFINISHED <- function(studies, genes = NULL, batch.correct = TRUE,
                                                tcga.only = TRUE) {
  ## Better use curatedTCGA bioconductor package!!!!
  ##



  mycgds <- cgdsr::CGDS("http://www.cbioportal.org/")

  GENE <- "CIITA"
  GENE <- "NLRC5"

  ## Gather data from all cancers
  all.X <- list()
  mystudy <- studies[1]
  for (mystudy in studies) {
    mystudy

    myprofiles <- cgdsr::getGeneticProfiles(mycgds, mystudy)[, 1]
    myprofiles

    ## prioritize datatypes
    pr.mrna <- grep("rna_seq_v2_mrna$|rna_seq_mrna$", myprofiles, value = TRUE)[1]

    pr.cna <- grep("_log2CNA$|_linear_CNA$", myprofiles, value = TRUE)[1]
    pr.gistic <- grep("_gistic$", myprofiles, value = TRUE)[1]
    pr.me <- grep("_methylation_hm450|_methylation_hm27", myprofiles, value = TRUE)[1]
    pr.mut <- grep("_mutations", myprofiles, value = TRUE)[1]

    all.cases <- cgdsr::getCaseLists(mycgds, mystudy)[, 1]
    all.cases
    if (!any(grepl("complete$", all.cases))) next

    caselist <- grep("complete$", all.cases, value = TRUE)
    cna <- counts <- cna.gistic <- me <- mut <- gx <- NULL
    counts <- cgdsr::getProfileData(mycgds, genes, pr.mrna, caselist)
    cna <- cgdsr::getProfileData(mycgds, GENE, pr.cna, caselist)

    cna.gistic <- cgdsr::getProfileData(mycgds, GENE, pr.gistic, caselist)
    me <- cgdsr::getProfileData(mycgds, GENE, pr.me, caselist)
    mut <- cgdsr::getProfileData(mycgds, GENE, pr.mut, caselist)
    mut <- 1 * !is.na(mut)
    gx <- log2(10 + as.matrix(counts))
    cna[is.na(cna)] <- NA
    if (grepl("linear", pr.cna)) cna <- log2(0.01 + 2 + cna) ## assume diploid


    if (!is.null(cna)) colnames(cna) <- paste0("CN:", colnames(cna))
    if (!is.null(cna.gistic)) colnames(cna.gistic) <- paste0("CNA:", colnames(cna.gistic))
    if (!is.null(me)) colnames(me) <- paste0("ME:", colnames(me))
    if (!is.null(mut)) colnames(mut) <- paste0("MT:", colnames(mut))



    xx <- list(gx, cna.gistic, me, mut)
    xx <- xx[sapply(xx, nrow) > 0]
    X <- do.call(cbind, xx)

    if (!is.null(X) && ncol(X) >= 4) {
      X <- X[, colMeans(is.na(X)) < 0.5, drop = FALSE]
      X <- X[rowMeans(is.na(X)) < 0.5, , drop = FALSE]
      all.X[[mystudy]] <- X
    }
  }
}
#' Run GCEA (Gene Concept Enrichment Analysis) using GSEA
#'
#' This function performs Gene Concept Enrichment Analysis (GCEA) using Gene Set Enrichment Analysis (GSEA).
#' It combines gene set analysis and concept analysis to identify enriched gene sets and concepts.
#'
#' @param X The gene expression matrix or pre-ranked gene list.
#' @param sets The gene set collection for gene set analysis.
#' @param concepts The gene concept collection for concept analysis.
#' @param Y The class labels or binary matrix indicating the phenotypes for phenotype comparison mode.
#' @param fdr The false discovery rate (FDR) threshold for selecting significant gene sets. Default is 1.
#' @param set.min The minimum gene set size. Default is 15.
#' @param set.max The maximum gene set size. Default is 500.
#' @param topgs The number of top gene sets to report. Default is 100.
#' @param nperm The number of permutations for estimating the statistical significance. Default is 1000.
#' @param rpt.label The label for the analysis report. Default is "analysis".
#' @param output.dir The output directory for saving the results. Default is NULL (results are not saved).
#' @param xmx The maximum memory allocation for the GSEA program in gigabytes (GB). Default is 10.
#' @param scoring The scoring scheme for ranking genes. Default is "weighted".
#' @param clean.files Logical value indicating whether to clean intermediate files. Default is TRUE.
#'
#' @return A list containing the results of the GCEA analysis, including gene set analysis results (gsea) and concept analysis results (gcea).
#'
#' @export
justGCEA <- function(X, sets, concepts, Y = NULL, fdr = 1, set.min = 15, set.max = 500,
                     topgs = 100, nperm = 1000, rpt.label = "analysis", output.dir = NULL,
                     xmx = 10, scoring = "weighted", clean.files = TRUE) {
  output.dir1 <- output.dir2 <- NULL
  if (!is.null(output.dir)) {
    output.dir1 <- paste(output.dir, "/gsea_out", sep = "")
    output.dir2 <- paste(output.dir, "/gcea_out", sep = "")
  }
  cat(">>> GCEA gene concept enrichment analysis <<<\n")
  cat("running gene set analysis\n")
  res1 <- justGSEA(X,
    gmt = sets, Y = Y, output.dir = output.dir1,
    topgs = topgs, nperm = nperm, rpt.label = rpt.label,
    xmx = xmx, scoring = scoring, fdr = fdr,
    set.min = set.min, set.max = set.max
  )
  cat("running concept analysis\n")
  res2 <- justGSEA(res1$NES,
    gmt = concepts, output.dir = output.dir2,
    topgs = topgs, nperm = nperm, rpt.label = rpt.label,
    xmx = xmx, scoring = scoring, fdr = fdr,
    set.min = set.min, set.max = set.max
  )
  res <- c()
  res$gsea <- res1
  res$gcea <- res2
  res
}
#' @export
cutGraph <- function(g, n = 2, k = 5, max.wt = 9999) {
  ## Cluster graph and cut crossing edges if requested.
  idx <- itercluster_louvain(g, n = n)
  g <- graph.cut_crossings(g, idx, max.wt = 9999)
  clust <- igraph::cluster_louvain(g)
  cmp <- clust$membership
  cmp.size <- table(cmp)
  cmp.top <- as.integer(Matrix::head(names(cmp.size)[order(-cmp.size)], k))
  igraph::induced_subgraph(g, which(cmp %in% cmp.top))
}
#' @title Dark Mode for Plotly Plots
#'
#' @description This function applies a dark mode theme to a plotly plot object.
#'
#' @param p A plotly plot object.
#' @param dim The number of dimensions of the plot (2 or 3).
#'
#' @details This function modifies the layout of the input plotly
#' plot object to apply a dark mode theme. It changes the background
#' color of the plot and paper to a dark color, and sets the font and
#' axis colors to a light gray color.
#'
#' @return A modified plotly plot object with the dark mode theme applied.
#'
#' @export
myplot_ly <- function(..., theme = "default") {
  ## 'Themed' plotly
  ##
  ##
  if (theme == "default") {
    p <- plotly::plot_ly(...)
  } else if (theme == "dark") {
    font.par <- list(
      color = "#FFF"
    )
    axis.par <- list(
      color = "#FFF",
      linecolor = "#FFF"
    )
    p <- plotly::plot_ly(...) %>%
      plotly::layout(
        plot_bgcolor = "rgb(10,30,50)",
        paper_bgcolor = "rgb(10,30,50)",
        xaxis = axis.par,
        yaxis = axis.par,
        title = font.par
      )
  }
  return(p)
}
#' Leading Edge Heatmap
#'
#' Generates a heatmap visualization of the leading edge analysis results.
#'
#' @param gsea An object containing GSEA results.
#' @param maxrow Maximum number of rows to display in the heatmap.
#' @param maxcol Maximum number of columns to display in the heatmap.
#' @param gsets A character vector specifying the gene sets to include in the analysis.
#' @param render The rendering method for the heatmap. Options are "gx.heatmap", "d3heatmap", or "heatmaply".
#' @param info.text Logical indicating whether to display additional information text.
#' @export
#' @return NULL if the leading edge matrix is too small (less than or equal to 1 row or 1 column),
#' otherwise, it returns the generated heatmap visualization.
gseaLeadingEdgeHeatmap <- function(gsea, maxrow = 60, maxcol = 60, gsets = NULL,
                                   render = "gx.heatmap", info.text = TRUE) {
  if (info.text) {
    cat("<br><h4>Leading Edge heatmap</h4>")
    cat("<br>Click here for the <a href='", gsea$LE_heatmap_clustered, "'>full Leading Edge heatmap</a>. ")
    cat("Complete Leading Edge analysis results are <a href='", gsea$LE_index, "'>here</a>")
  }
  LEmat <- gsea$LE_matrix
  if (!is.null(gsets)) {
    kk <- intersect(gsets, rownames(gsea$LE_matrix))
    LEmat <- LEmat <- gsea[kk, ]
  }
  rnk <- as.matrix(read.table(gsea$rnk_file, row.names = 1))[, 1]
  nes <- gsea$report[, "NES"]
  names(nes) <- gsea$report$NAME
  LEmat[which(is.na(LEmat))] <- 0
  LEmat <- LEmat[head(order(-rowMeans(LEmat != 0)), maxrow), , drop = FALSE]
  LEmat <- LEmat[, head(order(-colMeans(LEmat != 0)), maxcol), drop = FALSE]
  LEmat <- LEmat * abs(nes[rownames(LEmat)])
  LEmat <- t(t(LEmat) * (rnk[colnames(LEmat)]))
  LEmat <- abs(LEmat)**0.5 * sign(LEmat)
  LEmat[which(is.na(LEmat))] <- 0
  if (nrow(LEmat) <= 1 || ncol(LEmat) <= 1) {
    return(NULL)
  }

  if (render == "d3heatmap") {
    d3heatmap::d3heatmap(LEmat, yaxis_width = 600)
  } else if (render == "heatmaply") {
    heatmaply::heatmaply(LEmat)
  } else {
    playbase::gx.heatmap(LEmat - 1e-8,
      scale = "none",
      mar = c(8, 20), cexCol = 0.9, cexRow = 0.9,
      keysize = 0.4, key = FALSE
    )
  }
}
#' Plot copy number alteration heatmap
#'
#' @title Plot CNA heatmap
#'
#' @param ngs An NGS object containing copy number data
#' @param res CNA segmentation results from pgx.segmentCN()
#' @param annot Data frame with sample annotations
#' @param pca.filter Filter samples by PCA clustering (-1 to disable)
#' @param lwd Line width for chromosome lines
#' @param downsample Downsample CNA matrix for plotting (integer factor)
#' @param order.by Order samples by "clust" or "annot"
#' @param clip Clip copy number ratio values (0 to disable)
#' @param lab.cex Label cex size
#'
#' @return None. Plot is produced as a side-effect.
#'
#' @description
#' Generates a heatmap to visualize copy number alterations across samples.
#'
#' @details
#' This function takes CNA segmentation results from pgx.segmentCN() and
#' plots a heatmap of the copy number ratios. Samples can be ordered by
#' annotations or clustering. The full CNA matrix can be downsampled for
#' easier visualization. Chromosome lines and sample annotations are added.
#'
#' @export
pgx.plotCNAHeatmap <- function(ngs, res, annot = NA, pca.filter = -1, lwd = 1,
                               downsample = 10,
                               order.by = "clust", clip = 0, lab.cex = 0.6) {
  cna <- res$cna
  chr <- res$chr
  chr <- as.character(chr)
  pos <- res$pos

  ## ---------------------------------------------------------------------
  ## Downsample if needed
  ## ---------------------------------------------------------------------
  if (downsample > 1) {
    ## Downsample
    cat("downsampling CNA matrix...\n")
    n <- downsample
    jj <- as.vector(sapply(1:((nrow(cna) + n) / n), rep, n))[1:nrow(cna)]
    cna <- apply(cna, 2, function(x) tapply(x, jj, mean))
    gg <- tapply(rownames(res$cna), jj, paste, collapse = ",")
    rownames(cna) <- gg
    j1 <- which(!duplicated(jj))
    chr <- chr[j1]
    pos <- tapply(pos, jj, mean)
  }

  ## ---------------------------------------------------------------------
  ## take out small groups/chromsomes
  ## ---------------------------------------------------------------------
  ii <- which(chr %in% names(which(table(chr) > 3)))
  cna <- cna[ii, ]
  pos <- pos[ii]
  chr <- chr[ii]

  ## ensure order on chrpos
  ichr <- as.integer(sub("X", 23, sub("Y", 24, sub("chr", "", chr))))
  jj <- order(ichr, pos)
  cna <- cna[jj, ]
  chr <- chr[jj]
  pos <- pos[jj]

  cna <- cna - rowMeans(cna, na.rm = TRUE)
  cna <- cna / max(abs(cna), na.rm = TRUE)
  cna <- tanh(1.3 * cna)
  cna <- t(t(cna) - apply(cna, 2, median))

  if (pca.filter > 0) {
    k <- 20
    k <- pca.filter
    k <- ceiling(min(0.33 * ncol(cna), k))
    sv <- irlba::irlba(cna, nv = k)
    cna2 <- sv$u[, 1:k] %*% diag(sv$d[1:k]) %*% t(sv$v[, 1:k])
    colnames(cna2) <- colnames(cna)
    rownames(cna2) <- rownames(cna)
    cna <- cna2
  }

  hc <- NULL
  sv1 <- NULL
  if (order.by == "pc1") {
    ## by default order on SV1
    sv1 <- irlba::irlba(cna, nv = 1)$v[, 1]
    jj <- order(sv1)
    cna <- cna[, jj]
    sv1 <- sv1[jj]
  } else {
    ## order by hierarchical clustering
    jj <- Matrix::head(order(-apply(cna, 1, sd)), 1000)
    hc <- fastcluster::hclust(dist(t(cna[jj, ])), method = "ward.D2")
    cna <- cna[, hc$order]
  }

  ## create annation matrix
  ann.mat <- NULL
  if (!is.null(annot)) {
    if (is.na(annot)) {
      k <- c(grep("cell.type|tissue|cluster|group",
        colnames(ngs$samples),
        ignore.case = TRUE
      ), 1)[1]
    } else {
      k <- match(annot, colnames(ngs$samples))
    }
    k
    y <- as.character(ngs$samples[colnames(cna), k])
    ny <- length(setdiff(unique(y), NA))
    if (ny >= 2) {
      y[is.na(y)] <- "_"
      ann.mat <- model.matrix(~ 0 + y)
      colnames(ann.mat) <- sub("^y", "", colnames(ann.mat))
      rownames(ann.mat) <- colnames(cna)
    }
  }

  BLUERED2 <- colorRampPalette(c("blue3", "white", "red3"))

  ## ---------- do plotting ------------

  par(mgp = c(0.8, 0.4, 0))
  wa <- 0.1
  if (!is.null(ann.mat)) wa <- 0.05 + 0.016 * ncol(ann.mat)
  plotly::layout(matrix(1:3, 1, 3), widths = c(0.2, 0.7, wa))

  if (!is.null(hc)) {
    par(mar = c(8, 2, 12, 0))
    plot(as.dendrogram(hc),
      horiz = TRUE, leaflab = "none",
      yaxs = "i", xaxt = "n", yaxt = "n"
    )
  } else if (!is.null(sv1)) {
    par(mar = c(8, 3, 12, 0.3))
    barplot(sv1,
      horiz = TRUE, border = NA, col = "grey50", width = 0.1,
      space = 0, yaxs = "i", xaxt = "n"
    )
    mtext("PC1", side = 2, cex = 0.8)
  } else {
    frame()
  }

  ## main heatmap
  par(mar = c(8, 0.2, 12, 0))
  cna0 <- cna
  cna0 <- tanh(3 * cna0)
  cna0[which(abs(cna0) < clip)] <- NA
  Matrix::image(1:nrow(cna), 1:ncol(cna), cna0[, ],
    col = BLUERED2(16),
    ylab = "samples", xlab = "DNA copy number  (log2R)",
    yaxt = "n", yaxs = "i", xaxt = "n", xaxs = "i",
    zlim = c(-1, 1) * 1.0
  )

  ichr <- as.integer(sub("X", 23, sub("Y", 24, sub("chr", "", chr))))
  chrbrk <- which(diff(ichr) != 0)
  chrmid <- c(0, chrbrk) + diff(c(0, chrbrk, nrow(cna))) / 2
  abline(v = chrbrk, col = "grey50", lty = 1, lwd = lwd)
  chrlen <- length(unique(chr))
  j0 <- seq(1, chrlen, 2)
  j1 <- seq(2, chrlen, 2)
  mtext(unique(chr)[j0], side = 3, at = chrmid[j0], cex = lab.cex, line = 0.25)
  mtext(unique(chr)[j1], side = 3, at = chrmid[j1], cex = lab.cex, line = 0.9)

  if (!is.null(ann.mat)) {
    par(mar = c(8, 0.5, 12, 2))
    Matrix::image(1:ncol(ann.mat), 1:nrow(ann.mat), t(ann.mat),
      col = rev(grey.colors(2)), xlab = "", ylab = "",
      yaxt = "n", yaxs = "i", xaxt = "n", xaxs = "i"
    )
    mtext(colnames(ann.mat),
      side = 3, at = 1:ncol(ann.mat),
      las = 3, cex = lab.cex, line = 0.25
    )
  } else {
    frame()
  }

  ## done plotting
}
#' @export
hamming.min <- function(a, b) {
  min(rowSums(tagged.hamming(a, b)))
}
#' @title Get gene context from PubMed
#'
#' @description
#' Retrieves PubMed article information related to a gene and keyword context.
#'
#' @param gene Gene name or symbol.
#' @param keyword Optional keyword or phrases to search for context.
#'
#' @details
#' This function searches PubMed for articles related to the input \code{gene},
#' optionally filtered by \code{keyword} context phrases.
#'
#' It first searches PubMed for the gene name and synonyms. It then filters these
#' articles based on matches to the \code{keyword} terms.  The PubMed ids (PMIDs)
#' of matching articles are returned.
#'
#' Statistics are calculated to assess enrichment of the \code{keyword} within articles
#' matching the \code{gene}. These include a Fisher's exact test p-value, and word
#' co-occurrence statistics.
#'
#' @return
#' A list containing:
#' \itemize{
#' \item rifs: Character vector of PMIDs for articles matching the gene and keyword filter.
#' \item table: Contingency table used for calculating p-value.
#' \item p.value: Fisher's exact test p-value.
#' \item context: Word context scores.
#' }
#'
#' @export
pmid.getGeneContext <- function(gene, keyword) {
  gene1 <- c(gene, sub("([0-9])", "-\\1", gene))
  gene1 <- paste0("^", gene1, "$|^", gene1, "[-]")
  gene1 <- paste(gene1, collapse = "|")

  if (gene %in% biomaRt::keys(org.Hs.eg.db::org.Hs.egALIAS2EG)) {
    gname <- get(get(gene, org.Hs.eg.db::org.Hs.egALIAS2EG), org.Hs.eg.db::org.Hs.egGENENAME)
    gname <- gsub("[, -]", ".", gname)
    gene1 <- paste0(gene1, "|", gname)
  } else if (gene %in% biomaRt::keys(org.Mm.eg.db::org.Mm.egALIAS2EG)) {
    gname <- get(get(gene, org.Mm.eg.db::org.Mm.egALIAS2EG), org.Mm.eg.db::org.Mm.egGENENAME)
    gname <- gsub("[, -]", ".", gname)
    gene1 <- paste0(gene1, "|", gname)
  }

  rif.words <- colnames(GENERIF.MATRIX)
  ii <- grepl(gene1, rif.words, ignore.case = TRUE)
  match0 <- rowSums(GENERIF.MATRIX[, ii, drop = FALSE]) > 0
  match1 <- rep(1, length(match0))
  if (length(keyword) > 0) {
    i <- 1
    for (i in 1:length(keyword)) {
      jj <- grepl(keyword[i], rif.words, ignore.case = TRUE)
      match1 <- match1 & (rowSums(GENERIF.MATRIX[, jj, drop = FALSE]) > 0)
    }
  }
  sel <- ((match0 * match1) > 0)
  rif.hits <- rownames(GENERIF.MATRIX)[sel]
  rif.hits <- rif.hits[!duplicated(rif.hits)]

  ## calculate P-value for this keyword
  A <- table(gene = match0, keyword = match1)
  pv <- NA
  if (nrow(A) == 2 && ncol(A) == 2) {
    pv <- fisher.test(A, alternative = "greater")$p.value
  }

  context1 <- NULL
  match2 <- ((match0 * match1) > 0)
  m0 <- Matrix::colSums(match2 * (GENERIF.MATRIX != 0))
  m1 <- Matrix::colSums(GENERIF.MATRIX != 0)
  pp <- corpora::fisher.pval(m0, sum(match2) + 1, m1, nrow(GENERIF.MATRIX) + 1, alternative = "greater")
  pp <- sort(pp)
  qq <- p.adjust(pp)
  qq <- sort(qq)
  context1 <- Matrix::head(qq[qq < 1], 100)

  out <- list(rifs = rif.hits, table = A, p.value = pv, context = context1)
  return(out)
}
#' Share axis labels in plot grid
#'
#' @param plotList A list of ggplot objects to arrange in a grid
#' @param nrow Number of rows in the plot grid
#'
#' @return A combined ggplot object with shared axis labels
#'
#' @description Arranges multiple plots in a grid sharing common axis labels.
#'
#' @details This function takes a list of ggplot objects and arranges them in a grid
#' layout with a specified number of rows. It shares common x and y axis labels among
#' the plots to avoid repetition.
#'
#' Axis labels are removed from all plots except the bottom row and leftmost column.
#' The plots are then combined using cowplot::plot_grid() with shared axes.
#'
#' Useful for creating multi-panel figure grids with clean axis labels.
#'
#' @export
plot_grid.sharedAxisLabels <- function(plotList, nrow) {
  np <- length(plotList)
  np
  ncol <- ceiling(np / nrow)
  ncol
  ann.y <- which((0:(np - 1) %% ncol) != 0)
  ann.x <- which((0:(np - 1) %/% ncol) != (nrow - 1))
  for (i in ann.y) {
    plotList[[i]] <- plotList[[i]] + ggplot2::ylab("")
  }
  for (i in ann.x) {
    plotList[[i]] <- plotList[[i]] + ggplot2::xlab("")
  }
  cowplot::plot_grid(plotlist = plotList, nrow = nrow, labels = NA)
}
#' Signal to noise ratio test
#'
#' @title Signal to noise ratio test
#'
#' @description Performs a gene-wise signal to noise ratio test to detect differential expression.
#' Permutation test is used to compute p-values.
#'
#' @param X Numeric gene expression matrix with genes in rows and samples in columns.
#' @param y Factor of sample groups or classes. Must have exactly two levels.
#' @param ref.class Reference group name to use as baseline.
#' @param nperm Number of permutations to compute p-values.
#'
#' @details This function calculates the signal to noise ratio for each gene.
#' It then performs a permutation test by permuting the sample labels to generate a null distribution.
#' P-values are computed as the fraction of permuted statistics that are more extreme than the observed statistic.
#'
#' @return Named vector of p-values for each gene.
#'
#' @export
gx.snrtest <- function(X, y, ref.class, nperm = 200) {
  ## http://software.broadinstitute.org/gsea/doc/GSEAUserGuideFrame.html
  this.X <- X
  this.y <- y
  calc.snr <- function(this.X, this.y) {
    j1 <- which(this.y != ref.class)
    j0 <- which(this.y == ref.class)
    ma <- rowMeans(this.X[, j1], na.rm = TRUE)
    mb <- rowMeans(this.X[, j0], na.rm = TRUE)
    sa <- rowSums((this.X[, j1] - ma)**2 / (length(j1) - 1))**0.5
    sb <- rowSums((this.X[, j0] - mb)**2 / (length(j0) - 1))**0.5
    sa <- max(0.2 * abs(ma), sa, 0.2)
    sb <- max(0.2 * abs(mb), sb, 0.2)
    (ma - mb) / (sa + sb)
  }
  x0 <- calc.snr(X, y)
  pv <- x0 * 0
  i <- 1
  for (i in 1:nperm) {
    x1 <- calc.snr(X, sample(y))
    pv <- pv + 1 / nperm * (x1 > x0 * sign(x0))
  }
  pv <- pmax(pv, 1.0 / nperm)
  cat("\n")
  Matrix::head(sort(pv))
  pos.class <- setdiff(y, ref.class)[1]
  logFC <- ma - mb
  qv <- p.adjust(pv, method = "fdr")
  res <- data.frame(logFC, snr, pv, qv, ma, mb)
  colnames(res) <- c(
    "logFC", "SNR", "P.Value", "adj.P.Val",
    paste0("AveExpr.", pos.class),
    paste0("AveExpr.", ref.class)
  )
  res
}
#' @export
pgx.scTestDifferentialExpression <- function(counts, y, is.count = TRUE, samples = NULL,
                                             ref = NULL, zth = 0.25, cpm.total = 1e4,
                                             do.nzero = TRUE) {
  X1 <- NULL
  if (FALSE && is.null(cpm.total)) {
    cpm.total <- mean(Matrix::colSums(counts > 0)) ## mean number of expressing genes
    message("setting column sums to total = ", round(cpm.total, 2))
  }

  if (is.count && inherits(counts, "dgCMatrix")) {
    message("input matrix is counts (sparseMatrix)")
    ## normalize to CPM
    X1 <- counts
    X1@x <- cpm.total * X1@x / rep.int(Matrix::colSums(X1), diff(X1@p)) ## fast divide by columns sum
    ## compute log-expression
    X1@x <- log2(1 + X1@x)

    zth <- log2(1 + zth)
  } else if (is.count) {
    message("input matrix is counts (Matrix)")
    ## normalize to CPM
    X1 <- cpm.total * t(t(counts) / (1e-8 + Matrix::colSums(counts)))
    X1 <- log2(1 + X1)

    zth <- log2(1 + zth)
  } else {
    message("input matrix is normalized log-expression")
    X1 <- counts
    counts <- pmax(2**X1 - 1, 0) ## estimate counts
  }

  if (length(setdiff(unique(y), NA)) != 2) {
    stop("phenotype must have exactly two groups")
  }
  y <- factor(y)
  if (!is.null(ref)) y <- relevel(y, ref = ref)
  y0 <- levels(y)[1]
  y1 <- levels(y)[2]
  y0
  y1

  ## non-zero matrix (all zeros treated as NA)
  Z1 <- X1
  Z1[Z1 <= zth] <- NA ## notice threshold!

  ## t-test

  m1 <- matrixTests::row_t_welch(as.matrix(X1)[, y == y1], as.matrix(X1)[, y == y0]) ## all
  m2 <- matrixTests::row_t_welch(as.matrix(Z1)[, y == y1], as.matrix(Z1)[, y == y0]) ## non-zero


  pct.x <- m2$obs.x / m1$obs.x
  pct.y <- m2$obs.y / m1$obs.y
  pct.tot <- m2$obs.tot / m1$obs.tot
  pct.diff <- (pct.x - pct.y)

  ## both k1 and k2 cannot be zero...
  hack.err <- which(m2$obs.x == 0 & m2$obs.y == 0)
  if (length(hack.err) > 0) m2$obs.x[hack.err] <- 1


  pv <- corpora::fisher.pval(m2$obs.x, m1$obs.x, m2$obs.y, m1$obs.y)

  if (length(hack.err) > 0) m2$obs.x[hack.err] <- 0
  pct <- cbind(pct.x, pct.y, pct.diff, pct.tot, pct.pvalue = pv)

  kk <- c("obs.x", "obs.y", "obs.tot", "mean.x", "mean.y", "mean.diff", "pvalue")
  m1 <- m1[, kk]
  m2 <- m2[, kk]
  colnames(m2) <- paste0("nzero.", sub("^nzero[.]", "", colnames(m2)))

  b1 <- NULL
  if (!is.null(samples)) {
    ## geometric average (in log-space) by sample
    S1 <- tapply(1:ncol(X1), samples, function(i) rowMeans(X1[, i, drop = FALSE]))
    S1 <- do.call(cbind, S1)
    yb <- tapply(y, samples, function(y1) names(which.max(table(y1))))
    b1 <- matrixTests::row_t_welch(S1[, yb == y1], S1[, yb == y0]) ## all
    b1 <- b1[, kk]
    colnames(b1) <- paste0("sample.", colnames(b1))
  }

  b2 <- NULL
  if (!is.null(samples)) {
    ## sum counts then to average (in log-space) by sample
    C1 <- pmax(2**X1 - 1, 0)
    S2 <- tapply(1:ncol(C1), samples, function(i) rowMeans(C1[, i, drop = FALSE]))
    S2 <- do.call(cbind, S2)

    S2 <- log2(1 + S2)
    yb <- tapply(y, samples, function(y1) names(which.max(table(y1))))
    b2 <- matrixTests::row_t_welch(S2[, yb == y1], S2[, yb == y0]) ## all
    b2 <- b2[, kk]
    colnames(b2) <- paste0("sample2.", colnames(b2))
    remove(C1)
  }

  ## Group fold-change as used by Seurat: difference of
  ## log-average. Averages are done in linear count space. Note: no
  ## p-values are computed (not possible).
  grpx <- log2(1 + rowMeans(X1[, y == y1]))
  grpy <- log2(1 + rowMeans(X1[, y == y0]))
  fc0 <- grpx - grpy
  nz.grpx <- log2(1 + rowMeans(Z1[, y == y1], na.rm = TRUE))
  nz.grpy <- log2(1 + rowMeans(Z1[, y == y0], na.rm = TRUE))
  nz.fc0 <- nz.grpx - nz.grpy
  grp <- cbind(
    x = grpx, y = grpy, diff = fc0,
    nzero.x = nz.grpx, nzero.y = nz.grpy,
    nzero.diff = nz.fc0
  )
  colnames(grp) <- paste0("group.", colnames(grp))
  rownames(grp) <- rownames(X1)

  ## Bulk-like fold-change: sum up all counts in one group,
  ## normalize with CPM, then do difference in the log. Note: no
  ## p-values are computed (not possible).
  bulk.x <- rowSums(counts[, y == y1])
  bulk.y <- rowSums(counts[, y == y0])
  bulk.x <- bulk.x / sum(bulk.x) * cpm.total
  bulk.y <- bulk.y / sum(bulk.y) * cpm.total
  bulk.fc <- bulk.x - bulk.y
  grp2 <- cbind(x = bulk.x, y = bulk.y, diff = bulk.fc)
  colnames(grp2) <- paste0("bulk.", colnames(grp2))
  rownames(grp2) <- rownames(counts)


  df <- cbind(m1, pct, grp, grp2)
  if (do.nzero) df <- cbind(df, m2)
  if (!is.null(samples)) df <- cbind(df, b1, b2)


  P <- df[, grep("pvalue", colnames(df))]
  P[is.na(P)] <- 1
  meta.p <- apply(P, 1, function(p) metap::sumlog(p)$p)

  meta.q <- p.adjust(meta.p)
  df$meta.pvalue <- meta.p
  df$meta.qvalue <- meta.q

  ## simplify column names
  colnames(df) <- sub("mean.diff", "diff", colnames(df))

  df
}
#' @export
cosine_similarity.EXACT <- function(X) {
  vec.cos <- function(a, b) {
    j <- which(!is.na(a) & !is.na(b))
    as.numeric(sum(a[j] * b[j]) / sum(a[j]**2)**0.5 / sum(b[j]**2)**0.5)
  }
  M <- matrix(NA, ncol(X), ncol(X))
  for (i in 1:ncol(X)) {
    for (j in i:ncol(X)) {
      M[i, j] <- M[j, i] <- vec.cos(X[, i], X[, j])
    }
  }

  M
}
#' @title Matrix rownames to HUGO gene symbols
#'
#' @description
#' Converts the rownames of a matrix to official HUGO gene symbols.
#'
#' @param x A numeric matrix with rownames to convert.
#'
#' @return The input matrix with rownames converted to HUGO symbols.
#'
#' @details This function orders the rows by standard deviation, maps the
#' rownames to HUGO symbols using the playbase::alias2hugo function,
#' removes duplicates, and returns the resulting matrix with HUGO symbols
#' as rownames.
#'
#' @examples
#' \dontrun{
#' TODO
#' }
#'
#' @export
mat2hugo <- function(x) {
  x <- x[order(-apply(x, 1, sd, na.rm = TRUE)), ]
  hx <- playbase::alias2hugo(rownames(x))
  jj <- which(!duplicated(hx))
  x <- x[jj, ]
  rownames(x) <- hx[jj]
  return(x)
}
#' @title Compute Number of Significant Genes
#'
#' @param obj A pgx object
#' @param batch A character vector of batch parameters
#' @param contrast Character vector of contrasts to test
#' @param normalization Character vector of normalization methods to try
#' @param niter Number of iterations for resampling
#' @param resample Resampling fraction
#' @param show.progress Show progress
#'
#' @description Counts number of significant genes for each combination of batch correction
#'   parameters and normalization methods. Used to optimize batch correction.
#'
#' @details This function systematically goes through different combinations of batch correction
#'   parameters and normalization methods to correct the data in obj. For each combination, it performs
#'   differential expression analysis for the specified contrasts, and counts the number of significant
#'   genes. This can be used to optimize batch correction by finding the combination that maximizes
#'   the number of significant genes.
#'
#'   The batch parameters, normalization methods, number of iterations and resampling fraction can be
#'   specified. Progress is shown if show.progress is TRUE.
#'
#' @return Returns a data frame with the number of significant genes for each parameter/method combination.
#'
#' @export
pgx._runComputeNumSig <- function(ngs, parcomb, contrast, resample = -1,
                                  normalization = NORMALIZATION.METHODS,
                                  show.progress = 1) {
  k <- "cpm"
  numsig <- c()

  for (k in normalization) {
    aX <- NULL
    if (k == "nono") aX <- log(1 + ngs$counts)
    if (k == "cpm") aX <- edgeR::cpm(ngs$counts, log = TRUE)
    if (k == "TMM") aX <- log2(1 + normalizeTMM(ngs$counts))
    if (k == "RLE") aX <- log2(1 + normalizeRLE(ngs$counts))
    if (k == "quantile") aX <- limma::normalizeQuantiles(log2(1 + ngs$counts))
    if (k == "SVA") {
      mod1 <- model.matrix(~group, data = ngs$samples)
      mod0 <- cbind(mod1[, 1])
      logcpm <- edgeR::cpm(ngs$counts, log = TRUE) ## perform SVA on logCPM
      log <- capture.output({
        suppressWarnings(sv <- sva::sva(logcpm, mod1, mod0, n.sv = NULL)$sv)
        suppressWarnings(aX <- limma::removeBatchEffect(logcpm, covariates = sv, design = mod1))
      })
    }
    if (resample > 0) {
      jj <- sample(colnames(aX), ncol(aX) * resample)
      aX <- aX[, jj]
    }

    ## solve for all combinations
    pp <- parcomb[[1]]
    if (show.progress) dbg(paste0("[", k, "]")) ## show progress
    for (pp in parcomb) {
      if (show.progress) dbg(".") ## show progress
      Y <- ngs$samples[colnames(aX), ]
      bX <- NULL
      pp
      if (pp[1] != "no_correction") {
        model.formula <- formula(paste("~ 0 + ", paste(pp, collapse = " + ")))
        bvar <- model.matrix(model.formula, data = Y)
        group <- ngs$samples[colnames(aX), "group"]
        design <- model.matrix(~group)
        log <- capture.output({
          suppressWarnings(bX <- limma::removeBatchEffect(aX, covariates = bvar, design = design))
        })
      } else {
        bX <- aX
      }
      pp1 <- paste0(k, ":", paste(pp, collapse = "+"))
      pp1
      numsig[pp1] <- pgx._computeNumSig(ngs, X = bX, contrast = contrast)
    }
  }
  if (show.progress) dbg("\n") ## show progress
  return(numsig)
}
#' Generate contrasts within strata
#'
#' @param Y Factor vector indicating groups
#' @param strata Factor vector indicating strata
#' @param ref Reference group to contrast against
#'
#' @return Contrast matrix
#'
#' @description
#' Generates contrasts between groups within strata.
#'
#' @details
#' This function takes a group factor Y and strata factor, and generates contrasts between
#' the groups within each stratum.
#'
#' First, dummy variables are created for Y and the reference group is subtracted to create
#' contrasts. Then these are crossed with the strata dummy variables to encode comparisons
#' within each stratum.
#'
#' The result is a contrast matrix with columns encoding the group comparisons within
#' each stratum.
#'
#' @export
pgx.makeStratifiedContrasts <- function(Y, strata, ref) {
  S <- model.matrix(~ 0 + strata)
  colnames(S) <- sub("^strata", "", colnames(S))

  df <- data.frame(Y, strata)
  md <- makeDirectContrasts(df, c(ref, "others"))

  M <- md$contr.matrix
  C0 <- M[, grep("strata", colnames(M), invert = TRUE), drop = FALSE]
  C1 <- M[, grep("strata", colnames(M)), drop = FALSE]
  colnames(C1) <- paste0(gsub("strata:|_vs_others", "", colnames(C1)))
  contr.matrix <- c()
  i <- 1
  for (i in 1:ncol(C0)) {
    m1 <- (C1 == 1) * C0[, i]
    colnames(m1) <- paste0(colnames(C0)[i], "@", colnames(C1))
    contr.matrix <- cbind(contr.matrix, m1)
  }

  ## construct exp matrix
  grp <- md$group
  G <- model.matrix(~ 0 + grp)
  colnames(G) <- sub("^grp", "", colnames(G))
  G <- G[, rownames(contr.matrix)]
  exp.matrix <- G %*% contr.matrix
  rownames(exp.matrix) <- rownames(Y)

  ## check levels
  sel <- (Matrix::colSums(contr.matrix == -1) > 0 &
    Matrix::colSums(contr.matrix == +1) > 0)
  contr.matrix <- contr.matrix[, sel]
  exp.matrix <- exp.matrix[, sel]


  res <- list(
    contr.matrix = contr.matrix,
    group = md$group,
    exp.matrix = exp.matrix
  )
  return(res)
}
#' @title Plot Sample Clustering
#'
#' @description This function creates a plot of the sample clustering of a
#' data matrix using one of several available dimensionality reduction methods.
#'
#' @param x A numeric matrix or data frame containing the data to be clustered.
#' @param dim The number of dimensions to use for the dimensionality reduction.
#' @param method A character vector specifying the dimensionality reduction method to use.
#'   Possible values are "tsne", "umap", and "pca".
#' @param ntop The number of top features to use for the dimensionality reduction.
#' @param ... Additional arguments passed to the `plot` function.
#'
#' @details This function first calls the `pgx.clusterMatrix` function to perform
#' dimensionality reduction on the input data using the specified method. The
#' resulting low-dimensional representation of the data is then plotted using the
#' `plot` function from the base R graphics system.
#'
#' @return The function returns invisibly.
#'
#' @export
pgx.plotSampleClustering <- function(x, dim = 2,
                                     method = c("tsne", "umap", "pca"),
                                     ntop = 1000, ...) {
  method <- method[1]
  clust <- pgx.clusterMatrix(
    x,
    perplexity = NULL,
    ntop = ntop, dims = dim,
    row.center = TRUE, row.scale = FALSE,
    find.clusters = FALSE, kclust = 1,
    prefix = "C", clust.detect = "louvain",
    method = method
  )


  plot(clust$pos2d, ...)
}
#' @describeIn  pgx.getGEOseries Retrieve and process gene expression count
#' data from the supplementary files of a GEO dataset. It takes a GEO
#' accession ID, downloads the series matrix and annotation data from the
#' supplementary files, processes it into a count matrix, and returns the result.
#' @export
pgx.getGEOcounts.fromSuppl <- function(id) {
  ## Retrieve expression matrix, phenotype and probe annotation
  ## matrices for a certain GEO id.
  ##


  ## load series and platform data from GEO
  id

  gse <- try(GEOquery::getGEO(id, GSEMatrix = FALSE, getGPL = FALSE))

  if (inherits(gse, "try-error")) {
    cat("ERROR: GEOquery::getGEO() error\n")
    return(NULL)
  }

  supp_file <- gse@header$supplementary_file
  supp_file

  return(expr)
}
#' @export
prot.filterReverseContaminantOnly <- function(x, rev = T, con = T, only = T) {
  if (rev == T) {
    x <- x[x$Reverse != "+", ]
  }
  if (con == T & sum(colnames(x) == "Potential.contaminant") == 1) {
    x <- x[x$Potential.contaminant != "+", ]
  }
  if (con == T & sum(colnames(x) == "Contaminant") == 1) {
    x <- x[x$Contaminant != "+", ]
  }
  if (only == T) {
    x <- x[x$Only.identified.by.site != "+", ]
  }
  return(x)
}
#' @export
prot.imputeMissing <- function(X, method, downshift = 1.8, width = 0.3,
                               k = 10, q = 0.01, groups = NULL, zero.na = TRUE) {
  ## in count space!!!

  ## treat zeros as missing values
  if (zero.na) X[X == 0] <- NA

  if (method == "zero") {
    X[is.na(X)] <- 0
  }
  if (method == "min") {
    X[is.na(X)] <- min(X[!is.na(X) & X > 0])
  }
  if (method == "quantile") {
    X[is.na(X)] <- quantile(X[!is.na(X) & X > 0], probs = q)[[1]]
  }
  if (method == "gaussian") {
    logX <- log(1 + X)
    impX <- .RGimputation(logX, width = width, downshift = downshift, bycolumn = TRUE)
    X <- exp(pmax(impX, 0)) - 1
  }
  if (method == "group.median") {
    if (is.null(groups)) stop("'group.median' method needs groups")
    X <- prot.medianImpute(X, groups)
  }
  if (method == "nmf") {
    if (is.null(groups)) stop("'nmf' method needs groups")
    X <- prot.nmfImpute(X, groups, k = k)
  }
  if (method == "none") {
    ## nothing
  }
  return(X)
}
#' @describeIn eset.getPhenoData Extracts the title from the phenotype data of an ExpressionSet object.
#' @export
eset.getTitle <- function(eset) as.character(pData(phenoData(eset))$title)




#' @describeIn eset.getPhenoData parses the "Characteristics" metadata field from
#' an ExpressionSet into a data frame by splitting the text on a delimiter
#' like "," or ";".
#' @param ch "Characteristics" text string
#' @param split delimiter character to split on `c(",", ";", "\\|", "_", " ")`.
#' @export
eset.parseCharacteristicsInfo <- function(ch, split = ",") {
  terms <- sub(":.*$", "", trimws(strsplit(ch, split = ",")[[1]]))
  value <- t(sapply(ch, function(s) (strsplit(s, split = ",")[[1]])))
  value <- apply(value, 2, function(s) trimws(sub(".*[:]", "", s)))
  value <- data.frame(value, stringsAsFactors = FALSE, check.names = FALSE)
  rownames(value) <- NULL
  colnames(value) <- terms
  return(value)
}
#' @export
decompose_graph <- function(g, grouped = NULL, method = "louvain", k = NULL,
                            cex1 = 2, cex2 = 10, nlouvain = 1) {
  if (is.null(g$layout)) stop("graph must have layout")
  get.idx <- function(g, method, k) {
    if (method == "louvain") {
      idx <- rep(1, length(igraph::V(g)))
      for (k in 1:nlouvain) {
        idx_new <- rep(0, length(idx))
        i <- 1
        for (i in unique(idx)) {
          jj <- which(idx == i)
          g1 <- igraph::induced_subgraph(g, igraph::V(g)[jj])
          idx1 <- igraph::cluster_louvain(g1, weights = E(g1)$weight)$membership
          idx_new[jj] <- max(idx_new) + idx1
        }
        idx <- idx_new
        cat("louvain: n=", length(unique(idx)), "\n")
      } ## end of nlouvain iterations
    } else if (method == "hclust") {
      if (is.null(k)) stop("hclust needs k")

      D <- as.matrix(1.0 - g[, ])

      hc <- fastcluster::hclust(as.dist(D), method = "average")
      k1 <- max(min(k, nrow(D) / 2), 1)
      idx <- cutree(hc, k1)
      table(idx)
    } else if (method == "nclust") {
      if (is.null(k)) stop("nclust needs k")
      D <- as.matrix(1.0 - g[, ])


      hc <- nclust(as.dist(D), link = "average")
      k1 <- max(min(k, nrow(D) / 2), 1)
      idx <- cutree(hc, k1)
      table(idx)
    } else {
      stop("fatal. unknown method", method, "\n")
    }
    return(idx)
  }

  idx <- NULL
  group.idx <- NULL
  if (!is.null(grouped)) {
    grouped[is.na(grouped)] <- "(missing)"
    groups <- unique(grouped)
    table(grouped)
    ng <- length(igraph::V(g))
    idx <- rep(0, ng)
    group.idx <- list()
    i <- 1
    for (i in 1:length(groups)) {
      jj <- which(grouped == groups[i])
      sub.idx <- rep(1, length(jj))
      length(jj)
      if (length(jj) > 3) {
        subg <- igraph::induced_subgraph(g, jj)
        k0 <- max(ceiling(k * (length(jj) / ng)), 2)
        k0
        sub.idx <- get.idx(subg, method = method, k = k0)
        table(sub.idx)
      }
      table(sub.idx)
      sub.idx <- sub.idx + max(idx)
      idx[jj] <- sub.idx
      group.idx[[i]] <- unique(sub.idx)
    }
  } else {
    idx <- get.idx(g, method = method, k = k)
  }
  table(idx)

  ## resort indices
  idx <- as.integer(factor(idx, levels = order(-table(idx))))
  names(idx) <- igraph::V(g)$name
  table(idx)

  ## adjacency matrix
  M <- as.matrix(g[, ])
  diag(M) <- 1
  M[1:4, 1:4]

  ## downsample adjacency matrix
  ngrp <- length(unique(idx))
  ngrp
  subM <- tapply(1:ncol(M), idx, function(ii) rowMeans(M[, ii, drop = FALSE], na.rm = TRUE))
  subM <- lapply(subM, function(x) tapply(x, idx, mean))
  subM <- matrix(unlist(subM), ncol = ngrp, nrow = ngrp)
  sum(is.na(subM))
  rownames(subM) <- colnames(subM) <- paste0("cluster.", 1:nrow(subM))
  M[1:4, 1:4]
  subM[1:4, 1:4]

  ## create downsamples graph
  subG <- igraph::graph_from_adjacency_matrix(
    subM,
    weighted = TRUE, diag = FALSE, mode = "undirected"
  )
  igraph::V(subG)$name <- rownames(subM)
  subG$members <- tapply(idx, idx, names)
  names(subG$members) <- rownames(subM)


  avg.pos <- apply(g$layout[names(idx), ], 2, function(x) tapply(x, idx, mean))
  rownames(avg.pos) <- rownames(subM)
  subG$layout <- avg.pos

  subgraphs <- list()
  for (i in 1:length(igraph::V(subG))) {
    vv <- subG$members[[i]]
    child <- igraph::induced_subgraph(g, vv)
    child$layout <- g$layout[vv, ]
    subgraphs[[i]] <- child
  }

  group.idx2 <- rep(1, length(subgraphs))
  if (!is.null(grouped)) {
    group.idx2 <- sapply(1:length(subgraphs), function(i) {
      which(sapply(group.idx, function(x) (i %in% x)))
    })
  }

  subG <- igraph::simplify(subG)
  subG$group.idx <- group.idx2
  subG$membership <- idx
  subG$subgraphs <- subgraphs
  return(subG)
}
#' Convert expression matrix to contrast matrix
#'
#' @param exp.matrix Expression matrix with samples in columns
#'
#' @return A list containing:
#' \itemize{
#' \item{contr.matrix}{Contrast matrix}
#' \item{group}{Factor of sample groups}
#' \item{exp.matrix}{Original expression matrix}
#' }
#'
#' @description
#' Converts an expression matrix into a contrast matrix encoding comparisons between sample groups.
#'
#' @details
#' This function takes an expression matrix and generates a contrast matrix from it.
#'
#' It creates a group factor by pasting column names. The number of unique groups is used to name contrast columns as "groupN".
#'
#' The contrast matrix contains a row for each unique group, with 0/1/-1 values encoding the sample membership.
#'
#' The original group factor and expression matrix are also returned.
#'
#' @export
expmat2contrast <- function(exp.matrix) {
  group <- apply(exp.matrix, 1, paste, collapse = "_")
  n.group <- length(unique(group))
  group <- factor(group)
  if (ncol(exp.matrix) > 3) {
    levels(group) <- paste0("group", 1:n.group)
  }
  sel <- !duplicated(group)
  contr.matrix <- exp.matrix[sel, ]
  rownames(contr.matrix) <- group[sel]
  out <- list(contr.matrix = contr.matrix, group = group, exp.matrix = exp.matrix)
  return(out)
}
#' @export
tcosine.sparse.paral <- function(X, th = 0.01, k = 100, mc.cores = 4, blocksize = 100,
                                 verbose = 1) {
  ##
  ## Normalizing X then doing tcrossprod is same as doing
  ## tcosine_similarity on non-normalized X.
  ##
  X <- X / (1e-20 + Matrix::rowSums(X**2)**0.5)
  iter <- ceiling(nrow(X) / blocksize)
  sim <- c()
  i <- 1
  for (i in 1:iter) {
    down <- i * blocksize - blocksize + 1
    up <- min(i * blocksize, nrow(X))
    if (verbose) print(paste0(down, " - ", up))
    j <- 1
    s <- parallel::mclapply(down:up, function(j) {
      cr <- tcrossprod(X[j, ], X)[1, ]
      jj <- as.vector(which(cr >= th))
      jj <- Matrix::head(jj[order(cr[jj], decreasing = TRUE)], k)
      if (length(jj) > 0) {
        data.frame(from = rep(j, length(jj)), to = jj, val = as.numeric(cr[jj]))
      } else {
        data.frame(from = NA, to = NA, val = NA)[0, ]
      }
    }, mc.cores = mc.cores)
    sim[down:up] <- s
  }

  sim1 <- do.call(rbind, sim)
  dim(sim1)
  S <- Matrix::sparseMatrix(i = sim1[, 1], j = sim1[, 2], x = sim1[, 3])
  if (verbose) {
    sp <- round(100 * mean(S != 0, na.rm = TRUE), digits = 2)
    sp.na <- round(100 * mean(is.na(S), na.rm = TRUE), digits = 2)
    cat("tcosine.sparse: sparsity=", sp, "%\n")
    if (sp.na > 0) cat("tcosine.sparse: NA=", sp.na, "%\n")
  }
  return(S)
}
#' @export
knn.predict <- function(data, pred.var, x.var = NULL, K, samples = NULL, ntest = 0.33) {
  qvars <- names(data$Q)
  x.var0 <- x.var
  if (is.null(x.var) || "<all others>" %in% x.var) {
    x.var <- setdiff(qvars, pred.var)
  }
  x.var <- intersect(x.var, qvars)
  if (is.null(samples)) samples <- 1:nrow(data$Q[[1]])
  X <- do.call(cbind, data$Q[x.var])[samples, ]
  Y <- as.character(data$tcr[samples, pred.var])
  names(Y) <- rownames(data$tcr)[samples]

  if (inherits(X, "matrix")) {
    X <- Matrix::Matrix(X, sparse = TRUE)
  }

  ## random split into prediction and training samples
  cat("ntest=", ntest, "\n")
  if (ntest > 0 && ntest < 1) ntest <- round(ntest * nrow(X))
  cat("ntest=", ntest, "\n")
  cat("nrowX=", nrow(X), "\n")
  j1 <- sample(nrow(X), ntest) ## prediction
  j0 <- setdiff(1:nrow(X), j1)

  ## distance calculations for test-set (superfast!)
  pos <- dist.xy <- NULL

  if (sum(is.na(X)) > 0) {
    system.time(dist.xy <- 1 - cosine_similarity(t(X[j1, ]), t(X[j0, ])))
  } else {
    system.time(dist.xy <- 1 - qlcMatrix::cosSparse(t(X[j1, ]), t(X[j0, ])))
  }
  dim(dist.xy)

  ## Get best prediction from K-nearest neighbours
  knn <- apply(dist.xy, 1, function(x) Matrix::head(order(x), K))
  if (NCOL(knn) == 1) {
    pred.yy <- matrix(Y[j0][knn], ncol = 1)
  } else {
    pred.yy <- t(apply(knn, 2, function(i) Y[j0][i]))
  }
  colnames(pred.yy) <- paste0("nnb", 1:ncol(pred.yy))
  dim(pred.yy)
  pred.max <- apply(pred.yy, 1, function(x) names(sort(-table(x)))[1])
  names(pred.max) <- names(Y)[j1]

  ## Accuracies
  test.y <- Y[j1]
  acc <- sens <- specf <- NA
  acc <- mean(pred.max == test.y)
  acc

  res <- list(
    y.pred = pred.max, y = Y[j1],
    pred.var = pred.var, x.input = x.var0, x.var = x.var,
    train.idx = j0, test.idx = j1,
    accuracy = acc, sensitivity = sens, specificity = specf,
    pos = pos
  )
  return(res)
}
#' @title Filter Genes by Family
#'
#' @description This function filters a data frame of genes by a specified gene family.
#'
#' @param genes A data frame of genes to be filtered, with row names representing probe names.
#' @param family A character string specifying the name of the gene family to filter by.
#' @param ngs A list containing gene families, where `ngs$families` is a
#' named list of character vectors representing gene families.
#'
#' @details The function takes a data frame of genes `genes`, a gene family
#' name `family`, and a list of gene families `ngs` as input.
#' The function searches for rows in the `genes` data frame where the probe
#' name, short probe name, or gene name matches a gene in the specified gene family.
#' The row names of the matching rows are returned.
#'
#' @return A character vector representing the row names of the rows in the
#' `genes` data frame that match the specified gene family.
#'
#' @export
filterFamily <- function(genes, family, ngs) {
  gg <- ngs$families[[10]]
  gg <- ngs$families[[family]]
  ## check probe name, short probe name or gene name for match
  p0 <- (toupper(sub(".*:", "", rownames(genes))) %in% toupper(gg))
  p1 <- (toupper(rownames(genes)) %in% toupper(gg))
  p2 <- (toupper(as.character(genes$gene_name)) %in% toupper(gg))
  jj <- which(p0 | p1 | p2)
  rownames(genes)[jj]
}
#' Perform GSEA analysis for all contrasts
#'
#' This function performs GSEA (Gene Set Enrichment Analysis) analysis for all the contrasts specified in the design matrix.
#' It iterates over the contrasts, runs GSEA for each contrast, and saves the results in the specified output directory.
#'
#' @param X The expression matrix.
#' @param gmt The gene set collection in GMT format.
#' @param design The design matrix.
#' @param contr.matrix The contrast matrix.
#' @param output_dir The directory to save the GSEA results.
#' @param set.min The minimum size of gene sets to consider in the analysis.
#' @param set.max The maximum size of gene sets to consider in the analysis.
#' @param fdr The false discovery rate (FDR) threshold.
#' @param skip.done Logical value indicating whether to skip the contrasts that have already been analyzed and have existing results in the output directory.
#'
#' @export
#'
#' @return A list of GSEA analysis results for each contrast.
#'
gsea.fitAllContrasts <- function(X, gmt, design, contr.matrix, output_dir,
                                 set.min = 15, set.max = 500, fdr = 0.25, skip.done = TRUE) {
  exp.matrix <- (design %*% contr.matrix)
  length(gmt)
  comparisons <- colnames(contr.matrix)
  outputs <- list()
  k <- 1
  for (k in 1:length(comparisons)) {
    cat(">>> running GSEA for comparison: ", comparisons[k], " <<<\n")

    ## set params
    #
    comp.name <- comparisons[k]
    #
    comp.name

    ## check if already done
    sub_dir <- file.path(output_dir, comp.name)
    sub_dir
    if (dir.exists(sub_dir) && skip.done) {
      cat("comparison ", comp.name, "already done. skipping.\n")
      ## should we read the output file??
      next
    }

    sel <- which(exp.matrix[, comp.name] != 0)
    xx <- X[, sel]
    yy <- 1 * (exp.matrix[sel, comp.name] > 0)

    j1 <- which(contr.matrix[, comp.name] > 0)
    j0 <- which(contr.matrix[, comp.name] < 0)
    cls1 <- paste(rownames(contr.matrix)[j1], sep = "and")
    cls0 <- paste(rownames(contr.matrix)[j0], sep = "and")
    c(cls0, cls1)
    yy <- c(cls0, cls1)[1 + yy]
    yy
    ref <- cls0
    dim(xx)
    yy
    ref

    ## run GSEA
    outputs[[k]] <- run.GSEA(xx, yy, gmt,
      fdr = fdr,
      do.leading.edge = TRUE, clean.files = FALSE,
      output.dir = sub_dir,
      set.min = set.min, set.max = set.max, topgs = 100,
      ref.type = ref, permute = "gene_set"
    )
  } ## for-loop comparisons

  return(outputs)
}
#' @title Apply k-Nearest Neighbors Median Filter
#'
#' @description This function applies a k-nearest neighbors median filter to a vector of values.
#'
#' @param x A numeric vector of values to be filtered.
#' @param pos A matrix of positions for each element in `x`.
#' @param k An optional numeric value specifying the number of nearest neighbors to use for filtering.
#' The default value is 10.
#'
#' @details This function takes a numeric vector `x` and a matrix of positions `pos` as input and applies a k-nearest neighbors median filter to `x`.
#' The number of nearest neighbors used for filtering is specified by the `k` parameter.
#' For each element in `x`, the function finds its `k` nearest neighbors in `pos` and calculates the median of their values in `x`.
#' The resulting filtered values are returned as a numeric vector of the same length as `x`.
#'
#' @return A numeric vector of the same length as `x`, containing the filtered values.
#'
#' @export
knnMedianFilter <- function(x, pos, k = 10) {
  nb <- FNN::get.knn(pos[, ], k = k)$nn.index
  fx <- factor(x)
  mx <- matrix(fx[as.vector(nb)], nrow = nrow(nb), ncol = ncol(nb))
  x1 <- apply(mx, 1, function(x) names(which.max(table(x))))
  x1
}
#' @export
prot.readProteinGroups <- function(file, meta = NULL, sep = "\t", collapse.gene = TRUE,
                                   is.log = FALSE, use.LFQ = FALSE,
                                   filter.contaminants = TRUE) {
  ## Split data file
  message("reading proteinGroups file ", file)


  D <- data.table::fread(file, check.names = FALSE)
  D <- data.frame(D, check.names = FALSE)

  ## Filter contaminants
  contaminant.cols <- c("Reverse", "Only identified by site", "Potential contaminant")
  contaminant.cols <- intersect(contaminant.cols, colnames(D))
  contaminant.cols
  D$is.contaminant <- (rowSums(D[, contaminant.cols, drop = FALSE] == "+", na.rm = TRUE) >= 1)
  if (filter.contaminants) {
    D <- D[which(!D$is.contaminant), ]
  }

  ## parse gene annotation
  genes <- D[, c("Majority protein IDs", "Gene names", "Protein names")]
  colnames(genes) <- c("protein_id", "gene_name", "gene_title")
  gg <- as.character(genes$gene_name)
  gg <- sapply(gg, function(x) strsplit(x, split = ";")[[1]][1]) ## take just FIRST gene
  genes$gene_alias <- genes$gene_name
  genes$gene_name <- gg
  genes$is_contaminant <- D$is.contaminant

  ## give unique rownames
  rownames(D) <- paste0("tag", 1:nrow(D), ":", genes$gene_name) ## add one gene to tags
  rownames(genes) <- rownames(D)

  ## extract data blocks (use LFQ as intensity)
  counts <- D[, grep("^Intensity ", colnames(D), value = TRUE)]

  if (use.LFQ) {
    sel <- grep("^LFQ Intensity", colnames(D))
    if (length(sel) == 0) {
      stop("FATAL:: could not find LFQ Intensity data")
    }
    counts <- D[, sel, drop = FALSE]
  }
  if (is.log) {
    counts <- 2**as.matrix(counts)
  } else {
    counts <- as.matrix(sapply(counts, as.numeric)) ## from integer64
  }
  colnames(counts) <- sub("Intensity ", "", colnames(counts))
  colnames(counts) <- sub("LFQ intensity ", "", colnames(counts))
  summary(Matrix::colSums(counts, na.rm = TRUE))

  ## collapse by gene
  if (collapse.gene) {
    cat("collapsing to FIRST gene\n")
    gene <- genes$gene_name
    counts <- apply(counts, 2, function(x) tapply(x, gene, sum, na.rm = TRUE))
    prots <- tapply(as.character(genes$protein_id), gene, function(s) {
      paste(unique(unlist(strsplit(s, split = ";"))), collapse = ";")
    })
    genes <- genes[match(rownames(counts), genes$gene_name), ]
    genes$protein_id <- NULL
    genes$protein_ids <- prots[match(rownames(counts), names(prots))]
    rownames(genes) <- genes$gene_name
  }

  if (!is.null(meta) && is.character(meta)) {
    if (grepl("csv$", meta)) {
      meta <- read.csv(meta, header = TRUE)
    } else {
      meta <- read.delim(meta, header = TRUE, sep = "\t")
    }
  }
  if (is.null(meta)) {
    meta <- data.frame(sample.name = colnames(counts))
  }

  res <- c()
  res$samples <- meta
  res$genes <- genes
  res$counts <- as.matrix(counts)

  message("sample info template created but please complete ngs$samples")

  return(res)
}
#' Calculate Signal-to-Noise ratio for GSEA
#'
#' This function calculates the Signal-to-Noise ratio (SNR) for Gene Set Enrichment Analysis (GSEA).
#' It measures the difference in the means of two classes (0 and 1) relative to the standard deviation.
#' The SNR is used to rank genes and assess their relevance in the context of gene set enrichment.
#'
#' @param X The gene expression matrix.
#' @param Y The class labels or binary matrix indicating the two classes (0 and 1).
#'
#' @return A matrix of Signal-to-Noise ratios for each gene in the expression matrix.
#'
#' @export
gsea.snr <- function(X, Y) {
  if (NCOL(Y) == 1) Y <- matrix(Y, ncol = 1)
  if (ncol(X) != nrow(Y)) {
    stop("dimension mismatch")
  }
  S <- matrix(NA, nrow(X), ncol(Y))
  j <- 1
  for (j in 1:ncol(Y)) {
    yj <- Y[, j]
    mx <- cbind(
      rowMeans(X[, yj == 0, drop = FALSE], na.rm = TRUE),
      rowMeans(X[, yj == 1, drop = FALSE], na.rm = TRUE)
    )
    sx <- cbind(
      apply(X[, yj == 0, drop = FALSE], 1, sd, na.rm = TRUE),
      apply(X[, yj == 1, drop = FALSE], 1, sd, na.rm = TRUE)
    )
    sx <- pmax(pmax(sx, 0.2 * abs(mx)), 0.2)
    sx <- rowMeans(sx, na.rm = TRUE) * ncol(sx) ## robust "sum"
    S[, j] <- (mx[, 1] - mx[, 2]) / sx
  }
  colnames(S) <- colnames(Y)
  rownames(S) <- rownames(X)
  return(S)
}
#' @export
pgx.getTCGAdataset <- function(study, genes = NULL, matrix_file = NULL, from.h5 = TRUE,
                               datatype = "mrna") {
  ## For a specific TCGA study get the expression matrix and
  ## clinical data.
  ##

  ## check if H5 exists
  from.h5 <- (from.h5 && !is.null(matrix_file) && file.exists(matrix_file))



  mycgds <- cgdsr::CGDS("http://www.cbioportal.org/")

  all.studies <- sort(cgdsr::getCancerStudies(mycgds)[, 1])
  if (!all(study %in% all.studies)) {
    ss <- setdiff(study, all.studies)
    stop(ss, "is not in TCGA studies")
  }

  ## Gather data from all study
  X <- list()
  clin <- list()
  mystudy <- study[1]
  for (mystudy in study) {
    cat("getting TCGA data for", mystudy, "...\n")

    mystudy

    myprofiles <- cgdsr::getGeneticProfiles(mycgds, mystudy)[, 1]
    myprofiles

    ## datatypes
    datatype0 <- datatype
    if (datatype0 == "mrna") {
      datatype <- "_mrna"
      if (any(grepl("rna_seq_v2_mrna$", myprofiles))) datatype <- "rna_seq_v2_mrna"
      if (any(grepl("rna_seq_mrna$", myprofiles))) datatype <- "rna_seq_mrna"
      pr.datatype <- grep(paste0(datatype, "$"), myprofiles, value = TRUE)
      pr.datatype
      if (length(pr.datatype) == 0) {
        cat("WARNING:: could not find mRNA for", mystudy, "...\n")
        next()
      }
    } else {
      warning(paste("datatype", datatype0, "not yet supported"))
      return(NULL)
    }

    all.cases <- cgdsr::getCaseLists(mycgds, mystudy)[, 1]
    all.cases


    caselist <- grep(paste0(datatype, "$"), all.cases, value = TRUE)
    caselist
    samples <- NULL

    if (!is.null(genes) && !from.h5) {
      cat("downloading...\n")
      ## If only a few genes, getProfileData is a faster way
      ##
      expression <- t(cgdsr::getProfileData(mycgds, genes, pr.datatype, caselist))
      samples <- gsub("[.]", "-", colnames(expression))
      colnames(expression) <- samples
    } else {
      cat("extracting from locally stored H5 matrix...\n")
      ## For all genes, getProfileData cannot do and we use
      ## locally stored H5 TCGA data file from Archs4.
      ##
      xx <- cgdsr::getProfileData(mycgds, "---", pr.datatype, caselist)
      samples <- gsub("[.]", "-", colnames(xx))[3:ncol(xx)]

      rhdf5::h5closeAll()

      has.h5 <- file.exists(matrix_file)
      has.h5

      if (!has.h5) {
        stop("FATAL: could not find tcga_matrix.h5 matrix. Please download from Archs4.")
      } else {
        ## Retrieve information from locally stored H5 compressed data
        aa <- rhdf5::h5ls(matrix_file)[, 1:2]
        aa
        ii <- which(aa[, 1] == "/meta")[-1]
        lapply(ii, function(i) Matrix::head(rhdf5::h5read(matrix_file, paste0("/meta/", aa[i, 2]))))

        id1 <- rhdf5::h5read(matrix_file, "/meta/gdc_cases.samples.portions.submitter_id")
        id2 <- rhdf5::h5read(matrix_file, "/meta/gdc_cases.samples.submitter_id")
        id3 <- rhdf5::h5read(matrix_file, "/meta/gdc_cases.submitter_id")
        id2x <- substring(id2, 1, 15)

        h5.genes <- rhdf5::h5read(matrix_file, "/meta/genes")
        if (!is.null(genes)) h5.genes <- intersect(genes, h5.genes)
        samples <- intersect(samples, id2x)
        sample_index <- which(id2x %in% samples)
        gene_index <- 1:length(h5.genes)

        if (length(sample_index) == 0 || length(gene_index) == 0) {
          return(list(X = NULL, clin = NULL))
        }

        expression <- rhdf5::h5read(
          matrix_file, "data/expression",
          index = list(gene_index, sample_index)
        )
        rhdf5::H5close()
        colnames(expression) <- substring(id2[sample_index], 1, 15)
        rownames(expression) <- h5.genes
        expression <- expression[, order(-Matrix::colSums(expression))]
        expression <- expression[, samples]
      }
    }
    this.clin <- cgdsr::getClinicalData(mycgds, caselist)
    rownames(this.clin) <- gsub("[.]", "-", rownames(this.clin))
    this.clin <- this.clin[samples, , drop = FALSE]
    expression <- expression[, samples, drop = FALSE]
    X[[mystudy]] <- expression
    clin[[mystudy]] <- this.clin
  }

  res <- list(X = X, clin = clin)
  return(res)
}
#' @describeIn makeDirectContrasts Make direct contrasts from a factor nested function.
#' @export
makeDirectContrasts2 <- function(Y, ref, na.rm = TRUE) {
  makeDirectContrasts(Y = Y, ref = ref, na.rm = na.rm)
}
#' @export
prot.testTwoGroups <- function(X, group1, group2, method = "limma",
                               labels = NULL, gene = NULL) {
  out1 <- NULL ## important to avoid old results
  p1 <- p2 <- p3 <- NULL ## just for testing

  if (!is.null(labels) && all(c(group1, group2) %in% labels)) {
    group1 <- colnames(X)[which(labels %in% group1)]
    group2 <- colnames(X)[which(labels %in% group2)]
  }

  if (method == "geiger") {
    ## Your code
    out1 <- volcano(X[, group1], X[, group2], rownames(X))
    out1$Gene <- NULL
    colnames(out1) <- c("logFC", "P.Value")
  } else if (method %in% c("t.welch", "t.equalvar")) {
    ## faster t-test

    j1 <- which(colnames(X) %in% group1)
    j2 <- which(colnames(X) %in% group2)
    if (method == "t.welch") {
      out0 <- matrixTests::row_t_welch(X[, j1], X[, j2])
    }
    if (method == "t.equalvar") {
      out0 <- matrixTests::row_t_equalvar(X[, j1], X[, j2])
    }
    out0$pvalue[is.na(out0$pvalue)] <- 1
    out0$qvalue <- p.adjust(out0$pvalue, method = "fdr")
    out1 <- out0[, c("mean.diff", "pvalue", "qvalue")]
    colnames(out1) <- c("logFC", "P.Value", "adj.P.Val")
  } else if (method == "limma" && is.null(labels)) {
    ## See e.g. https://bioconductor.org/help/course-materials/2010/BioC2010/limma2.pdf

    jj <- which(colnames(X) %in% c(group1, group2))
    y <- 1 * (colnames(X)[jj] %in% group1)
    design <- model.matrix(~y)
    colnames(design) <- c("Intercept", "group1_vs_group2")
    fit <- limma::eBayes(limma::lmFit(X[, jj], design), trend = TRUE, robust = TRUE)
    out1 <- limma::topTable(fit, coef = 2, sort.by = "none", number = Inf)
    out1 <- out1[, c("logFC", "P.Value", "adj.P.Val")]
  } else if (method == "limma" && !is.null(labels)) {
    ## Same as above but we retain all samples in the model

    design <- model.matrix(~ 0 + labels)
    colnames(design) <- sub("labels", "", colnames(design))
    rownames(design) <- colnames(X)

    level1 <- unique(labels[match(group1, colnames(X))])
    level2 <- unique(labels[match(group2, colnames(X))])


    levels <- colnames(design)
    contr <- matrix(1 * (levels %in% level1) - 1 * (levels %in% level2), ncol = 1)
    rownames(contr) <- levels
    contr
    contr <- contr[match(colnames(design), rownames(contr)), , drop = FALSE]
    rownames(contr) <- colnames(design)
    colnames(contr)[1] <- "contrast1"
    contr[is.na(contr)] <- 0
    fit1 <- limma::lmFit(X[, ], design)
    fit2 <- limma::contrasts.fit(fit1, contr)
    fit2 <- limma::eBayes(fit2, trend = TRUE, robust = TRUE)
    out1 <- limma::topTable(fit2, coef = 1, sort.by = "none", number = Inf)
    out1 <- out1[, c("logFC", "P.Value", "adj.P.Val")]
  } else if (method == "msms.edgeR") {
    s <- colnames(X)
    grp <- c(NA, "group1", "group2")[1 + 1 * (s %in% group1) + 2 * (s %in% group2)]
    pd <- data.frame(sample = colnames(X), group = grp)
    rownames(pd) <- colnames(X)
    fd <- data.frame(proteins = rownames(X), gene = rownames(X))
    rownames(fd) <- rownames(X)
    jj <- which(!is.na(grp))
    e <- MSnSet(2**X[, jj], fd, pd[jj, ])
    null.f <- "y~1"
    alt.f <- "y~group"
    div <- apply(exprs(e), 2, sum, na.rm = TRUE)
    out1 <- msms.edgeR(e, alt.f, null.f, div = div, fnm = "group")
    out1$adj.P.Val <- p.adjust(out1$p.value, method = "fdr")
    out1$LogFC <- -out1$LogFC ## seems switched
    out1 <- out1[, c("LogFC", "p.value", "adj.P.Val")]
    colnames(out1) <- c("logFC", "P.Value", "adj.P.Val")
  } else {
    stop("ERROR:: unknown method")
  }

  ## set NaN p.values to 1??
  out1$P.Value[which(is.nan(out1$P.Value))] <- 1

  ## add FDR q.value
  if (!("adj.P.Val" %in% names(out1))) {
    out1$adj.P.Val <- p.adjust(out1$P.Value, method = "fdr")
  }

  ## good practice to add group means
  avg <- cbind(rowMeans(X[, group1], na.rm = TRUE), rowMeans(X[, group2], na.rm = TRUE))

  colnames(avg) <- paste0("avg.", c("group1", "group2"))
  out1 <- cbind(out1, avg)
  if (!is.null(gene)) {
    out1 <- cbind(gene, out1)
  }

  return(out1)
}
#' @export
pgx.createSeuratObject <- function(counts, aggr.csv = NULL,
                                   project = "SeuratProject", max.cells = 2000) {
  obj <- Seurat::CreateSeuratObject(counts, min.cells = 5, project = project)

  if (!is.null(aggr.csv)) {
    aggr <- read.csv(aggr.csv)
    sample.idx <- as.integer(sub(".*-", "", colnames(counts)))
    pheno <- aggr[sample.idx, c("library_id", "phenotype")]

    rownames(pheno) <- colnames(counts)


    obj@meta.data <- cbind(obj@meta.data, meta.data)
  }

  Seurat::DefaultAssay(obj) <- "RNA"
  obj@project.name <- project

  ## QC filtering of cells
  message("performing QC filtering")
  obj$percent.mito <- Seurat::PercentageFeatureSet(obj, pattern = "^mt-|^Mt-")
  obj$percent.ribo <- Seurat::PercentageFeatureSet(obj, pattern = "^RP[LS]|^Rp[ls]")
  summary(obj$nFeature_RNA)
  obj <- subset(obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2000 & percent.mito < 20)

  ## Equalizing libraries
  if ("library_id" %in% colnames(obj@meta.data)) {
    ncells <- median(table(obj$library_id))
    message("library_id parameter found. Equalizing cells to ", ncells)
    sel <- unlist(tapply(1:ncol(obj), obj$library_id, head, ncells))
  }

  if ("batch" %in% colnames(obj@meta.data)) {
    message("batch parameter found. integrating batches using MNN.")
    split <- Seurat::SplitObject(obj, split.by = "batch")
    for (i in 1:length(split)) {
      split[[i]] <- Seurat::NormalizeData(split[[i]])
      split[[i]] <- Seurat::FindVariableFeatures(split[[i]])
    }
    anchors <- Seurat::FindIntegrationAnchors(
      split,
      dims = 1:30, verbose = FALSE
    )
    genes <- rownames(counts)
    integrated <- Seurat::IntegrateData(
      anchorset = anchors, features.to.integrate = genes,
      dims = 1:30, verbose = FALSE
    )
    obj <- integrated
    Seurat::DefaultAssay(obj) <- "integrated"
  } else {
    obj <- Seurat::NormalizeData(obj)
    obj <- Seurat::FindVariableFeatures(obj)
    Seurat::DefaultAssay(obj) <- "RNA"
  }

  if (ncol(obj) > max.cells) {
    message("Subsampling cells to ", max.cells)
    obj <- subset(obj, cells = sample(Seurat::Cells(obj), max.cells))
  }

  ## Dimensionality reductions
  message("Calculating dimensionality reductions...")
  obj <- Seurat::FindVariableFeatures(obj)

  obj <- Seurat::ScaleData(obj, vars.to.regress = c("percent.mito"))
  obj <- Seurat::RunPCA(obj, npcs = 30, verbose = FALSE)
  obj <- Seurat::RunUMAP(obj, reduction = "pca", dims = 1:30, dim.embed = 2)
  obj <- Seurat::RunTSNE(obj, reduction = "pca", dims = 1:30, dim.embed = 2)
  obj <- Seurat::FindNeighbors(obj, reduction = "pca", dims = 1:30)
  obj <- Seurat::FindClusters(obj)

  ## Pre-calculate markers

  message("Finished!")
  obj
}
#' Differential expression analysis for two factorial designs
#'
#' @title Differential expression analysis using limma for two factorial designs
#'
#' @description Performs differential expression analysis using limma for data with two factorial designs.
#' Designed for experiments with two factors, each with two or more levels.
#'
#' @param X Numeric matrix of gene expression values (genes in rows, samples in columns).
#' @param factors Data frame with factor columns defining the experimental design.
#' @param fdr FDR threshold for identifying differentially expressed genes.
#' @param lfc Log2 fold change cutoff for differential expression.
#' @param trend Logical indicating whether to fit a mean-variance trend.
#' @param ref Character vector of reference levels to use as baseline in contrasts.
#' @param compute.means Logical indicating whether to append mean expression values.
#'
#' @details This function handles the limma model matrix and contrasts for two factorial designs.
#' The \code{factors} data frame specifies the experimental design.
#' Multiple contrasts are performed between factor level combinations.
#' Results are filtered by FDR and fold change thresholds.
#'
#' @return Data frame containing limma analysis results for each contrast.
#'
#' @export
gx.limma.two.factorial <- function(X, factors, fdr = 0.05, lfc = 0.20, trend = FALSE,
                                   ref = c(
                                     "ctrl", "ctr", "control", "dmso", "nt", "0", "0h", "0hr",
                                     "non", "no", "not", "neg", "negative", "ref", "veh", "vehicle",
                                     "wt", "wildtype", "untreated", "normal", "false", "healthy"
                                   ),
                                   compute.means = TRUE) {
  ## LIMMA

  cat("Two-factorial LIMMA\n")
  cat("analyzing", ncol(X), "samples\n")
  X <- X[!(rownames(X) %in% c(NA, "", "NA")), ]

  ## create factors
  fct <- vector("list", ncol(factors))
  names(fct) <- colnames(factors)
  j <- 1
  for (j in 1:ncol(factors)) {
    vv <- sort(unique(as.character(factors[, j])))
    vv
    if (sum(vv %in% ref) > 0) {
      vv0 <- vv[which(vv %in% ref)]
      vv <- c(sort(vv0), setdiff(vv, vv0))
    }
    fct[[j]] <- factor(factors[, j], levels = vv)
  }

  ## check
  a <- fct[[1]]
  b <- fct[[2]]
  v1 <- levels(a)
  v2 <- levels(b)
  cat("factors 1:", paste(v1), "\n")
  cat("factors 2:", paste(v2), "\n")
  if (!(length(fct) == 2 && length(v1) == 2 && length(v2) == 2)) {
    cat("gx.limma2:: fatal error: only 2-factorial with 2 levels implemented\n")
    return(NULL)
  }

  ## setup LIMMA
  design <- model.matrix(~ a * b)

  colnames(design)[2] <- paste(rev(v1), collapse = "vs")
  colnames(design)[3] <- paste(rev(v2), collapse = "vs")
  colnames(design)[4] <- paste(colnames(factors)[1], colnames(factors)[2], sep = "*")

  ## perform fitting
  fit0 <- limma::lmFit(X, design)
  cc0 <- paste(v1[1], ".", paste(rev(levels(fct[[2]])), collapse = "vs"), sep = "")
  cc1 <- paste(v1[2], ".", paste(rev(levels(fct[[2]])), collapse = "vs"), sep = "")
  cont.matrix <- cbind(c(0, 0, 1, 0), c(0, 0, 1, 1), diff = c(0, 0, 0, 1))
  colnames(cont.matrix) <- c(cc0, cc1, "Diff")
  rownames(cont.matrix) <- colnames(fit0$coefficients)
  fit1 <- limma::contrasts.fit(fit0, cont.matrix)
  fit2 <- limma::eBayes(fit1, trend = trend)

  ## extract toptable
  top1 <- limma::topTable(fit2, coef = colnames(cont.matrix)[1], number = nrow(X))
  top2 <- limma::topTable(fit2, coef = colnames(cont.matrix)[2], number = nrow(X))
  top3 <- limma::topTable(fit2, coef = colnames(cont.matrix)[3], number = nrow(X))
  top1 <- top1[rownames(X), ]
  top2 <- top2[rownames(X), ]
  top3 <- top3[rownames(X), ]

  ## only significant
  kk <- rownames(X)
  sig <- limma::decideTests(fit2, p.value = fdr, lfc = lfc)
  limma::vennDiagram(sig, cex = 0.8)
  title(sub = paste("fdr=", fdr, sep = ""))
  if (fdr < 1) {
    sig.up <- which(sig[, 1] > 0 & sig[, 2] > 0 & sig[, 3] == 0)
    sig.down <- which(sig[, 1] < 0 & sig[, 2] < 0 & sig[, 3] == 0)
    length(sig.up)
    length(sig.down)
    kk <- rownames(sig)[c(sig.up, sig.down)]
    top1 <- top1[kk, ]
    top2 <- top2[kk, ]
    top3 <- top3[kk, ]
    sig <- sig[kk, ]
  }
  cat("found", nrow(sig), "significant at fdr=", fdr, "and minimal FC=", lfc, "\n")

  ## compute means if requested
  X.m <- NULL
  if (compute.means && nrow(sig) > 0) {
    ff <- paste(factors[, 1], factors[, 2], sep = ".")
    X.m <- t(apply(X[rownames(sig), ], 1, function(x) tapply(x, ff, mean)))
  }

  ## fold-change
  logfc <- cbind(top1$logFC, top2$logFC, top3$logFC)
  colnames(logfc) <- colnames(cont.matrix)
  rownames(logfc) <- rownames(top1)
  jj <- order(-abs(rowMeans(logfc[, 1:2])))
  top1 <- top1[jj, ]
  top2 <- top2[jj, ]
  top3 <- top3[jj, ]
  sig <- sig[jj, ]
  logfc <- logfc[jj, ]
  if (!is.null(X.m)) X.m <- X.m[jj, ]

  ## pq-values
  pv <- cbind(top1$P.Value, top2$P.Value, top3$P.Value)
  qv <- cbind(top1$adj.P.Val, top2$adj.P.Val, top3$adj.P.Val)
  tt <- cbind(top1$t, top2$t, top3$t)
  rownames(pv) <- rownames(qv) <- rownames(top1)
  colnames(pv) <- colnames(qv) <- colnames(sig)

  ## simple summary
  ss0 <- data.frame(
    logFC = rowMeans(logfc[, 1:2]),
    t = rowMeans(tt[, 1:2]),
    P.Value = rowMeans(pv[, 1:2]),
    adj.P.Val = rowMeans(qv[, 1:2]),
    AveExpr = X.m
  )

  ## results
  res <- c()
  res$fdr <- fdr
  res$means <- X.m
  res$limma <- list(top1, top2, top3)
  names(res$limma) <- colnames(cont.matrix)
  res$sig <- sig
  res$logFC <- logfc
  res$p.value <- pv
  res$q.value <- qv
  res$summary <- ss0

  return(res)
}
#' @title Create a hive plot visualization for an NGS object
#'
#' @description
#' Creates a hive plot to visualize associations between omics datasets
#' in an NGS object.
#'
#' @param ngs An NGS object containing multiple omics datasets.
#' @param pheno The phenotype data column name.
#' @param level The omics data level (1, 2, 3).
#' @param ntop Number of top features to include.
#' @param main Plot title.
#' @param axis.lab Axis labels c("GX", "CN", "ME").
#' @param bkgnd Background color.
#' @param rx Hive plot x-radius.
#' @param cex Text size factor.
#' @param slen Feature name substring length.
#'
#' @details
#' This function takes an NGS object containing gene expression (GX), copy number (CN),
#' and methylation (ME) data. It calculates correlation statistics between features in
#' each dataset vs. the phenotype. It selects the top \code{ntop} features and creates
#' a hive plot visualization, with hubs representing GX/CN/ME and edges representing
#' correlations between features across datasets related to the phenotype.
#'
#' The plot is customizable via parameters like \code{main}, \code{bkgnd}, \code{cex}, etc.
#'
#' @return
#' A hive plot grob object is returned, no value.
#'
#' @export
ngs.hiveplot <- function(ngs, pheno, level = 1, ntop = 400, main = "", axis.lab = c("GX", "CN", "ME"),
                         bkgnd = "white", rx = 2.2, cex = 1, slen = -1) {
  res <- omx$zstats[[level]]
  res <- cbind(res, omx$stats[[pheno]][[level]])
  res <- res[!grepl("^PHENO", rownames(res)), ]
  rownames(res) <- gsub("^.*:", "", rownames(res)) ## strip prefix
  if (slen > 0) rownames(res) <- substr(rownames(res), 1, slen)
  ann.cex <- cex * 0.7

  hpd <- omx.makeHivePlotData_(res,
    rho.min = 0.15, ntop = ntop, rx = rx,
    cxi = (0.05 + ann.cex / 7)
  )
  write.csv(hpd$nodes.ann, file = "/tmp/annode.csv", row.names = FALSE)
  grid::grid.newpage()
  axlab.col <- ifelse(bkgnd == "black", "grey90", "grey15")
  HiveR::plotHive(hpd,
    np = FALSE, ch = 1.4, bkgnd = bkgnd,
    axLabs = c(paste0(main, "\n", axis.lab[1]), axis.lab[2], axis.lab[3]),
    axLab.pos = c(0.4, 0.33, 0.33),
    axLab.gpar = grid::gpar(col = axlab.col, cex = cex * 1.3),
    anNodes = "/tmp/annode.csv",
    anNode.gpar = grid::gpar(cex = 0.7 * cex, col = axlab.col, lwd = 0.50)
  )
}
#' Make contrasts between clusters
#'
#' @param clusters A vector of cluster assignments
#' @param min.freq The minimum cluster frequency required to include in contrasts
#' @param full Logical indicating whether to include all pairwise contrasts between clusters
#' @param by.sample Logical indicating whether to cross clusters with samples to get sample x cluster contrasts
#'
#' @return A contrast matrix
#'
#' @description
#' Creates a contrast matrix to compare clusters against each other.
#'
#' @details
#' This function takes a vector of cluster assignments and generates contrasts between the clusters.
#'
#' It first creates indicator columns for each cluster.
#' Then it subtracts the diagonal to generate comparisons of each cluster versus the others.
#' Clusters below min.freq are excluded.
#'
#' If full=TRUE, all pairwise contrasts between clusters are also included.
#'
#' If by.sample=TRUE, the cluster indicators are crossed with the samples
#' to generate sample x cluster interaction contrasts.
#'
#' @export
makeClusterContrasts <- function(clusters, min.freq = 0.01, full = FALSE,
                                 by.sample = FALSE) {
  idx <- sort(unique(as.character(clusters)))
  m1 <- model.matrix(~ 0 + idx)
  colnames(m1) <- sub("^idx", "", colnames(m1))
  rownames(m1) <- colnames(m1)
  colnames(m1) <- paste0(colnames(m1), "_vs_others")
  m1 <- t(t(m1 == 1) / Matrix::colSums(m1 == 1) - t(m1 == 0) / Matrix::colSums(m1 == 0))

  diag(m1) <- 1
  if (full == TRUE) {
    m2 <- makeFullContrasts(unique(clusters))
    m2 <- m2[rownames(m1), ]
    m1 <- cbind(m1, m2)
  }
  if (by.sample) {
    design <- model.matrix(~ 0 + clusters)
    colnames(design) <- sub("^clusters", "", colnames(design))
    rownames(design) <- names(clusters)
    design <- design[, rownames(m1)]
    m1 <- design %*% m1
  }
  return(m1)
}
#' @describeIn read.csv3 debug BAK
#' @export
dbg.BAK <- function(...) {
  if (exists("DEBUG") && DEBUG) {
    msg <- sapply(list(...), paste, collapse = " ")
    message(cat(paste0("[DBG] ", sub("\n$", "", paste(msg, collapse = " ")), "\n")))
  }
}
#' Get gene sets from playbase data
#'
#' @param pattern Pattern to match gene set names
#' @param lib.dir Directory containing custom gene set files
#' @param custom_families_file Custom gene families file name. Default "custom-families.gmt".
#'
#' @return A list of gene sets matching the pattern
#'
#' @details This function extracts gene sets from the playbase data package.
#' It returns gene sets matching the provided pattern.
#' Custom additional gene sets can be included from the lib.dir if provided.
#'
#' @export
getGSETS_playbase.SAVE <- function(pattern, lib.dir, custom_families_file = "custom-families.gmt") {
  # Note: this function needs to be refactored since lib.dir will not exist
  # get gene symbols
  GENE.SYMBOL <- unlist(as.list(org.Hs.eg.db::org.Hs.egSYMBOL))
  # get f1 and families
  FAMILIES <- pgx.getGeneFamilies(GENE.SYMBOL, min.size = 10, max.size = 9999)
  fam.file <- file.path(lib.dir, custom_families_file)
  if (file.exists(fam.file)) {
    custom.gmt <- read.gmt(file.path(lib.dir, custom_families_file), add.source = TRUE)
    names(custom.gmt)
    FAMILIES <- c(FAMILIES, custom.gmt)
  }
  FAMILIES[["<all>"]] <- GENE.SYMBOL
  f1 <- FAMILIES
  names(f1) <- paste0("FAMILY:", names(f1))
  names(f1) <- sub("FAMILY:<all>", "<all>", names(f1))

  # get GSETS
  GSETS <- c(playdata::GSETS, f1)

  # to get iGSETS
  GSET.GENES <- sort(unique(unlist(GSETS))) ## slow...
  iGSETS <- parallel::mclapply(GSETS, function(a) match(a, GSET.GENES)) ## slow...
  gs <- grep(pattern, names(iGSETS), value = TRUE)
  names(iGSETS) <- names(GSETS)

  lapply(iGSETS[gs], function(i) GSET.GENES[i])
}
#' @export
cutGraph0 <- function(g, k = 3, cut = FALSE) {
  ## Cluster graph and cut crossing edges if requested.
  clust <- igraph::cluster_louvain(g)
  if (cut) {
    g <- igraph::delete_edges(g, igraph::E(g)[crossing(clust, g)])
  }
  cmp <- clust$membership
  cmp.size <- table(cmp)
  cmp.top <- Matrix::head(names(cmp.size)[order(-cmp.size)], k)
  igraph::induced_subgraph(g, which(cmp %in% cmp.top))
}
#' @export
pgx.readDatasetProfiles <- function(pgx.dir, file = "datasets-allFC.csv", verbose = TRUE) {
  if (!dir.exists(pgx.dir)) {
    stop(paste("[readDatasetProfiles1] FATAL ERROR : folder", pgx.dir, "does not exist"))
  }
  fn <- file.path(pgx.dir, file)
  fn
  if (!file.exists(fn)) {
    stop("FATAL : could not find profiles matrix. please create first with initDatasetFolder().\n")

    return()
  } else {
    if (verbose) message("[readDatasetProfiles1] Found existing dataset profiles matrix")
  }

  allFC <- fread.csv(file = file.path(pgx.dir, file), row.names = 1, check.names = FALSE)
  allFC <- as.matrix(allFC)
  if (verbose) message("[readDatasetProfiles1] dataset profiles matrix : dim=", dim(allFC))
  return(allFC)
}
#' Clean names for GSEA analysis
#'
#' This function cleans the names of a vector or list for use in GSEA (Gene Set Enrichment Analysis) analysis.
#' It replaces special characters and consecutive underscores in the names with a single underscore.
#'
#' @param s A vector or list with names to be cleaned.
#'
#' @export
#'
#' @return A character vector or list with cleaned names.
#'
gsea.clean_names <- function(s) {
  gsub("__|___", "_", gsub("[-+: *?!$,'|\\.]", "_", names(s)))
}
#' Write data to a GCT file
#'
#' This function writes data to a GCT (Gene Cluster Text) file format.
#' The GCT format is commonly used to store gene expression data in a tabular format.
#'
#' @param X The data matrix or data frame to be written to the GCT file.
#'          Rows represent genes and columns represent samples.
#' @param file The file path to save the GCT file.
#'
#' @export
#'
#' @return None
#'
write.gct <- function(X, file) {
  gX <- data.frame(NAME = rownames(X), DESCRIPTION = NA, X)
  write("#1.2", file = file)
  suppressWarnings(write(paste(nrow(X), ncol(X), sep = "\t"), file = file, append = TRUE))
  suppressWarnings(write.table(format(gX[, ], digits = 3),
    file = file, append = TRUE,
    row.names = FALSE, quote = FALSE, sep = "\t"
  ))
}
#' @title Prepare count data for edgeR differential analysis
#'
#' @description This function takes raw count data and prepares it for differential
#' expression analysis with edgeR. It handles filtering, normalization, and batch correction.
#'
#' @param counts Data frame with raw counts, or named list with counts, samples, and genes components.
#' @param samples Data frame with sample metadata, required if counts is just a matrix.
#' @param genes Data frame with gene metadata, required if counts is just a matrix.
#' @param normalization Normalization method to use, default is "none".
#' @param filter Logical indicating whether to filter low counts.
#' @param prior.cpm Filtering threshold for counts per million. Default is 0.
#' @param remove.batch Logical indicating whether to estimate and remove batch effects.
#'
#' @return A DGEList object ready for analysis with edgeR.
#'
#' @details This function handles conversion of raw count data into a properly formatted
#' DGEList object for differential analysis with edgeR. It checks that the sample and gene
#' metadata matches the counts. Filtering, normalization, and batch correction can be applied.
#'
#' The returned DGEList contains the normalized, filtered count matrix, along with sample and
#' gene data frames, ready for analysis with edgeR functions like glmQLFit() and glmQLFTest().
#'
#' @export
ngs.cookForEDGER <- function(counts, samples = NULL, genes = NULL, normalization = "none",
                             filter = TRUE, prior.cpm = 0, remove.batch = TRUE) {
  if (all(c("counts", "samples", "genes") %in% names(counts))) {
    samples <- counts$samples
    genes <- counts$genes
    counts <- counts$counts
  }
  #
  if (!all(colnames(counts) == rownames(samples))) stop("samples do not match")
  if (!all(rownames(counts) == rownames(genes))) stop("genes do not match")
  if (is.null(samples)) stop("need samples specified")
  if (is.null(genes)) stop("need genes specified")

  ## ------------------------------------------------------------
  ## Now create an DGEList object (see also tximport Vignette)
  ## ------------------------------------------------------------


  cooked <- edgeR::DGEList(round(counts), group = NULL) ## we like integer counts...
  cooked$samples$group <- NULL
  cooked$samples <- cbind(cooked$samples, samples)
  if (!is.null(genes)) cooked$genes <- genes

  ## filter out non-expressed genes (using edgeR standard function)
  if (filter) {
    keep <- edgeR::filterByExpr(cooked)
    table(keep)
    cooked <- cooked[keep, ]
  }

  ## normalized for RNA composition (TMM)
  cooked <- edgeR::calcNormFactors(cooked, method = "TMM")

  ## ------------------------------------------------------------------
  ## prior count regularization
  ## ------------------------------------------------------------------
  ## based on the above, we add a prior count
  if (prior.cpm > 0) {
    cat("adding prior counts at prior.cpm=", prior.cpm, "\n")
    CPM.FACTORS <- Matrix::colSums(counts) / 1e6
    prior.counts <- (prior.cpm * CPM.FACTORS)
    summary(prior.counts)
    cooked$counts <- t(t(cooked$counts) + prior.counts)
  }

  if (normalization != "none") {
    ## remove batch-effects with LIMMA. Be sure to include batch in the
    ## model to avoid overfitting.
    has.batch <- ("batch" %in% colnames(cooked$samples))
    if (has.batch && remove.batch == TRUE) {
      cat("Found 'batch' column in sample table. Correcting for batch effects...\n")
      #
      design1 <- model.matrix(~group, data = cooked$samples)
      batch1 <- cooked$samples[, "batch"]
      xnorm <- limma::removeBatchEffect(xnorm, batch = batch1, design = design1)
    } else {
      ##
    }

    ## quantile normalize. be careful, may introduce artifacts (e.g. clipping)
    if (normalization == "quantile") {
      xnorm <- limma::normalizeQuantiles(xnorm)
    }

    ## we need to undo log and normalizations for further analysis???
    cooked$counts <- 2**(xnorm) ## undo log2
    cooked$samples$norm.factors <- 1 ## we assumed its already used
    cooked$samples$lib.size <- round(Matrix::colSums(cooked$counts)) ## update lib.size
  }

  return(cooked)
}
#' @export
getGseaTable <- function(path) {
  ff <- dir(path, full.names = TRUE)
  report_name <- ff[grep("gsea_report.txt$", ff)]
  if (is.null(report_name) || length(report_name) == 0) {
    return(NULL)
  }
  R <- read.csv(report_name, sep = "\t", check.names = FALSE)
  R <- R[order(-abs(R$NES)), ]
  colnames(R)[1] <- "NAME"
  R$NES <- round(R$NES, digits = 3)
  kk <- grep("NAME|SIZE|NES|NOM|FDR|LEADING", colnames(R), ignore.case = TRUE)
  R <- R[, kk]
  return(R)
}
#' @export
pgx.computeMultiOmicsGSE <- function(X, gmt, omx.type,
                                     method = NULL, center = TRUE) {
  if (is.null(omx.type)) {
    omx.type <- gsub("[:=].*", "", rownames(X))
  }
  omx.types <- setdiff(unique(omx.type), c("MIR", ""))
  omx.types

  sx <- list()
  for (tp in omx.types) {
    x1 <- X[which(omx.type == tp), ]
    rownames(x1) <- sub(":.*", "", rownames(x1))
    sx[[tp]] <- pgx.computeGeneSetExpression(x1, gmt, method = method, center = center)
    sx[[tp]] <- lapply(sx[[tp]], function(x) {
      rownames(x) <- paste0(tp, "=", rownames(x))
      x
    })
  }

  ## concatenate all omx-types
  cx <- sx[[1]]
  for (j in 1:length(sx[[1]])) {
    cx[[j]] <- do.call(rbind, lapply(sx, "[[", j))
  }

  return(cx)
}
#' @export
silac.fitWeibull <- function(obj, protein, samples = NULL, samples2 = NULL) {
  if (!is.null(samples) && inherits(samples[1], "character") &&
    samples[1] %in% names(obj$groups)) {
    samples <- obj$groups[[samples]]
  }
  if (!is.null(samples2) && inherits(samples2[1], "character") &&
    samples2[1] %in% names(obj$groups)) {
    samples2 <- obj$groups[[samples2]]
  }
  if (is.null(samples)) samples <- rownames(obj$samples)
  timeN <- obj$samples[samples, "SILAC"]
  timeN <- as.integer(sub("h$", "", timeN))
  timeN


  ## Nice examples are TCF7, ENO1, GAPDH, FOXO1, B2M, CD5, SQSTM1, CD3E, CXCR4, FOXP1, RPS9 as examples

  tag <- rownames(obj$proteins)[which(obj$proteins$Gene.names == protein)]
  if (length(tag) > 1) cat("warning:: protein has multiple tags!\n")
  tag <- tag[1] ## really???
  q1 <- obj$LFQ.ratio[tag, samples]

  plot(timeN, q1,
    pch = 17, col = "blue",
    ylim = c(0, 1), ylab = "", xlab = "Time [h]", main = protein,
    cex = 1.4, xaxt = "n", yaxt = "n", cex.lab = 1.5, cex.main = 1.4
  )
  axis(1, at = c(0, 6, 12, 24, 48), las = 0, cex.axis = 1.5)
  axis(2, at = c(0, 0.25, 0.5, 0.75, 1), las = 2, cex.axis = 1.5)

  abline(h = 0.5, lty = 3)
  abline(h = 0.25, lty = 3)
  abline(h = 0.75, lty = 3)
  fit1 <- fit.weibull2(timeN, q1)
  lines(fit1$x, fit1$y, lty = 1, col = "blue")
  tt <- paste0("t50= ", round(fit1$t50, digits = 2), "h")

  fit2 <- NULL
  if (!is.null(samples2)) {
    timeM <- obj$samples[samples2, "SILAC"]
    timeM <- as.integer(sub("h$", "", timeM))
    timeM

    q2 <- obj$LFQ.ratio[tag, samples2]
    points(timeM, q2, pch = 19, col = "#EF4136")
    fit2 <- fit.weibull2(timeM, obj$LFQ.ratio[tag, samples2])
    lines(fit2$x, fit2$y, lty = 1, col = "#EF4136")
    tt <- c(tt, paste0("t50= ", round(fit2$t50, digits = 2), "h"))
  }

  legend("bottomright", legend = tt, pch = c(17, 19), col = c("blue", "#EF4136"), cex = 1, bty = "n")

  res <- list(fit = fit1, fit2 = fit2)
  return(res)
}
#' @title Normalize count data
#'
#' @description
#' Normalizes a matrix of RNA-seq count data using different methods.
#'
#' @param counts Matrix of count data, with genes in rows and samples in columns.
#' @param method Normalization method to use. Options are "TMM", "RLE", "upperquartile", or "none".
#' @param log Logical, if TRUE apply log2 transformation after normalization.
#' @param prior Count offset to add prior to normalization. Default 1.
#'
#' @details
#' This function normalizes a matrix of RNA-seq count data using different normalization methods:
#'
#' - "TMM": Trimmed Mean of M-values method from edgeR
#' - "RLE": Relative Log Expression method
#' - "upperquartile": Upper quartile normalization
#' - "none": No normalization
#'
#' Counts are offset by prior value before normalization. A log2 transformation can also be applied after normalization.
#'
#' @return
#' The normalized count matrix.
#'
#' @examples
#' \dontrun{
#' counts <- matrix(rnbinom(10000, mu = 10, size = 1), 100, 100)
#' normalized <- pgx.countNormalization(counts, method = "TMM", log = TRUE)
#' }
#' @export
pgx.removeBatchEffect <- function(X, batch, model.vars = NULL,
                                  method = c("ComBat", "BMC", "limma", "MNN", "fastMNN")) {
  ## treat as factor variable
  method <- method[1]
  batch0 <- as.character(batch)
  batch0[is.na(batch0)] <- "NA" ## NA as separate group??
  if (method == "MNN") {
    matlist <- tapply(1:ncol(X), batch0, function(i) X[, i, drop = FALSE])

    suppressWarnings(out <- do.call(
      scran::mnnCorrect,
      c(matlist, pc.approx = TRUE)
    ))
    new.X <- do.call(cbind, out$corrected)
    colnames(new.X) <- unlist(lapply(matlist, colnames))
    X <- new.X[, colnames(X)]
  } else if (method == "fastMNN") {
    d <- min(50, ncol(X) / 2)
    matlist <- tapply(1:ncol(X), batch0, function(i) X[, i, drop = FALSE])
    out <- do.call(fastMNN, c(matlist, d = d))
    cor.exp <- tcrossprod(out$rotation[, ], out$corrected)
    rownames(cor.exp) <- rownames(X)
    colnames(cor.exp) <- unlist(lapply(matlist, colnames))
    X <- cor.exp[, colnames(X)]
  } else if (method == "limma") {
    X <- limma::removeBatchEffect(X, batch = batch0)
  } else if (method == "ComBat") {
    X <- sva::ComBat(X, batch = batch0)
  } else if (method == "BMC") {
    ## batch mean center
    matlist <- tapply(1:ncol(X), batch0, function(i) X[, i, drop = FALSE])
    matlist <- lapply(matlist, function(x) (x - Matrix::rowMeans(x, na.rm = TRUE)))
    new.X <- do.call(cbind, matlist)
    new.X <- new.X + Matrix::rowMeans(X, na.rm = TRUE)
    X <- new.X[, colnames(X)]
  } else {
    dbg("ERROR! uknown method\n")
  }
  return(X)
}
#' @title Plot omics network from PGX object
#'
#' @param pgx PGX object containing network
#' @param gene Gene name to highlight
#' @param reduced Use reduced network
#' @param levels Data levels to include
#' @param contrast Contrast to color by
#' @param layout Network layout algorithm
#' @param colorcluster Color by module
#' @param hilight Additional genes to highlight
#'
#' @return A network plot object
#'
#' @description
#' Visualize the omics network from a PGX analysis
#'
#' @details
#' This function extracts the omics network stored in a PGX object.
#' It generates a network graph, coloring and sizing nodes by the
#' specified \code{contrast}.
#'
#' The \code{reduced} network can be plotted for a simplified view.
#' \code{levels} determines which data levels are included.
#'
#' \code{gene} and \code{hilight} highlight specific genes.
#' \code{colorcluster} colors nodes by module.
#' \code{layout} specifies the network layout algorithm.
#'
#' @export
pgx.plotOmicsNetwork <- function(pgx, gene = NULL, reduced = NULL, levels = c("gene", "geneset"),
                                 contrast = NULL, layout = NULL, colorcluster = FALSE,
                                 hilight = NULL) {
  ## for col2hex


  has.graph <- all(c("omicsnet", "omicsnet.reduced") %in% names(pgx))
  has.graph
  if (!has.graph) {
    return(NULL)
  }


  gr <- pgx$omicsnet
  if (is.null(gr)) {
    return(NULL)
  }
  if (is.null(reduced) && is.null(gene)) gr <- pgx$omicsnet.reduced

  if (!is.null(gene)) {
    gene0 <- paste0("{gene}", gene)
    k <- which(igraph::V(gr)$name %in% c(gene, gene0))
    nb <- names(igraph::neighbors(gr, igraph::V(gr)[k]))
    vv <- unique(c(gene0, nb))
    gr <- igraph::induced_subgraph(gr, vv)
  }
  gr <- igraph::induced_subgraph(gr, which(igraph::V(gr)$level %in% levels))
  if (is.null(gr)) {
    return(NULL)
  }

  ## ------------- get fold-change for node color and size ------------------
  fc0 <- gr$foldchange[igraph::V(gr)$name, , drop = FALSE]
  if (is.null(contrast)) {
    fc <- rowMeans(fc0**2, na.rm = TRUE)**0.5
  } else {
    if (!(contrast %in% colnames(fc0))) stop("unknown contrast")
    fc <- fc0[, contrast]
  }
  fc[is.na(fc)] <- 0
  fc <- fc / max(abs(fc), na.rm = TRUE)
  fc.cex <- (0.01 + abs(fc))**0.66

  ## defaults graph parameters
  vlabel <- igraph::V(gr)$name
  if ("label" %in% igraph::vertex_attr_names(gr)) vlabel <- igraph::V(gr)$label
  vlabel0 <- vlabel
  vklr <- c("blue3", "grey70", "red3")[2 + sign(fc)]

  if (colorcluster) {
    vklr <- rep(rainbow(16), 99)[igraph::V(gr)$cluster]
  }
  lab.cex <- 1

  ee <- igraph::get.edgelist(gr)
  ee <- igraph::get.edges(gr, igraph::E(gr))

  ew <- 1 + 5 * sqrt(fc.cex[ee[, 1]] * fc.cex[ee[, 2]])
  ew <- 1 + 5 * abs(igraph::E(gr)$weight)
  vsel <- rep(0, length(fc))
  esel <- rep(0, nrow(ee))

  ## ------------------ highlight selection with labels
  if (!is.null(hilight) && length(hilight)) {
    sel <- hilight

    mm <- igraph::V(gr)$name
    if (!is.null(gr$members)) {
      mm <- gr$members[igraph::V(gr)$name]
    }
    mm <- sapply(mm, function(s) sub(".*\\}", "", s))
    vlabel <- sapply(mm, function(x) intersect(x, sel))
    vlabel <- sapply(vlabel, paste, collapse = "\n")

    sel <- which(vlabel != "")
    if (length(sel) > 0) {
      vsel[sel] <- 1
      lab.cex[sel] <- 1 + 18 * (fc.cex[sel] / max(fc.cex[sel], na.rm = TRUE))
    }

    jj <- which(vsel[ee[, 1]] == 1 | vsel[ee[, 2]] == 1)
    esel[jj] <- 1
    ew[jj] <- 2.4 * ew[jj]

    nnb <- unique(unlist(sapply(sel, neighbors, graph = gr)))
    is.nb <- (1:length(sel) %in% nnb)
    vlabel[which(vsel == 0)] <- NA

    vklr[which(vsel == 0)] <- "grey60"
    lab.cex[which(vsel == 0)] <- 1
  }


  vklr <- substring(gplots::col2hex(vklr), 1, 7)
  names(vklr) <- igraph::V(gr)$name
  igraph::V(gr)$label <- vlabel ## filtered labels
  igraph::V(gr)$title <- gsub("\n", "<br>", vlabel0) ## tooltip has complete names

  igraph::V(gr)$size <- 40 * fc.cex
  igraph::V(gr)$color <- paste0(vklr, ifelse(vsel == 1, "99", "55"))


  igraph::E(gr)$color <- paste0(vklr[ee[, 1]], ifelse(esel == 1, "99", "55"))
  igraph::E(gr)$width <- 1 * (2 + 5 * (ew / max(ew)))

  if (!is.null(layout)) {
    layout.fun <- match.fun(layout)
    tmp.gr <- gr
    igraph::E(tmp.gr)$weight <- abs(igraph::E(tmp.gr)$weight)**0.2

    pos <- layout.fun(tmp.gr)
    remove(tmp.gr)
    rownames(pos) <- igraph::V(gr)$name
  }

  visdata <- visNetwork::toVisNetworkData(gr, idToLabel = FALSE)
  pos <- pos[igraph::V(gr)$name, ]
  pos[, 2] <- -pos[, 2]

  ## ------------------ plot using visNetwork (zoomable) -----------------
  graph <- visNetwork::visNetwork(
    nodes = visdata$nodes, edges = visdata$edges,
    height = "1200px", width = "1600px"
  ) %>%
    visNetwork::visNodes(font = list(size = 14)) %>%
    visNetwork::visEdges(hidden = FALSE, width = 2, color = list(opacity = 0.9)) %>%
    visNetwork::visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T)) %>%
    ## visHierarchicalLayout(direction = "LR") %>%
    ## visInteraction(hideEdgesOnDrag = TRUE) %>%
    visNetwork::visIgraphLayout(layout = "layout.norm", layoutMatrix = pos)
  graph
}
#' @export
silac.plotProteins <- function(obj, proteins, samples = NULL) {
  if (!is.null(samples) && inherits(samples[1], "character") &&
    samples[1] %in% names(obj$groups)) {
    samples <- obj$groups[[samples]]
  }

  ss <- match(proteins, obj$proteins$Gene.names)
  ss <- setdiff(ss, NA)
  pp <- rownames(obj$proteins)[ss]
  Q <- obj$LFQ.ratio
  if (is.null(samples)) samples <- colnames(Q)
  plot.df <- data.frame(Q[pp, samples])
  row.names(plot.df) <- proteins
  PL <- t(plot.df)
  boxplot(t(plot.df), ylim = c(0, 1), col = "darkgrey", las = 3)
  stripchart(list(PL[, 1], PL[, 2], PL[, 3], PL[, 4]), vertical = T, add = T, pch = 1, method = "jitter", cex = 1.5)
}
#' @title Build PMID annotation matrix
#'
#' @description
#' Builds an annotation matrix mapping PubMed IDs to gene symbols
#' based on mappings in org.Hs.eg.db.
#'
#' @details
#' This function retrieves PubMed ID to Entrez Gene ID mappings from
#' org.Hs.eg.db and converts them to a sparse matrix mapping PMIDs to
#' gene symbols. It filters to PMIDs associated with <=10 genes.
#'
#' It collapses duplicate PMID mappings and builds a sparse matrix
#' with rows as PMIDs, columns as gene symbols, and entries indicating
#' an association between a PMID and symbol.
#'
#' @return
#' A sparse matrix with rows as PMIDs, columns as gene symbols, and
#' entries indicating a mapping between a PMID and symbol.
#'
#' @export
pmid.buildMatrix <- function() {
  pmid <- as.list(org.Hs.eg.db::org.Hs.egPMID2EG)
  symbol <- as.list(org.Hs.eg.db::org.Hs.egSYMBOL)
  eg <- names(symbol)
  symbol <- sapply(symbol, "[", 1)
  names(symbol) <- as.character(eg)
  ngene <- sapply(pmid, length)
  pmid <- pmid[which(ngene <= 10)]

  ## collapse duplicates
  pmid.gg <- sapply(pmid, paste, collapse = ",")
  idx <- tapply(names(pmid), pmid.gg, function(x) x)
  pmid <- pmid[sapply(idx, "[", 1)]
  idx <- lapply(idx, function(x) paste0("PMID:", x))
  names(pmid) <- sapply(idx, paste, collapse = ",")

  ## build PMID2SYMBOL matrix

  idx0 <- parallel::mclapply(1:length(pmid), function(i) cbind(i, which(eg %in% pmid[[i]])),
    mc.cores = NCORE()
  )
  idx <- do.call(rbind, idx0)
  P <- Matrix::sparseMatrix(
    i = idx[, 1], j = idx[, 2], x = rep(1, nrow(idx)),
    dims = c(length(idx0), length(eg))
  )
  rownames(P) <- names(pmid)
  colnames(P) <- symbol
  P <- P[, which(Matrix::colSums(P) > 0)]
  return(P)
}
#' @export
getGseaOutputDir <- function(path) {
  ## untangle gsea subfolder
  gsea_dir <- dir(path)[grep("\\.Gsea\\.", dir(path))]
  gsea_dir <- file.path(path, gsea_dir)
  gsea_dir
  if (length(gsea_dir) == 0) {
    cat("Missing results")
    return(NULL)
  }
  return(gsea_dir)
}
#' GSEA Barplot
#'
#' Generates a barplot visualization of GSEA scores.
#'
#' @param scores A numeric vector of GSEA scores.
#' @param names A character vector of names for the scores.
#' @param xlab The label for the x-axis.
#' @param xlim The limits for the x-axis.
#' @param cex.text The size of the text labels.
#' @param main The main title of the plot.
#' @param n The number of scores to display in the barplot.
#' @export
#' @return barplot
gsea.barplot <- function(scores, names = NULL, xlab = "score", xlim = NULL,
                         cex.text = 1, main = "enrichment", n = 16) {
  if (!is.null(names)) names(scores) <- names
  scores <- rev(Matrix::head(scores[order(-abs(scores))], n))
  col1 <- c("lightskyblue1", "rosybrown1")[1 + 1 * (scores > 0)]
  if (min(scores, na.rm = TRUE) == 0 && max(scores, na.rm = TRUE) == 0) {
    xlim <- c(0, 1)
  }
  barplot(abs(scores),
    horiz = TRUE, las = 1, width = 5 / 6, col = col1,
    border = NA, xlab = xlab, names.arg = rep("", length(scores)),
    xlim = xlim, cex.axis = 0.9
  )
  title(main)
  gs <- names(scores)
  mx <- max(abs(scores))
  text(0.02 * mx, 1:length(scores) - 0.45, gs,
    adj = 0, offset = 0,
    cex = cex.text
  )
}
#' Delete pgx entry in datasets-info table in files (WIP)
#'
#' @export
pgxinfo.deletePgx <- function(pgx.dir, pgxname,
                              delete.fc = FALSE) {
  info.file <- file.path(pgx.dir, "datasets-info.csv")
  pgxname <- sub("[.]pgx$", "", pgxname)
  pgxinfo <- read.csv(info.file, row.names = 1)
  pgxinfo

  info_datasets <- sub("[.]pgx$", "", pgxinfo$dataset)
  sel <- which(info_datasets == pgxname)
  sel
  if (!is.null(sel) && length(sel) > 0) {
    pgxinfo <- pgxinfo[-sel, , drop = FALSE]
    write.csv(pgxinfo, file = info.file)
  }

  ## Should we also delete the entry in allFC and sigdb? This will
  ## take some overhead, so if not needed better skip.
  allfc.file <- file.path(pgx.dir, "datasets-allFC.csv")
  tsne.file <- file.path(pgx.dir, "datasets-tsne.csv")
  h5.file <- file.path(pgx.dir, "datasets-sigdb.h5")

  ## delete columns from allFC file
  if (delete.fc && file.exists(allfc.file)) {
    allFC <- data.table::fread(allfc.file, check.names = FALSE) ## HEADER!!!
    allFC <- as.matrix(allFC, rownames = 1)
    del <- grep(paste0("\\[", pgxname, "\\]"), colnames(allFC))
    del
    if (length(del) > 0) {
      allFC <- allFC[, -del, drop = FALSE]
      allFC <- round(allFC, digits = 4)
      allFC1 <- data.frame(gene = rownames(allFC), allFC, check.names = FALSE)
      data.table::fwrite(allFC1, allfc.file) ## HEADER!!!
    }
  }

  ## delete rows from t-SNE file
  if (delete.fc && file.exists(tsne.file)) {
    tsne <- data.table::fread(tsne.file, check.names = FALSE) ## HEADER!!!
    tsne <- as.matrix(tsne, rownames = 1)
    del <- grep(paste0("\\[", pgxname, "\\]"), rownames(tsne))
    if (length(del)) {
      tsne <- tsne[-del, , drop = FALSE]
      write.csv(tsne, file = tsne.file)
    }
  }

  ## delete dataset from H5 file
  if (delete.fc && file.exists(h5.file)) {
    sigdb.removeDataset(h5.file, pgxname)
  }
}
#' @title Compute enrichment for drug combinations
#'
#' @param obj PGX object containing gene expression data
#' @param X Gene expression matrix
#' @param xdrugs Character vector of drug names
#' @param ntop Number of top drugs to combine per sample
#' @param nsample Number of samples to sample combinations
#' @param nprune Number of genes to use for GSEA
#' @param contrasts Contrasts to compute enrichment
#' @param res.mono Single drug enrichment results to use
#'
#' @return List with enrichment results for each drug combination
#'
#' @description Computes enrichment of drug combinations using single drug enrichment results.
#'
#' @details This function takes single drug enrichment results and computes enrichment for
#' combinations by combining the individual drug profiles. It samples combinations from the
#' top drugs per sample then runs GSEA to test for enrichment.
#'
#' The output is a list containing the enrichment results for each sampled drug combination.
#' This allows searching for synergistic combinations with higher enrichment than individual drugs.
#'
#' @export
pgx.computeComboEnrichment <- function(obj, X, xdrugs,
                                       ntop = 10, nsample = 20, nprune = 250,
                                       contrasts = NULL, res.mono = NULL) {
  if ("gx.meta" %in% names(obj)) {
    F <- sapply(obj$gx.meta$meta, function(x) x$meta.fx)
    rownames(F) <- rownames(obj$gx.meta$meta[[1]])
    ## check if multi-omics
    is.multiomics <- any(grepl("\\[gx\\]|\\[mrna\\]", rownames(F)))
    is.multiomics
    if (is.multiomics) {
      jj <- grep("\\[gx\\]|\\[mrna\\]", rownames(F))
      F <- F[jj, , drop = FALSE]
    }
    rownames(F) <- toupper(sub(".*:|.*\\]", "", rownames(F)))
    F <- F[order(-rowMeans(F**2)), , drop = FALSE]
    F <- F[!duplicated(rownames(F)), , drop = FALSE]
  }

  if (is.null(contrasts)) contrasts <- colnames(F)
  contrasts <- intersect(contrasts, colnames(F))

  ## calculate average drug profile
  if (is.null(res.mono)) {
    cat("Calculating single drug enrichment using GSEA ...\n")
    er.mono <- pgx.computeDrugEnrichment(
      obj, X, xdrugs,
      methods = "GSEA",
      nprune = nprune, contrast = NULL
    )
    er.mono <- er.mono[["GSEA"]]
  } else {
    cat("Using passed single drug enrichment results...\n")
    if ("GSEA" %in% names(res.mono)) {
      er.mono <- res.mono[["GSEA"]]
    }
  }

  ## determine top-combinations

  top.mono.up <- apply(er.mono$X, 2, function(x) Matrix::head(order(-x), ntop))
  top.mono.dn <- apply(er.mono$X, 2, function(x) Matrix::head(order(x), ntop))
  top.combo.up <- apply(top.mono.up, 2, function(idx) list(combn(idx, 2)))
  top.combo.dn <- apply(top.mono.dn, 2, function(idx) list(combn(idx, 2)))
  top.combo.up <- unlist(top.combo.up, recursive = FALSE)
  top.combo.dn <- unlist(top.combo.dn, recursive = FALSE)
  top.combo.up <- lapply(top.combo.up, function(d) apply(d, 2, function(j) paste(sort(j), collapse = "-")))
  top.combo.dn <- lapply(top.combo.dn, function(d) apply(d, 2, function(j) paste(sort(j), collapse = "-")))
  top.combo <- unique(c(unlist(top.combo.up), unlist(top.combo.dn)))

  ## -------------- sample pairs from original mono-matrix
  sample.pairs <- list()
  k <- 1
  for (k in 1:length(top.combo)) {
    cmbn.idx <- as.integer(strsplit(top.combo[k], split = "-")[[1]])
    cmbn.idx
    cmbn <- sort(rownames(er.mono$X)[cmbn.idx])
    p1 <- sample(which(xdrugs == cmbn[1]), nsample, replace = TRUE)
    p2 <- sample(which(xdrugs == cmbn[2]), nsample, replace = TRUE)
    pp <- cbind(p1, p2)
    sample.pairs[[k]] <- pp
  }
  sample.pairs <- do.call(rbind, sample.pairs)

  ## --------------- now create combination matrix X
  comboX <- apply(sample.pairs, 1, function(ii) rowMeans(X[, ii], na.rm = TRUE))

  combo.drugs <- apply(sample.pairs, 1, function(ii) paste(sort(xdrugs[ii]), collapse = "+"))

  Matrix::tail(sort(table(combo.drugs)))
  sum(table(combo.drugs) >= 15)
  sel.combo <- names(which(table(combo.drugs) >= 15))
  jj <- which(combo.drugs %in% sel.combo)
  comboX <- comboX[, jj, drop = FALSE]
  combo.drugs <- combo.drugs[jj]
  colnames(comboX) <- paste0(combo.drugs, "_combo", 1:ncol(comboX))

  cat("Calculating drug-combo enrichment using GSEA ...\n")
  res.combo <- pgx.computeDrugEnrichment(
    obj,
    X = comboX, xdrugs = combo.drugs, methods = "GSEA", nprune = nprune,
    contrast = NULL
  )
  res.combo <- res.combo[["GSEA"]]

  return(res.combo)
}
#' Plot mitochondrial and ribosomal counts
#'
#' @description This function plots the mitochondrial and ribosomal counts in a given count matrix.
#'
#' @param counts A numeric matrix of count data, where rows represent features and columns represent samples.
#' @param percentage Logical, whether to plot the counts as percentages of total counts.
#'
#' @return A bar plot of mitochondrial and ribosomal counts.
#'
#' @export
pgx.plotMitoRibo <- function(counts, percentage = TRUE) {
  tot.counts <- Matrix::colSums(counts, na.rm = TRUE)
  sel.mt <- grep("^mt-", rownames(counts), ignore.case = TRUE)
  sel.rb <- grep("^rp[ls]", rownames(counts), ignore.case = TRUE)
  mito.counts <- Matrix::colSums(counts[sel.mt, , drop = FALSE], na.rm = TRUE)
  ribo.counts <- Matrix::colSums(counts[sel.rb, , drop = FALSE], na.rm = TRUE)
  other.counts <- tot.counts - mito.counts - ribo.counts

  df <- cbind(ribo = ribo.counts, mito = mito.counts)
  if (percentage) df <- round((df / tot.counts) * 100, digits = 2)
  barplot(t(df), beside = FALSE, las = 3)
}
#' @export
prot.plotVolcano <- function(res, use.q = TRUE, psig = 0.05, lfc = 1, pmin = 1e-12,
                             pcex = 0.7, cex = 0.7, plot.type = "default", ...) {
  pdata <- cbind(res[, c("logFC", "P.Value")], res$gene)
  colnames(pdata) <- c("EFFECTSIZE", "P", "Gene")
  if (use.q) pdata$P <- res$adj.P.Val
  sig <- (pdata$P < psig & abs(pdata$EFFECTSIZE) > lfc)
  highlight <- pdata$Gene[which(sig)]
  pdata$P <- pmin + pdata$P

  if (plot.type == "volcanoly") {
    volcanoly(pdata[, ],
      snp = "Gene", highlight = highlight,
      effect_size_line = c(-1, 1) * lfc,
      genomewideline = -log10(pig)
    )
  } else {
    klr <- c("black", "red")[1 + 1 * sig]
    ylab <- "significance (-log10p)"
    if (use.q) ylab <- "significance (-log10q)"
    plot(pdata$EFFECTSIZE, -log10(pdata$P),
      xlab = "effect size (log2FC)", ylab = ylab, ...,
      pch = 20, cex = pcex, col = klr
    )
    abline(v = c(-1, 1) * lfc, lty = 2, col = "grey70")
    abline(h = -log10(psig), lty = 2, col = "grey70")
    jj <- which(sig)
    if (length(jj)) {
      text(pdata$EFFECTSIZE[jj], -log10(pdata$P)[jj],
        pdata$Gene[jj],
        cex = cex, adj = 0.5,
        offset = 0.4, pos = 3
      )
    }
  }
}
#' Create stratified contrasts design matrix
#'
#' @param data Data frame containing variables and strata column
#' @param vars Character vector of variable names
#' @param strata Name of stratification variable column
#' @param ref Reference level for strata
#'
#' @return A design matrix for stratified contrasts
#'
#' @description
#' Generates a design matrix for stratified contrasts from a data frame.
#'
#' @details
#' This function takes a data frame containing variables of interest and a stratification
#' variable. It constructs a model matrix encoding contrasts between a reference level
#' and other levels of the stratification variable. The contrasts are computed separately
#' within each stratum.
#'
#' The result is a design matrix with columns for the variables of interest, followed by
#' columns encoding the stratified contrasts. This matrix can be used as a design matrix in
#' differential expression analysis.
#'
#' @export
pgx.makeStratifiedContrastsDF <- function(data, vars, strata, ref) {
  dstrata <- data[, strata]
  S <- model.matrix(~ 0 + dstrata)
  colnames(S) <- sub("^dstrata", "", colnames(S))

  df <- cbind(data[, vars, drop = FALSE], strata = dstrata)
  md <- makeDirectContrasts(df, c(ref, "others"))

  M <- md$contr.matrix
  C0 <- M[, grep("strata", colnames(M), invert = TRUE), drop = FALSE]
  C1 <- M[, grep("strata", colnames(M)), drop = FALSE]

  colnames(C1) <- paste0(gsub("strata:|_vs_others", "", colnames(C1)))
  contr.matrix <- c()
  i <- 1
  for (i in 1:ncol(C0)) {
    m1 <- (C1 == 1) * C0[, i]
    colnames(m1) <- paste0(colnames(C0)[i], "@", colnames(C1))
    contr.matrix <- cbind(contr.matrix, m1)
  }

  ## construct exp matrix
  grp <- md$group
  G <- model.matrix(~ 0 + grp)
  colnames(G) <- sub("^grp", "", colnames(G))
  G <- G[, rownames(contr.matrix)]
  exp.matrix <- G %*% contr.matrix
  rownames(exp.matrix) <- rownames(data)

  ## check levels
  sel <- (Matrix::colSums(contr.matrix == -1) > 0 &
    Matrix::colSums(contr.matrix == +1) > 0)
  contr.matrix <- contr.matrix[, sel]
  exp.matrix <- exp.matrix[, sel]

  res <- list(
    contr.matrix = contr.matrix,
    group = md$group,
    exp.matrix = exp.matrix
  )
  return(res)
}
#' @export
merge_similarities <- function(S.list, sparse = NULL) {
  length(S.list)
  if (is.null(sparse)) sparse <- is(S.list[[1]], "sparseMatrix")
  if (length(S.list) > 1 && sparse) {
    ## merge similarities (i.e. multiply) but allow for missing values
    S.list <- lapply(S.list, function(m) Matrix::Matrix(m, sparse = TRUE))
    cat("merge similarities...(sparse)\n")
    idx <- Matrix::which(S.list[[1]] != 0, arr.ind = TRUE)
    x <- S.list[[1]][idx]
    p <- rep(1, length(x))
    i <- 2
    for (i in 2:length(S.list)) {
      not.na <- !is.na(S.list[[i]][, 1])
      jj <- which(not.na[idx[, 1]] & not.na[idx[, 2]])
      x[jj] <- x[jj] * S.list[[i]][idx[jj, ]]
      p[jj] <- p[jj] + 1
    }
    table(p)
    x <- x**(1 / p)
    table(x != 0)
    jj <- which(x != 0)
    n <- nrow(S.list[[1]])
    Q <- Matrix::sparseMatrix(
      i = idx[jj, 1], j = idx[jj, 2], x = x[jj],
      use.last.ij = TRUE, dims = c(n, n)
    )
    ## set NA??
    jj <- which(Matrix::rowMeans(Q == 0) == 1)
    if (length(jj) > 0) {
      Q[jj, 1] <- NA
      Q[1, jj] <- NA
    }
  } else if (length(vars) > 1 && !sparse) {
    ## merge similarities (i.e. multiply) but allow for missing values
    cat("merge similarities... (full)\n")
    S.list <- lapply(S.list, function(m) Matrix::Matrix(m, sparse = FALSE))
    Q <- as.matrix(S.list[[1]])
    p <- matrix(1, nrow(Q), ncol(Q))
    p <- 1 * (!is.na(Q))
    i <- 2
    for (i in 2:length(S.list)) {
      s1 <- as.matrix(S.list[[i]])
      not.na <- !is.na(S.list[[i]][, 1])
      j0 <- which(!not.na)
      s1[j0, ] <- 1
      s1[, j0] <- 1
      Q <- Q * s1
      j1 <- which(not.na)
      p[j1, j1] <- p[j1, j1] + 1
    }
    Q <- Q**(1 / p)
  } else {
    Q <- S.list[[1]]
  }
  rownames(Q) <- colnames(Q) <- rownames(S.list[[1]])
  return(Q)
}
#' Detect batch parameters in expression matrix
#'
#' @title Detect batch parameters in expression matrix
#'
#' @description Detects potential batch effect parameters in the column names of
#' an expression matrix.
#'
#' @param Y Expression matrix with samples in columns.
#'
#' @details This function checks the column names of the expression matrix Y for
#' common batch effect terms like 'batch', 'cell.line', 'patient', etc. It returns
#' the column names that match these terms, which typically indicate a batch effect parameter.
#'
#' @return A character vector of detected batch parameter names.
#'
#' @export
pgx.detect_batch_params <- function(Y) {
  grep("batch|cell.line|patient|mouse|sample|repl|strain", colnames(Y), value = TRUE)
}
#' Get Gene Set Tables
#'
#' This function retrieves gene set tables from a specified directory path and organizes the data into a list.
#'
#' @param path The path to the directory containing the gene set tables.
#'
#' @return A list containing the gene set tables and associated metadata.
#' @export
getGeneSetTables <- function(path) {
  dd <- dir(path, full.names = TRUE)
  dd
  tables <- vector("list", length(dd))
  d <- dd[1]
  i <- 1
  for (i in 1:length(dd)) {
    d <- dd[i]
    table_dir <- file.path(d, "tables")
    if (!file.exists(table_dir)) next
    ff <- dir(table_dir, full.names = TRUE)
    ff.short <- dir(table_dir, full.names = FALSE)
    ff.name <- gsub("output-|-results.csv|-test.csv", "", ff.short)
    tables[[i]] <- lapply(ff, read.csv, row.names = 1, check.names = FALSE)
    names(tables[[i]]) <- ff.name
    names(tables)[i] <- gsub(".*/", "", d)
  }

  gsets <- sort(unique(unlist(lapply(tables[[1]], rownames))))
  length(gsets)

  meta <- list()
  i <- 1
  for (i in 1:length(dd)) {
    pv <- c()
    fx <- c()
    j <- 1
    for (j in 1:length(tables[[i]])) {
      x <- tables[[i]][[j]]
      pv.col <- grep("p.value|^p$|p-val|pval|nom.p.val|nom p-val|p.val", tolower(colnames(x)))[1]
      pv <- cbind(pv, x[match(gsets, rownames(x)), pv.col])
      fx.col <- grep("sign|nes|logfc|fc", tolower(colnames(x)))[1]
      fx <- cbind(fx, (x[match(gsets, rownames(x)), fx.col]))
      colnames(pv)[ncol(pv)] <- names(tables[[i]])[j]
      colnames(fx)[ncol(fx)] <- names(tables[[i]])[j]
    }
    Matrix::head(pv)
    rownames(pv) <- gsets
    rownames(fx) <- gsets
    meta[[i]] <- data.frame(geneset = gsets, fx = fx, p = pv)
  }
  names(meta) <- names(tables)
  p

  res <- list(tables = tables, meta = meta)
  return(res)
}
#' @export
matrix.prod <- function(..., na.value = 1) {
  matlist <- list(...)
  if (inherits(matlist[[1]], "list")) matlist <- matlist[[1]]
  p <- Reduce("+", lapply(matlist, function(x) 1 * !is.na(x)))
  i <- 1
  for (i in 1:length(matlist)) {
    if (any(is.na(matlist[[i]]))) {
      jj <- Matrix::which(is.na(matlist[[i]]), arr.ind = TRUE)
      matlist[[i]][jj] <- na.value
    }
  }
  xprod <- Reduce("*", matlist)
  xprod <- xprod**(1 / p)
  jj <- Matrix::which(p == 0, arr.ind = TRUE)
  if (length(jj) > 0) {
    xprod[jj] <- NA
  }
  return(xprod)
}
#' @export
pgx.correlateSignature.matrix <- function(fc, refmat, nsig = 100, ntop = 1000, nperm = 10000) {
  ##
  ##
  ##
  ##

  if (is.null(names(fc))) stop("fc must have names")

  ## mouse... mouse...
  names(fc) <- toupper(names(fc))

  ## or instead compute correlation on top100 fc genes (read from file)

  rn <- rownames(refmat)
  cn <- colnames(refmat)

  ## ---------------------------------------------------------------
  ## Compute simple correlation between query profile and signatures
  ## ---------------------------------------------------------------
  gg <- intersect(rn, names(fc))
  fc1 <- sort(fc[gg])
  gg <- unique(names(c(Matrix::head(fc1, nsig), Matrix::tail(fc1, nsig))))


  G <- refmat[gg, , drop = FALSE]

  ## rank correlation??
  rG <- apply(G[gg, ], 2, rank, na.last = "keep")
  rfc <- rank(fc[gg], na.last = "keep")


  rG[is.na(rG)] <- 0
  rfc[is.na(rfc)] <- 0
  rho <- stats::cor(rG, rfc, use = "pairwise")[, 1]

  remove(G, rG, rfc)

  ## --------------------------------------------------
  ## test all signature on query profile using fGSEA
  ## --------------------------------------------------


  sel <- Matrix::head(names(sort(-abs(rho))), ntop)
  notx <- setdiff(sel, colnames(refmat))
  sel <- intersect(sel, colnames(refmat))
  X <- refmat[, sel, drop = FALSE]
  X[is.na(X)] <- 0
  orderx <- apply(X, 2, function(x) {
    idx <- order(x)
    list(DN = head(idx, 100), UP = rev(Matrix::tail(idx, 100)))
  })
  sig100.dn <- sapply(orderx, "[[", "DN")
  sig100.dn <- apply(sig100.dn, 2, function(i) rn[i])
  sig100.up <- sapply(orderx, "[[", "UP")
  sig100.up <- apply(sig100.up, 2, function(i) rn[i])

  ## ---------------------------------------------------------------
  ## combine up/down into one (unsigned GSEA test)
  ## ---------------------------------------------------------------
  gmt <- rbind(sig100.up, sig100.dn)
  gmt <- unlist(apply(gmt, 2, list), recursive = FALSE)
  names(gmt) <- colnames(X)

  suppressMessages(suppressWarnings(
    res <- fgsea::fgseaSimple(gmt, abs(fc), nperm = nperm)
  ))

  ## ---------------------------------------------------------------
  ## Combine correlation+GSEA by combined score (NES*rho)
  ## ---------------------------------------------------------------
  jj <- match(res$pathway, names(rho))
  res$rho <- rho[jj]
  res$R2 <- rho[jj]**2
  res$score <- res$R2 * res$NES
  res <- res[order(res$score, decreasing = TRUE), ]

  return(res)
}
#' @export
jaccard_similarity <- function(m) {
  A <- tcrossprod(m)
  im <- which(A > 0, arr.ind = TRUE, useNames = F)
  b <- rowSums(m)
  Aim <- A[im]
  x <- Aim / (b[im[, 1]] + b[im[, 2]] - Aim)
  sp <- Matrix::sparseMatrix(i = im[, 1], j = im[, 2], x = x, dims = dim(A))
  rownames(sp) <- colnames(sp) <- rownames(m)
  sp
}
#' @export
graph_from_pos.DEPRECATED <- function(pos, trh) {
  if (ncol(pos) > 3 || NCOL(pos) == 1) {
    stop("positions must be 2 or 3 columns\n")
  }
  d <- apply(pos, 1, function(x) Matrix::colSums((t(pos) - x)**2)) ## huge??
  w <- 1 / (1 + d / mean(d))
  g <- igraph::graph_from_adjacency_matrix(w, mode = "undirected", weighted = TRUE, diag = FALSE)
  g <- igraph::subgraph.edges(g, which(igraph::E(g)$weight > trh))
  g$layout <- pos
  return(g)
}
#' @export
gseaSnapshot <- function(gsets, gsea_dir) {
  enplots <- dir(gsea_dir, pattern = "enplot_")
  enplots0 <- dir(gsea_dir, pattern = "enplot_", full.names = TRUE)
  kk <- match(gsets, gsub("enplot_|_[0-9]*.png$", "", enplots))
  imgs <- lapply(enplots0[kk], function(p) {
    grid::rasterGrob(as.raster(png::readPNG(p)),
      interpolate = TRUE
    )
  })
  gridExtra::grid.arrange(grobs = imgs, ncol = 5)
}
#' @export
clustering.score <- function(xy, labels, kinter = NULL) {
  ## Defined as the ratio between inter-distance between clusters
  ## and intra-distance within the clusters.
  ##
  nlab <- length(unique(labels))
  groups <- tapply(1:nrow(xy), labels, list)
  groups

  xy <- scale(xy)
  xy[is.na(xy)] <- 0
  xy.mean <- t(sapply(groups, function(g) colMeans(xy[g, , drop = FALSE])))
  inter.dist <- pos2dist(xy.mean)
  intra.dist <- sapply(groups, function(g) mean(pos2dist(xy[g, , drop = FALSE])))
  intra.dist

  diag(inter.dist) <- NA

  if (!is.null(kinter)) {
    d1 <- inter.dist * NA
    for (i in 1:nrow(inter.dist)) {
      jj <- Matrix::head(order(inter.dist[i, ]), kinter)
      d1[i, jj] <- d1[jj, i] <- inter.dist[i, jj]
    }
    d1
    inter.dist <- d1
  }
  inter.dist
  intra.dist

  intra.dist0 <- intra.dist + 0.1 * mean(intra.dist, na.rm = TRUE)
  S <- outer(intra.dist0, intra.dist0, FUN = "+")**0.5
  inter.t <- inter.dist / S
  inter.t
  score <- mean(inter.t, na.rm = TRUE)
  return(score)
}
#' Plot expression of two genes across cell types
#'
#' @description
#' Generates a scatter plot showing the expression levels of two genes across inferred cell types.
#'
#' @param pgx A PGX object containing single-cell expression data and inferred cell types.
#' @param gene1 First gene name or ID to plot.
#' @param gene2 Second gene name or ID to plot.
#' @param cex Point size scaling factor. Default is 1.
#' @param col Point color. Default is "grey60".
#' @param lab.unit Units for axis labels, eg 'log2(TPM)'. Default is NULL.
#' @param cex.names Text size for sample labels. Default is 1.
#' @param samples Vector of sample names to highlight. Default is NULL.
#' @param k Number of contour levels if drawing density contours. Default is 11.
#'
#' @details
#' This function takes a PGX object containing single-cell expression data and inferred cell types.
#' It extracts the expression values for the two specified genes, and generates a scatter plot
#' with each cell as a point colored by its annotated cell type. Contour lines indicate the density distribution.
#' Sample names can be highlighted, and axis labels customized.
#'
#' @return
#' A scatter plot is generated showing the expression distribution of the two genes across cell types.
#'
#' @export
pgx.cytoPlot <- function(pgx, gene1, gene2, cex = 1, col = "grey60",
                         lab.unit = NULL, cex.names = 1, samples = NULL, k = 11) {
  ## some pretty colors

  my.cols <- rev(RColorBrewer::brewer.pal(k, "RdYlBu"))


  if (is.null(samples)) {
    samples <- colnames(pgx$X)
  }
  samples <- intersect(samples, colnames(pgx$X))
  x1 <- pgx$X[gene1, samples]
  x2 <- pgx$X[gene2, samples]
  x1 <- x1 + 1e-3 * rnorm(length(x1))
  x2 <- x2 + 1e-3 * rnorm(length(x2))
  names(x1) <- samples
  names(x2) <- samples
  m1 <- mean(x1)
  m2 <- mean(x2)

  ## select samples in different quadrants
  j1 <- samples[which(x1 < m1 & x2 > m2)]
  j2 <- samples[which(x1 > m1 & x2 < m2)]
  j3 <- samples[which(x1 > m1 & x2 > m2)]
  j4 <- samples[which(x1 < m1 & x2 < m2)]

  z1 <- z2 <- z3 <- z4 <- NULL
  if (length(j1) > 1) z1 <- MASS::kde2d(x1[j1], x2[j1], n = 50)
  if (length(j2) > 1) z2 <- MASS::kde2d(x1[j2], x2[j2], n = 50)
  if (length(j3) > 1) z3 <- MASS::kde2d(x1[j3], x2[j3], n = 50)
  if (length(j4) > 1) z4 <- MASS::kde2d(x1[j4], x2[j4], n = 50)



  xlab1 <- paste(gene1, lab.unit, collapse = "  ")
  ylab1 <- paste(gene2, lab.unit, collapse = "  ")
  plot(x1, x2, xlab = xlab1, ylab = ylab1, col = col, pch = 19, cex = cex)
  abline(h = mean(x1), v = mean(x2), lwd = 1, lty = 2)



  if (length(j1) > 10) contour(z1, drawlabels = FALSE, nlevels = k, col = my.cols, add = TRUE, lwd = 2)
  if (length(j2) > 10) contour(z2, drawlabels = FALSE, nlevels = k, col = my.cols, add = TRUE, lwd = 2)
  if (length(j3) > 10) contour(z3, drawlabels = FALSE, nlevels = k, col = my.cols, add = TRUE, lwd = 2)
  if (length(j4) > 10) contour(z4, drawlabels = FALSE, nlevels = k, col = my.cols, add = TRUE, lwd = 2)

  N <- length(x1)
  d1 <- 0.02 * max(x1)
  d2 <- 0.04 * max(x2)
  legend("topright", paste(round(100 * length(j3) / N, 2), "%"),
    cex = 1.2, col = "gray50",
    bty = "n", xpd = TRUE
  )
  legend("topleft", paste(round(100 * length(j1) / N, 2), "%"),
    cex = 1.2, col = "gray50",
    bty = "n", inset = c(-0.05, 0), xpd = TRUE
  )
  legend("bottomright", paste(round(100 * length(j2) / N, 2), "%"),
    cex = 1.2, col = "gray50",
    bty = "n", xpd = TRUE
  )
  legend("bottomleft", paste(round(100 * length(j4) / N, 2), "%"),
    cex = 1.2, col = "gray50",
    bty = "n", inset = c(-0.05, 0), xpd = TRUE
  )

  if (!is.null(pgx$deconv)) {
    inferred.celltype <- pgx$deconv[[1]][["meta"]]

    lab1 <- Matrix::head(names(sort(-Matrix::colSums(inferred.celltype[j1, , drop = FALSE]))), 3)
    pos1 <- apply(cbind(x1, x2)[j1, , drop = FALSE], 2, median)
    text(pos1[1], pos1[2], paste(lab1, collapse = "\n"), cex = 0.9 * cex.names, xpd = TRUE)

    lab2 <- Matrix::head(names(sort(-Matrix::colSums(inferred.celltype[j2, , drop = FALSE]))), 3)
    pos2 <- apply(cbind(x1, x2)[j2, , drop = FALSE], 2, median)
    text(pos2[1], pos2[2], paste(lab2, collapse = "\n"), cex = 0.9 * cex.names, xpd = TRUE)

    lab3 <- Matrix::head(names(sort(-Matrix::colSums(inferred.celltype[j3, , drop = FALSE]))), 3)
    pos3 <- apply(cbind(x1, x2)[j3, , drop = FALSE], 2, median)
    text(pos3[1], pos3[2], paste(lab3, collapse = "\n"), cex = 0.9 * cex.names, xpd = TRUE)
  }
}
#' @export
pgx.scFilterOutliers <- function(counts, a = 2.5, plot = FALSE) {
  ## QC filter of (single) cells like Seurat
  ## See https://satijalab.org/seurat/v3.2/pbmc3k_tutorial.html

  ## --------------------------------------------------------------
  ## calculate percentages
  ## --------------------------------------------------------------
  mt.genes <- grep("^MT-", rownames(counts), ignore.case = TRUE, value = TRUE)
  rb.genes <- grep("^RP[SL]", rownames(counts), ignore.case = TRUE, value = TRUE)
  percent.mito <- Matrix::colSums(counts[mt.genes, ]) / Matrix::colSums(counts) * 100
  percent.ribo <- Matrix::colSums(counts[rb.genes, ]) / Matrix::colSums(counts) * 100
  nfeature <- Matrix::colSums(counts > 0)
  ncounts <- Matrix::colSums(counts)

  if (plot) {
    log.nfeature <- log10(1 + nfeature)
    log.ncounts <- log10(1 + ncounts)
    log.mito <- log10(1 + percent.mito)
    log.ribo <- log10(1 + percent.ribo)
    nfeature.th <- mean(log.nfeature) + c(-a, a) * sd(log.nfeature)
    ncounts.th <- mean(log.ncounts) + c(-a, a) * sd(log.ncounts)
    mito.th <- mean(log.mito) + c(-a, a) * sd(log.mito)
    ribo.th <- mean(log.ribo) + c(-a, a) * sd(log.ribo)
    mito.th
    ribo.th

    par(mfrow = c(2, 2))
    hist(log.nfeature, breaks = 100)
    abline(v = nfeature.th, col = "red")
    hist(log.ncounts, breaks = 100)
    abline(v = ncounts.th, col = "red")

    hist(log.mito, breaks = 100)
    abline(v = mito.th, col = "red")
    hist(log.ribo, breaks = 100)
    abline(v = ribo.th, col = "red")
  }

  selectInlier <- function(x, a = 2.5) {
    xmin <- mean(x) - a * sd(x)
    xmin <- max(xmin, 0.01 * mean(x))
    xmax <- mean(x) + a * sd(x)
    x > xmin & x < xmax
  }

  selectInlier <- function(x, a = 2.5) {
    x <- log10(1 + x)
    xmin <- mean(x) - a * sd(x)
    xmin <- max(xmin, 0.01 * mean(x))
    xmax <- mean(x) + a * sd(x)
    (x > xmin & x < xmax)
  }


  ## sel <- nfeature < nfeature.th & ncounts < ncounts.th &

  sel <- selectInlier(nfeature, a) &
    selectInlier(ncounts, a) &
    selectInlier(percent.mito, a) &
    selectInlier(percent.ribo, a)

  counts <- counts[, sel]
  counts
}
