# This file is part of the Omics Playground project.
# Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
# Batch-correction diagnostic and comparison plots.
# This file owns visualization of covariates, embeddings, and correction scores.
# Plotting functions consume computed results and do not alter correction policy.

#' @export
bc.plotCovariateHeatmap <- function(bc.res) {
  ## bc <- detectBatchEffects(X, samples, pheno, contrasts = NULL,
  ##                          params = c("statistical", "technical", "pca"),
  ##                          p.pca = 0.5, p.pheno = 0.05,
  ##                          k.pca = 10, nv = 1, xrank = NULL)
  B <- bc.res$covariates_plus
  rho <- cor(apply(B, 2, rank))
  colnames(rho) <- rep("", ncol(rho))
  gx.heatmap(rho,
    sym = TRUE, mar = c(1, 15), keysize = 0.4, cexCol = 0.0001,
    scale = "none", key = FALSE
  )
}

#' @export
bc.plotResults <- function(X, xlist, pos, pheno, samples = NULL, scores = NULL,
                           type = "umap", nmax = 1000, cex = 1, text.cex = 1,
                           ncol = NULL, par = TRUE) {
  ## samples=NULL;scores = NULL;type='umap';nmax=1000;cex=1;text.cex = 1
  if (par) {
    if (is.null(ncol)) {
      ncol <- ceiling(sqrt(length(xlist)))
    }
    nr <- ceiling(length(xlist) / ncol)
    par(mfrow = c(nr, ncol))
  }

  methods <- names(xlist)
  if (!is.null(pos)) {
    methods <- intersect(methods, names(pos))
  }

  if (!is.null(scores)) {
    methods <- intersect(methods, rownames(scores))
    m.score <- scores[methods, "score"]
    methods <- methods[order(-m.score)]
    scores <- scores[methods, ]
    xlist <- xlist[methods]
    pos <- pos[methods]
  }

  cex1 <- ifelse(length(pheno) > 20, 3, 4)
  cex1 <- ifelse(length(pheno) > 100, 2.5, cex1)
  cex1 <- ifelse(length(pheno) > 400, 2, cex1)
  cex1 <- ifelse(length(pheno) > 1000, 1, cex1)
  cex1 <- 0.7 * cex * cex1

  if (tolower(type) == "umap") {
    par(mar = c(2.4, 3, 2, 1))
    i <- 1
    for (m in methods) {
      plot(pos[[m]], col = factor(pheno), pch = 20, cex = cex1)
      title(main = m, cex.main = 1.4, line = 0.5)
      ## tt <- paste("score = ",round(res$results[m,"score"],2))
      ## legend("topleft", legend=tt, cex=1.4)
    }
  }

  if (type == "heatmap") {
    par(mar = c(2, 3, 1.8, 1))
    i <- 1
    for (m in methods) {
      xx <- xlist[[m]]
      xx <- head(xx[order(-apply(xx, 1, sd, na.rm = TRUE)), ], nmax)
      xx <- xx - rowMeans(xx, na.rm = TRUE)
      xx <- abs(xx)**0.5 * sign(xx)
      gx.imagemap(xx, main = m, cex.main = 1.4, cex = 0)
      mtext("samples", 1, line = 0.5, las = 1)
      mtext("genes", 2, line = 0.5, las = 3)
    }
  }

  if (tolower(type) %in% c("pc", "pc2")) {
    if (is.null(samples)) message("samples in NULL!")
    B <- pgx.computeTechnicalEffects(X, nv = 1)
    bcat <- sub("[.].*", "", colnames(B))
    colnames(B) <- paste0(bcat, ":", colnames(B)) ## for collapsing
    if (!is.null(samples)) B <- cbind(B, samples)
    horiz <- ifelse(tolower(type) == "pc2", TRUE, FALSE)

    plist <- list()
    i <- 1
    for (m in methods) {
      xx <- xlist[[m]]
      plist[[m]] <- pgx.PC_correlation(xx, B,
        nv = 3, stat = "F",
        plot = TRUE, main = m, expand = FALSE, collapse = TRUE,
        horiz = horiz, text.cex = text.cex
      )
    }

    gridExtra::grid.arrange(grobs = plist, ncol = ncol, padding = unit(0.0, "line"))
  }

  if (type == "hist") {
    par(mar = c(3, 3, 4, 2))
    i <- 1
    for (m in methods) {
      xx <- xlist[[m]]
      hist(xx, breaks = 100, main = m, cex.main = 1.8)
    }
  }

  if (type == "scores" && !is.null(scores)) {
    m <- 1
    plt <- list()
    ylabs <- c(
      "score" = "overall score",
      "genes" = "significant genes",
      "gsets" = "significant genesets",
      "avg.fc" = "average abs.logFC",
      "avg.sd" = "average SD",
      "r.genes" = "gene.coverage",
      "r.gsets" = "gset coverage",
      "SNR" = "signal-to-noise",
      "pc1.ratio" = "PC1 ratio",
      "silhouette" = "silhoutte score"
    )
    for (i in 1:ncol(scores)) {
      nn <- colnames(scores)[i]
      xx <- scores[, i]
      plt[[i]] <- plot_ggbarplot(
        t(xx),
        ylab = ylabs[nn],
        srt = 45,
        legend.cex = 1.2 * text.cex,
        label.cex = 1.2 * text.cex,
        axis.cex = 1.2 * text.cex,
        group.name = ""
      ) +
        ggplot2::theme(
          plot.margin = ggplot2::margin(t = 0, r = 4, b = 0, l = 8, "pt"),
          plot.title = ggplot2::element_text(size = 15 * text.cex)
        ) +
        ggplot2::xlab("") + ggplot2::ggtitle(nn)
    }

    gridExtra::grid.arrange(grobs = plt, ncol = ncol, padding = unit(0.0, "line"))
  }
}

#' @export
bc.CovariateAnalysisPlot <- function(bc.results, k = 1:3, par = TRUE, col = 1) {
  bc <- bc.results
  if (par == TRUE) par(mfrow = c(2, 2), mar = c(4, 4, 3, 2))
  pp <- intersect(rownames(bc$p.pca), rownames(bc$p.values))
  pp <- grep("^pca", pp, value = TRUE, invert = TRUE)
  k <- k[which(k <= ncol(bc$p.pca))]
  pxx <- bc$p.pca[pp, k, drop = FALSE]
  py <- bc$p.value[pp, 2]
  x1 <- -log10(1e-04 + py)
  col <- rep(col, length(k))
  for (i in c(0, k)) {
    if (i == 0) {
      plot.new()
      abline(h = 0.5, lty = 2)
      abline(v = 0.5, lty = 2)
      mtext("correlation with PC", side = 2, line = 1.3, cex = 0.85)
      mtext("correlation with phenotype", side = 1, line = 1.3, cex = 0.85)
      ##      axis(side=1, tick='n', cex.axis=0.001)
      ##      axis(side=2, tick='n', cex.axis=0.001)
      text(
        x = 0.2, y = 0.80, adj = 0.5,
        labels = "strong batch-effects\nor\nstratification factors"
      )
      text(
        x = 0.75, y = 0.80, adj = 0.5,
        labels = "well designed model-parameters\nor\nstrong confouders"
      )
      text(
        x = 0.2, y = 0.20, adj = 0.5,
        labels = 'nuisance parameters\nor\n"noise"'
      )
      text(
        x = 0.75, y = 0.20, adj = 0.5,
        labels = "weak model-parameters\nor\nweak confouders"
      )
    } else {
      y1 <- -log10(1e-04 + pxx[, i])
      ylim <- c(-0.1 * max(y1, na.rm = TRUE), 1.1 * max(y1, na.rm = TRUE))
      col1 <- col[i]
      plot(x1, y1,
        pch = 20, cex = 1.5, col = col1,
        xlab = "significance with phenotype (-log10p)",
        ylab = "significance with PC (-log10p)",
        xlim = c(-0.4, 4.4), ylim = ylim
      )
      title(paste0("PC", i), cex.main = 1.4)
      text(x1, y1, pp, pos = c(1:4), cex = 1.3, col = col1)
    }
  }
}
