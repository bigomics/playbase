##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

## =========================================================================
## Plotting functions for WGCNA analysis (module heatmaps, dendrograms,
## eigengene graphs, power analysis, etc).
## =========================================================================


#' @seealso WGCNAplus::plotConsensusOverlapHeatmap
#' @export
wgcna.plotConsensusOverlapHeatmap <- function(net1, net2,
                                              setLabels = NULL,
                                              lab.line = c(8, 8),
                                              plotDendro = FALSE,
                                              setpar = TRUE) {
  WGCNAplus::plotConsensusOverlapHeatmap(
    net1 = net1, net2 = net2, setLabels = setLabels,
    lab.line = lab.line, plotDendro = plotDendro, setpar = setpar
  )
}

#' @export
wgcna.plotPreservationSummaries <- function(pres, setpar = TRUE) {
  # Create a simple bar plot of Zsummary:
  Z <- pres$Zsummary
  ntest <- ncol(Z)

  if (setpar) {
    par(mfrow = c(3, ntest), mar = c(5, 5, 4, 1))
  }
  xylist <- list(
    c("moduleSize", "Zsummary.pres"),
    c("moduleSize", "medianRank.pres"),
    c("Zsummary.pres", "medianRank.pres")
  )

  for (xy in xylist) {
    for (k in colnames(Z)) {
      X <- data.frame(
        color = substring(names(pres$moduleSize), 3, 99),
        moduleSize = pres$moduleSize,
        Zsummary.pres = pres$Zsummary[, k],
        medianRank.pres = pres$medianRank[, k]
      )
      xvar <- xy[1]
      yvar <- xy[2]
      ylim <- c(0, max(X[, yvar]))
      if (yvar == "medianRank.pres") ylim <- rev(ylim)
      plot(
        X[, xvar],
        X[, yvar],
        pch = 21,
        cex = 2,
        bg = X$color,
        ylim = ylim,
        xlab = xvar,
        ylab = yvar
      )
      title(yvar, cex.main = 1.4, line = 2.2)
      sub <- paste(k, "vs.", "reference")
      title(sub, cex.main = 1, line = 0.9)
      if (yvar == "Zsummary.pres") abline(h = c(2, 10), lty = 2)
    }
  }
}

#' @export
wgcna.plotPreservationModuleTraits <- function(pres,
                                               subplots = c("zsummary", "consmt", "wt.consmt"),
                                               order.by = "name",
                                               setpar = TRUE, rm.na = FALSE) {
  if (all(is.numeric(subplots))) {
    subplots <- c("zsummary", "consmt", "wt.consmt")[subplots]
  }

  if (setpar) {
    par(mfrow = c(2, 2), mar = c(14, 12, 4, 2))
  }

  ## compute consensus
  Zsummary <- pres$Zsummary

  cR <- pres$modTraits
  ydim <- sapply(pres$layers, function(w) nrow(w$datTraits))
  consZ <- wgcna.computeConsensusMatrix(cR, ydim, psig = 1, consfun = "gmean")
  ## consZ <- consZ[rownames(cR[[1]]), colnames(cR[[1]])]

  ## match
  ii <- intersect(rownames(Zsummary), rownames(consZ))
  Zsummary <- Zsummary[ii, , drop = FALSE]
  consZ <- consZ[ii, , drop = FALSE]

  ## order
  order.method <- "clust"
  if (order.by == "name") {
    ii <- order(rownames(Zsummary))
    Zsummary <- Zsummary[ii, , drop = FALSE]
    consZ <- consZ[ii, , drop = FALSE]
  }
  if (order.by == "zsummary") {
    ii <- order(-rowMeans(Zsummary**2))
    Zsummary <- Zsummary[ii, , drop = FALSE]
    consZ <- consZ[ii, , drop = FALSE]
  }
  if (order.by == "clust") {
    consZ1 <- consZ
    consZ1[is.na(consZ1)] <- 0
    ii <- fastcluster::hclust(dist(consZ1))$order
    jj <- fastcluster::hclust(dist(t(consZ1)))$order
    Zsummary <- Zsummary[ii, , drop = FALSE]
    consZ <- consZ[ii, jj, drop = FALSE]
  }

  ## --------------------------------------
  ## Zsummary heatmap
  ## --------------------------------------
  if ("zsummary" %in% subplots) {
    WGCNA::labeledHeatmap(
      Matrix = Zsummary,
      xLabels = colnames(Zsummary),
      yLabels = rownames(Zsummary),
      ySymbols = rownames(Zsummary),
      colors = tail(WGCNA::blueWhiteRed(100), 50),
      colorLabels = TRUE,
      setStdMargins = FALSE
    )
    title("Module preservation (Zsummary)", line = 1.2, cex.main = 1.2)
  }

  ## --------------------------------------
  ## Consensus Module-Trait
  ## --------------------------------------
  validcol <- function(R) {
    which(colMeans(is.na(R)) < 1 &
      matrixStats::colSds(R, na.rm = TRUE) > 0.01)
  }

  if ("consmt" %in% subplots) {
    clim <- max(abs(consZ), na.rm = TRUE)
    cval <- seq(-clim, clim, length.out = 201)
    ii <- which(cval >= min(consZ, na.rm = TRUE) & cval <= max(consZ, na.rm = TRUE))
    col2 <- WGCNA::blueWhiteRed(201)[ii]
    jj <- 1:ncol(consZ)
    if (rm.na) jj <- validcol(consZ)
    WGCNA::labeledHeatmap(
      Matrix = consZ[, jj, drop = FALSE],
      xLabels = colnames(consZ)[jj],
      yLabels = rownames(consZ),
      ySymbols = rownames(consZ),
      colors = col2,
      colorLabels = TRUE,
      setStdMargins = FALSE
    )
    title("Consensus Module-Traits", line = 1.2, cex.main = 1.2)
  }

  ## --------------------------------------
  ## preservation-weighted Consensus Module-Trait
  ## --------------------------------------
  if ("wt.consmt" %in% subplots) {
    wz <- rowMeans(Zsummary**2, na.rm = TRUE)
    wz <- wz / max(wz)
    consW <- consZ * wz[rownames(consZ)]

    clim <- max(abs(consW), na.rm = TRUE)
    cval <- seq(-clim, clim, length.out = 201)
    ii <- which(cval >= min(consW, na.rm = TRUE) & cval <= max(consW, na.rm = TRUE))
    col2 <- WGCNA::blueWhiteRed(201)[ii]

    jj <- 1:ncol(consW)
    if (rm.na) jj <- validcol(consW)
    WGCNA::labeledHeatmap(
      Matrix = consW[, jj, drop = FALSE],
      xLabels = colnames(consW)[jj],
      yLabels = rownames(consW),
      ySymbols = rownames(consW),
      colors = col2,
      colorLabels = TRUE,
      setStdMargins = FALSE
    )
    title("Preservation-weighted Consensus\nModule-Traits", line = 1, cex.main = 1.2)
  }
}


## =========================================================================
## PLOTTING FUNCTIONS
## =========================================================================


#' @export
wgcna.plotTopModules <- function(wgcna, trait, nmax = 16, setpar = TRUE) {
  MEx <- wgcna$net$MEs
  Y <- wgcna$datTrait
  kk <- intersect(rownames(MEx), rownames(Y))
  MEx <- MEx[kk, ]
  Y <- Y[kk, , drop = FALSE]

  rho <- cor(MEx, Y, use = "pairwise")
  sel <- order(-abs(rho[, trait]))
  sel <- head(sel, nmax)
  n <- length(sel)
  nr <- ceiling(sqrt(n))
  nc <- ceiling(n / nr)

  yclass <- sapply(as.data.frame(Y), class)
  is.binary <- apply(Y, 2, function(x) all(x %in% c(TRUE, FALSE, 0, 1, NA)))
  yclass[which(is.binary)] <- "factor"

  if (setpar == 1) par(mfrow = c(nr, nc), mgp = c(2.6, 0.85, 0), mar = c(4, 4, 2.5, 1))
  if (setpar == 2) par(mfrow = c(nc, nr), mgp = c(2.6, 0.85, 0), mar = c(4, 4, 2.5, 1))

  yclass <- sapply(as.data.frame(Y), class)
  is.binary <- apply(Y, 2, function(x) all(x %in% c(TRUE, FALSE, 0, 1, NA)))
  yclass[which(is.binary)] <- "logical"
  yclass

  i <- sel[1]
  for (i in head(sel, nmax)) {
    x <- Y[, trait]
    y <- MEx[, i]
    label <- colnames(MEx)[i]
    col <- substring(label, 3, 99)
    col1 <- adjustcolor(col, alpha.f = 0.5)

    if (yclass[trait] == "factor") {
      boxplot(y ~ x,
        main = label, col = col1,
        xlab = trait, ylab = "ME score"
      )
      points(1 + x + 0.04 * rnorm(length(x)), y,
        pch = 21, bg = col1, lwd = 0.5
      )
    }

    if (yclass[trait] == "numeric") {
      plot(x, y,
        main = label,
        pch = 21, cex = 1.1, col = 1, bg = col1, lwd = 0.25,
        xlab = trait, ylab = "ME score"
      )
      abline(h = 0, lty = 2, lwd = 0.5)
      r <- cor(x, y, use = "pairwise")
      if (abs(r) > 0.3) {
        abline(lm(y ~ x), col = 1, lwd = 0.6)
        legend("bottomright", legend = paste("r=", round(r, 3)))
      }
    }
  }
}


#' Plot top modules most correlated with trait for multi expression
#' data.
#'
#' @export
wgcna.plotTopModules_multi <- function(multi, trait, nmax = 16, collapse = FALSE,
                                       plotlib = "base", setpar = TRUE) {
  if (!"MEs" %in% names(multi)) {
    multi$MEs <- lapply(multi$net$multiMEs, function(m) m$data)
  }

  ## we 'just' concatenate the ME matrices
  multi$MEs <- lapply(multi$MEs, as.matrix)
  MEx <- as.matrix(mofa.merge_data(multi$MEs))
  Y <- lapply(multi$MEs, function(m) multi$datTraits[rownames(m), , drop = FALSE])
  Y <- mofa.merge_data(Y)

  batch <- sub(":.*", "", rownames(Y))
  names(batch) <- rownames(Y)
  nbatch <- length(multi$MEs)

  ## select top modules
  rho <- cor(MEx, Y, use = "pairwise")
  sel <- names(sort(-abs(rho[, trait])))
  sel <- head(sel, nmax)

  Y <- type.convert(data.frame(Y, check.names = FALSE), as.is = FALSE)
  if (collapse) {
    Y <- collapseTraitMatrix(Y)
    trait <- sub("=.*", "", trait)
  }

  n <- length(sel)
  nr <- ceiling(sqrt(n))
  nc <- ceiling(n / nr)
  if (setpar == 1) par(mfrow = c(nr, nc), mgp = c(2.6, 0.85, 0), mar = c(2.5, 4, 2.5, 1))
  if (setpar == 2) par(mfrow = c(nc, nr), mgp = c(2.6, 0.85, 0), mar = c(2.5, 4, 2.5, 1))

  yclass <- sapply(as.data.frame(Y), class)
  is.binary <- apply(Y, 2, function(x) all(x %in% c(TRUE, FALSE, 0, 1, NA)))
  yclass[which(is.binary)] <- "logical"
  yclass

  i <- sel[1]
  for (i in head(sel, nmax)) {
    x <- Y[, trait]
    y <- MEx[, i]
    label <- i
    col <- substring(label, 3, 99)
    col1 <- adjustcolor(col, alpha.f = 0.66)
    col2 <- col1

    yclass[trait]
    if (yclass[trait] %in% c("factor", "logical", "binary")) {
      if (yclass[trait] %in% c("logical", "binary")) {
        x <- (x == 1)
      }
      df <- data.frame(x = x, y = y, group = factor(batch))

      par(mgp = c(2.4, 0.9, 0))
      aa <- c(0.15, 0.55)
      col1 <- sapply(aa, function(a) adjustcolor(col, alpha.f = a))

      nb <- nbatch
      atx <- c(seq(1, 2, length.out = nb), seq(4, 5, length.out = nb))
      atx <- unlist(sapply(1:nbatch, function(i) i + c(-0.15, 0.15), simplify = FALSE))
      atmid <- 1:nbatch

      boxplot(
        df$y ~ df$x + df$group,
        # df$y ~ df$group + df$x,
        at = atx,
        xlim = c(0.6, nbatch + 0.4),
        cols = col1,
        main = label,
        col = col1,
        boxwex = 0.24,
        xaxt = "n",
        xlab = "",
        ylab = "ME score"
      )

      mtext(levels(df$group), side = 1, line = 0.6, cex = 1.0, at = atmid)
      bb <- c("FALSE", "TRUE")
      legend("topright",
        legend = bb, fill = col1,
        cex = 0.8, y.intersp = 0.82, title = trait, title.cex = 1.1
      )
    } ## end of if factor

    if (yclass[trait] %in% c("numeric", "integer")) {
      df <- data.frame(x = x, y = y, group = factor(batch))
      par(mgp = c(2.4, 0.9, 0))

      aa <- seq(0.55, 0.15, length.out = nbatch)
      col1 <- sapply(aa, function(a) adjustcolor(col, alpha.f = a))
      names(col1) <- sort(unique(batch))

      nb <- nbatch
      colx <- col1[as.integer(factor(batch))]

      plot(
        df$x, df$y,
        pch = 21,
        cex = 1,
        lwd = 0.4,
        bg = colx,
        main = label,
        xlab = trait,
        ylab = "ME score"
      )
      bb <- names(multi$MEs)
      legend("topright",
        legend = bb, fill = col1, cex = 0.9,
        y.intersp = 0.85
      )

      ## add regression lines
      b <- df$group[1]
      col2 <- adjustcolor(col1, red.f = 0.5, green.f = 0.5, blue.f = 0.5)
      names(col2) <- names(col1)
      for (b in unique(df$group)) {
        ii <- which(df$group == b)
        abline(lm(df$y[ii] ~ df$x[ii]), lwd = 1, lty = 1, col = col2[1])
      }
    } ## end of if continuous
  }
}


#' Plot top modules most correlated with trait for multi expression
#' data.
#'
#' @export
wgcna.plotModuleScores <- function(res, trait,
                                   multi = FALSE, nmax = 16,
                                   collapse.trait = FALSE,
                                   plotlib = "base", setpar = TRUE) {
  Y <- res$datTrait

  ## For multi we 'just' concatenate the ME matrices
  if (multi) {
    MEx <- do.call(rbind, lapply(res$MEs, as.matrix))
    me.samples <- lapply(res$MEs, rownames)
    batch <- max.col(sapply(me.samples, function(s) rownames(Y) %in% s))
    batch <- names(res$MEs)[batch]
    names(batch) <- rownames(Y)
    nbatch <- length(res$MEs)
  } else {
    MEx <- res$net$MEs
    batch <- ""
    nbatch <- 1
  }

  ## align
  kk <- intersect(rownames(MEx), rownames(Y))
  MEx <- MEx[kk, ]
  Y <- Y[kk, ]
  if (!is.null(batch)) batch <- batch[kk]

  sel.modules <- colnames(MEx)
  if (nmax > 0) {
    ## select top modules
    rho <- cor(MEx, Y, use = "pairwise")
    sel.modules <- names(sort(-abs(rho[, trait])))
    sel.modules <- head(sel.modules, nmax)
  }
  ncol <- ceiling(sqrt(length(sel.modules)))

  if (collapse.trait) {
    Y <- type.convert(data.frame(Y, check.names = FALSE), as.is = FALSE)
    Y <- collapseTraitMatrix(Y)
    trait <- sub("=.*", "", trait)
  }

  module <- as.vector(sapply(sel.modules, rep, nrow(MEx)))
  dfx <- data.frame(
    sample = rownames(MEx),
    trait = Y[, trait],
    score = as.vector(unlist(MEx[, sel.modules])),
    module = module,
    group = batch
  )

  xtype <- class(type.convert(Y[, trait], as.is = TRUE))
  xtype
  if (xtype != "numeric") {
    dfx$trait <- factor(dfx$trait)
    if (nbatch == 1) {
      ggplot2::ggplot(
        dfx,
        ggplot2::aes(x = factor(trait), y = score, fill = trait)
      ) +
        # ggplot2::aes(y = score, x = trait)) +
        ggplot2::geom_boxplot() +
        ggplot2::xlab(trait) +
        ggplot2::ylab("ME score") +
        ggplot2::facet_wrap(~module, ncol = ncol) +
        ggplot2::theme_bw(base_size = 18)
    }

    if (nbatch > 1) {
      ggplot2::ggplot(
        dfx,
        ggplot2::aes(x = group, y = score, fill = trait)
      ) +
        ggplot2::geom_boxplot() +
        ggplot2::xlab(trait) +
        ggplot2::ylab("ME score") +
        ggplot2::facet_wrap(~module, ncol = ncol) +
        ggplot2::theme_bw(base_size = 18)
    }
  } else {
    dfx$trait <- as.numeric(dfx$trait)
    ggplot2::ggplot(
      dfx,
      # ggplot2::aes(x = trait, y = score, color = group)) +
      ggplot2::aes(x = trait, y = score, color = group)
    ) +
      ggplot2::geom_point(size = 0.6) +
      ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 0.6) +
      ggplot2::xlab(trait) +
      ggplot2::ylab("ME score") +
      ggplot2::facet_wrap(~module, ncol = ncol, scales = "free") +
      ggplot2::theme_bw(base_size = 18)
  }
}

#'
#'
#' @export
wgcna.plotTraitCorrelationBarPlots <- function(res, trait, multi = FALSE,
                                               colored = TRUE, beside = TRUE,
                                               main = NULL, cex.main = 1.3,
                                               setpar = TRUE) {
  if (setpar) {
    nr <- ceiling(sqrt(length(trait)))
    nc <- ceiling(length(trait) / nr)
    par(mfrow = c(nr, nc))
  }
  p <- trait[1]
  for (p in trait) {
    groups <- NULL
    if (multi) {
      mt <- res$modTraits
      groups <- names(mt)
      m1 <- sapply(mt, function(x) x[, p])
    } else {
      m1 <- res$modTraits[, p]
    }
    colnames(m1) <- paste0(p, " (", colnames(m1), ")")
    me.col <- grey.colors(2)
    if (colored) {
      me.col <- sub("ME", "", rownames(m1))
      me.col <- rbind(me.col, me.col)
      aa <- seq(0.7, 0.25, length.out = nrow(me.col))
      for (i in 1:nrow(me.col)) {
        me.col[i, ] <- adjustcolor(me.col[i, ], alpha.f = aa[i])
      }
    }

    if (beside) {
      barplot(t(m1),
        las = 3, beside = TRUE, col = me.col,
        ylab = "trait correlation (rho)"
      )
      tt <- p
      if (!is.null(main)) tt <- main
      title(tt, cex.main = cex.main)
      if (length(groups) > 1) {
        legend("topright", legend = groups, fill = grey.colors(length(groups)))
      }
    } else {
      me.col <- NULL
      if (colored) me.col <- sub("ME", "", rownames(m1))
      for (i in 1:ncol(m1)) {
        barplot(m1[, i],
          las = 3, beside = TRUE, col = me.col,
          ylab = "trait correlation (rho)"
        )
        tt <- colnames(m1)[i]
        if (!is.null(main)) tt <- main
        title(tt, cex.main = cex.main)
      }
    }
  }
}

#' @seealso WGCNAplus::plotTOM
#' @export
wgcna.plotTOM <- function(wgcna, justdata = FALSE, block = NULL,
                          legend = TRUE, downsample = NULL) {
  WGCNAplus::plotTOM(
    wgcna = wgcna, justdata = justdata, block = block,
    legend = legend, downsample = downsample
  )
}

#' @seealso WGCNAplus::plotDendroAndColors
#' @export
wgcna.plotDendroAndColors <- function(wgcna, main = NULL,
                                      extra.colors = NULL,
                                      show.kme = FALSE,
                                      show.traits = FALSE,
                                      show.contrasts = FALSE,
                                      show.tom = FALSE,
                                      show.mat = NULL,
                                      clust = TRUE,
                                      use.tree = 0,
                                      block = 1,
                                      rm.na = TRUE,
                                      sd.wt = 0,
                                      nmax = -1,
                                      marAll = c(0.4, 5, 2, 0),
                                      setLayout = TRUE,
                                      cex = 1,
                                      layout.ncols = 1,
                                      layout.heights = c(0.3, 1, 1),
                                      tom.subsample = 400,
                                      ...) {
  WGCNAplus::plotDendroAndColors(
    wgcna = wgcna, main = main, extra.colors = extra.colors,
    show.kme = show.kme, show.traits = show.traits,
    show.contrasts = show.contrasts, show.tom = show.tom,
    show.mat = show.mat, clust = clust, use.tree = use.tree, block = block,
    rm.na = rm.na, sd.wt = sd.wt, nmax = nmax, marAll = marAll,
    setLayout = setLayout, cex = cex, layout.ncols = layout.ncols,
    layout.heights = layout.heights, tom.subsample = tom.subsample, ...
  )
}


#'
#'
#' @export
wgcna.plotMultiDendroAndColors <- function(multi_wgcna,
                                           block = 1,
                                           extra.colors = NULL,
                                           show.kme = FALSE,
                                           show.traits = FALSE,
                                           show.contrasts = FALSE,
                                           show.mat = NULL,
                                           show.tom = FALSE,
                                           clust = TRUE,
                                           use.tree = 0,
                                           rm.na = TRUE,
                                           sd.wt = 0,
                                           nmax = -1,
                                           main = NULL,
                                           colorHeight = 0.5,
                                           marAll = c(0.4, 5, 1, 0.2),
                                           cex = 1) {

  nw <- length(multi_wgcna)
  par(cex = cex)

  if (is.null(main)) {
    main <- names(multi_wgcna)
  }

  for (k in 1:nw) {
    wgcna.plotDendroAndColors(
      multi_wgcna[[k]],
      marAll = marAll,
      show.traits = show.traits,
      show.contrasts = show.contrasts,
      show.kme = show.kme,
      show.mat = show.mat,
      show.tom = show.tom,
      use.tree = use.tree,
      clust = clust,
      sd.wt = sd.wt,
      nmax = nmax,
      setLayout = (k==1),
      layout.ncols = nw,
      main = main[k],
      cex = cex
    )
  }
}


#'
#'
#' @export
wgcna.plotDendroAndTraitCorrelation <- function(wgcna,
                                                traits = NULL,
                                                show.traits = TRUE,
                                                show.contrasts = TRUE,
                                                show.kme = FALSE,
                                                main = NULL,
                                                block = NULL,
                                                rm.na = TRUE,
                                                use.tree = 0,
                                                marAll = c(0.2, 8, 2, 0.2),
                                                setLayout = TRUE,
                                                ...) {
  message("DEPRECATED: please use wgcna.plotDendroAndColors")

  ## if consensus output do this
  is.cons <- ("class" %in% names(wgcna) && wgcna$class == "cons")
  is.cons2 <- (all(c("layers", "zlist") %in% names(wgcna)))
  if (is.cons || is.cons2) {
    message("object is consensus result")
    wgcna.plotDendroAndTraitCorrelation_cons(
      cons = wgcna,
      traits = traits,
      main = main,
      rm.na = rm.na,
      show.traits = show.traits,
      show.contrasts = show.contrasts,
      marAll = marAll,
      use.tree = use.tree,
      setLayout = setLayout,
      ...
    )
    return()
  }

  moduleColors <- cbind(wgcna$net$colors)
  if (NCOL(moduleColors) == 1) colnames(moduleColors) <- "Module"
  colors <- moduleColors

  if (show.traits) {
    X <- wgcna$datExpr
    Y <- wgcna$datTraits
    sel <- grep("_vs_", colnames(Y), invert = TRUE)
    if (length(sel)) {
      traitSig <- cor(X, Y[, sel], use = "pairwise")
      if (rm.na) {
        sel <- colMeans(is.na(traitSig)) < 1
        traitSig <- traitSig[, sel, drop = FALSE]
      }
      traitColors <- rho2bluered(traitSig)
      colors <- cbind(colors, 0, traitColors)
    }
  }
  if (show.contrasts) {
    X <- wgcna$datExpr
    Y <- wgcna$datTraits
    sel <- grep("_vs_", colnames(Y))
    if (length(sel)) {
      traitSig <- cor(X, Y[, sel], use = "pairwise")
      if (rm.na) {
        sel <- colMeans(is.na(traitSig)) < 1
        traitSig <- traitSig[, sel, drop = FALSE]
      }
      traitColors <- rho2bluered(traitSig)
      colors <- cbind(colors, 0, traitColors)
    }
  }

  geneTree <- wgcna$net$dendrograms[[1]]
  geneTree$labels <- names(wgcna$net$colors)
  colors <- colors[geneTree$labels, ]

  if (is.null(main)) main <- "Gene Dendrogram, Modules and Trait Correlation"

  WGCNA::plotDendroAndColors(
    geneTree,
    colors = colors,
    colnames(colors),
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    marAll = marAll,
    main = main,
    setLayout = setLayout,
    ...
  )
}

#' wgcna.plotDendroAndTraits for Consensus output
#'
#'
wgcna.plotDendroAndTraitCorrelation_cons <- function(cons,
                                                     show.traits = TRUE,
                                                     show.contrasts = TRUE,
                                                     traits = NULL,
                                                     main = NULL,
                                                     rm.na = TRUE,
                                                     use.tree = 0,
                                                     marAll = c(0.2, 8, 2, 0.2),
                                                     setLayout = TRUE,
                                                     ...) {
  message("DEPRECATED: please use wgcna.plotDendroAndColors")

  if (0) {
    show.traits <- TRUE
    traits <- NULL
    main <- NULL
    rm.na <- TRUE
    use.tree <- 0
    marAll <- c(0.2, 8, 2, 0.2)
    setLayout <- TRUE
  }

  ## quick hack to use wgcna.plotDendroAndTraitCorrelation_multi()
  multi <- c(list(Consensus = cons), cons$layers)
  use.tree0 <- use.tree
  if (use.tree %in% 0:99) use.tree <- as.integer(use.tree)
  if (is.character(use.tree)) {
    use.tree <- match(use.tree, names(multi))
  } else {
    use.tree <- as.integer(use.tree) + 1
  }
  if (is.na(use.tree)) {
    message("ERROR: invalid class(use.tree) = ", class(use.tree0))
    message("ERROR: invalid use.tree = ", use.tree0)
    return(NULL)
  }

  wgcna.plotDendroAndTraitCorrelation_multi(
    multi,
    show.traits = show.traits,
    show.contrasts = show.contrasts,
    traits = traits,
    main = main,
    rm.na = rm.na,
    use.tree = use.tree,
    marAll = marAll,
    setLayout = setLayout,
    ...
  )
}


#' @export
wgcna.plotDendroAndTraitCorrelation_multi <- function(multi,
                                                      show.traits = TRUE,
                                                      show.contrasts = TRUE,
                                                      traits = NULL,
                                                      main = NULL,
                                                      rm.na = TRUE,
                                                      use.tree = 1,
                                                      marAll = c(0.2, 8, 2, 0.2),
                                                      setLayout = TRUE,
                                                      ...) {
  message("DEPRECATED: please use wgcna.plotDendroAndColors")

  ## module colors
  colors <- sapply(multi, function(m) m$net$colors)

  if (show.traits || show.contrasts) {
    traitSig <- list()
    nsets <- length(multi)
    i <- 1
    for (k in names(multi)) {
      if (k == "Consensus") next
      w <- multi[[k]]
      Y <- w$datTraits
      sel1 <- sel2 <- NULL
      if (show.traits) sel1 <- grep("_vs_", colnames(Y), invert = TRUE)
      if (show.contrasts) sel2 <- grep("_vs_", colnames(Y))
      sel <- c(sel1, sel2)
      if (!is.null(traits)) sel <- intersect(sel, traits)
      X <- w$datExpr
      kk <- intersect(rownames(X), rownames(Y))
      traitSig[[k]] <- cor(X[kk, ], Y[kk, sel], use = "pairwise")
    }

    if (rm.na) {
      for (i in 1:length(traitSig)) {
        sel <- colMeans(is.na(traitSig[[i]])) < 1
        traitSig[[i]] <- traitSig[[i]][, sel, drop = FALSE]
      }
    }

    ## prepend datatype/set name
    for (k in names(traitSig)) {
      colnames(traitSig[[k]]) <- paste0(k, ":", colnames(traitSig[[k]]))
    }

    traitSig2 <- c()
    for (i in 1:length(traitSig)) {
      traitSig2 <- cbind(traitSig2, traitSig[[i]])
      if (i < length(traitSig)) traitSig2 <- cbind(traitSig2, 0)
    }
    traitColors <- rho2bluered(traitSig2, f = 0.95)
    ii <- which(colnames(traitColors) == "")
    if (length(ii)) traitColors[, ii] <- "#FFFFFF"
    if (is.null(colors)) {
      colors <- traitColors
    } else {
      colors <- cbind(colors, 0, traitColors)
    }
  }

  message("using tree of layer: ", names(multi)[use.tree])
  geneTree <- multi[[use.tree]]$net$dendrograms[[1]]

  if (is.null(main)) main <- "Gene Dendrogram, Modules and Trait Correlation"

  WGCNA::plotDendroAndColors(
    geneTree,
    colors = colors,
    colnames(colors),
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    marAll = marAll,
    main = main,
    ...
  )
}

#' Purple-grey-yellow color ramp
#' @seealso WGCNAplus::purpleGreyYellow
#' @export
purpleGreyYellow <- function(n) {
  WGCNAplus::purpleGreyYellow(n)
}

#' Converts correlation values [-1;1] to blue-white-red colors. Good
#' for creating color labels for labeledHeatmaps that expect colors.
#' @seealso WGCNAplus::rho2bluered
rho2bluered <- function(R, a = 1, f = 0.95) {
  WGCNAplus::rho2bluered(R, a = a, f = f)
}


#' @export
wgcna.plotModuleTraitHeatmap <- function(wgcna, setpar = TRUE, cluster = FALSE,
                                         multi = FALSE, main = NULL, justdata = FALSE,
                                         transpose = FALSE, colorlabel = TRUE,
                                         show = c("both", "traits", "contrasts")[1],
                                         nmax = -1, tmax = -1,
                                         text = TRUE, pstar = TRUE) {
  if (!multi) wgcna <- list(wgcna)
  MEs <- lapply(wgcna, function(w) as.matrix(w$net$MEs))
  MEs <- wgcna.mergeME(MEs)

  Y <- wgcna[[1]]$datTraits
  sel <- 1:ncol(Y)
  if (show == "traits") sel <- grep("_vs_", colnames(Y), invert = TRUE)
  if (show == "contrasts") sel <- grep("_vs_", colnames(Y))
  Y <- Y[, sel, drop = FALSE]

  moduleTraitCor <- cor(MEs, Y, use = "pairwise.complete")

  # nSamples <- nrow(wgcna[[1]]$datExpr)
  nSamples <- t(!is.na(MEs)) %*% (!is.na(Y))

  if (nmax > 0) {
    sel <- head(order(-apply(abs(moduleTraitCor), 1, max, na.rm = TRUE)), nmax)
    moduleTraitCor <- moduleTraitCor[sel, , drop = FALSE]
    nSamples <- nSamples[sel, , drop = FALSE]
  }
  if (tmax > 0) {
    sel <- head(order(-apply(abs(moduleTraitCor), 2, max, na.rm = TRUE)), tmax)
    moduleTraitCor <- moduleTraitCor[, sel, drop = FALSE]
    nSamples <- nSamples[, sel, drop = FALSE]
  }

  if (transpose) {
    moduleTraitCor <- t(moduleTraitCor)
    nSamples <- t(nSamples)
  }

  wgcna.plotLabeledCorrelationHeatmap(
    R = moduleTraitCor,
    nSamples = nSamples,
    setpar = setpar,
    cluster = cluster,
    text = text,
    main = main,
    justdata = justdata,
    colorlabel = colorlabel,
    pstar = pstar
  )
}


#' Plot cluster dendrogram with eigengenes and traits.
#'
#' @export
wgcna.plotEigenGeneClusterDendrogram <- function(wgcna = NULL,
                                                 ME = NULL,
                                                 add_traits = TRUE,
                                                 horiz = FALSE,
                                                 setMargins = TRUE,
                                                 method = "wgcna",
                                                 showlabels = TRUE,
                                                 plot = TRUE,
                                                 multi = FALSE,
                                                 main = NULL) {
  # Matrix with eigengenes and traits
  if (is.null(wgcna) && is.null(ME)) {
    stop("ERROR: wgcna or ME must be given")
  }

  if (is.null(ME)) {
    if (multi) {
      ME <- lapply(wgcna, function(w) as.matrix(w$net$MEs))
      ME <- wgcna.mergeME(ME)
      Y <- wgcna[[1]]$datTraits
    } else {
      ME <- wgcna$net$MEs
      Y <- wgcna$datTraits
    }

    if (length(add_traits) == 1 && is.logical(add_traits) && add_traits == TRUE) {
      ME <- wgcna.mergeME(ME, Y)
    } else if (length(add_traits) > 0 && !is.logical(add_traits)) {
      sel <- intersect(add_traits, colnames(Y))
      if (length(sel)) ME <- wgcna.mergeME(ME, Y[, sel])
    }
  }

  impME <- svdImpute2(as.matrix(ME))
  ME <- WGCNA::orderMEs(impME)

  if (NCOL(ME) <= 2) ME <- cbind(ME, ME) ## error if ncol(ME)<=2 !!!!
  if (is.null(main)) main <- "Eigengene Dendrogram"

  hc <- NULL
  if (method == "wgcna") {
    ## plot dendrogram with WGCNA function
    WGCNA::plotEigengeneNetworks(
      ME, main,
      setMargins = setMargins,
      marDendro = c(0, 4, 2, 0),
      plotHeatmaps = FALSE
    )
  } else {
    ## plot dendrogram with hclust function
    if (setMargins && horiz) par(mar = c(4, 4, 4, 8))
    if (setMargins && !horiz) par(mar = c(8, 4, 4, 1))
    hc <- fastcluster::hclust(as.dist(1 - cor(ME)), method = "average")
    if (plot) {
      save.labels <- hc$labels
      if (!showlabels) hc$labels <- rep("", ncol(ME))
      plot(as.dendrogram(hc), horiz = horiz, main = main)
      hc$labels <- save.labels
    }
  }
  invisible(hc)
}


#' Plot the adjacency correlation heatmap matrix of eigengenes with or
#' without traits. This can show how traits cluster together with the
#' eigengenes.
#'
#' @export
wgcna.plotEigenGeneAdjacencyHeatmap <- function(wgcna,
                                                add_traits = TRUE,
                                                traits = NULL,
                                                add_me = TRUE,
                                                marx = 1, main = NULL,
                                                multi = FALSE,
                                                phenotype = NULL,
                                                colorlabel = TRUE,
                                                text = FALSE,
                                                pstar = TRUE,
                                                power = 1,
                                                setMargins = TRUE,
                                                mar1 = c(5.6, 4.5, 1.8, 0),
                                                mar2 = c(8, 10, 4, 2),
                                                cex.lab = 0.8,
                                                cex.text = 0.6,
                                                plotDendro = TRUE,
                                                plotHeatmap = TRUE,
                                                dendro.horiz = TRUE,
                                                dendro.width = 0.3,
                                                dendro.labels = TRUE,
                                                nmax = -1,
                                                fixclust = FALSE,
                                                mask.intra = FALSE,
                                                justdata = FALSE) {
  if (0) {
    add_traits <- TRUE
    traits <- NULL
    marx <- 1
    main <- NULL
    multi <- FALSE
    phenotype <- NULL
    colorlabel <- TRUE
    text <- FALSE
    pstar <- TRUE
    setMargins <- TRUE
    mar1 <- c(5.5, 5, 1.6, 1)
    mar2 <- c(8, 10, 4, 2)
    cex.lab <- 0.8
    cex.text <- 0.7
    plotDendro <- TRUE
    plotHeatmap <- TRUE
    dendro.horiz <- TRUE
    dendro.width <- 0.3
    dendro.labels <- TRUE
    nmax <- -1
    fixclust <- FALSE
    mask.intra <- FALSE
    justdata <- FALSE
    add_me <- TRUE
  }

  if (!multi) wgcna <- list(gx = wgcna)

  # Matrix with eigengenes and traits
  ME <- NULL
  if (add_me) {
    ME <- lapply(wgcna, function(w) as.matrix(w$net$MEs))
    ME <- wgcna.mergeME(ME)
  }
  Y <- wgcna[[1]]$datTraits

  if (add_traits) {
    sel <- colnames(Y)
    if (!is.null(traits)) {
      sel <- intersect(traits, sel)
    }
    if (is.null(ME)) {
      ME <- Y[, sel, drop = FALSE]
    } else {
      ME <- wgcna.mergeME(ME, Y[, sel, drop = FALSE])
    }
  }

  if (!add_traits && !is.null(phenotype)) {
    if (is.null(ME)) {
      ME <- Y[, phenotype, drop = FALSE]
    } else {
      ME <- wgcna.mergeME(ME, Y[, phenotype, drop = FALSE])
    }
  }

  if (NCOL(ME) <= 2) ME <- cbind(ME, ME) ## error if ncol(ME)<=2 !!!!

  ## Compute eigengene correlation matrix. Repeat 'power' times for
  ## higher order adjacency.
  power <- round(power)
  R <- ME
  for (i in 1:power) {
    tt <- cortest(R, R)
    R <- tt$rho
    nSamples <- tt$n
  }

  ## If phenotype is given, we condition the heatmap using the
  ## correlation to the phenotype.
  if (!is.null(phenotype)) {
    ## proper sign in case of inhibitor layer (like miRNA)
    layersign <- rep(1, length(wgcna))
    names(layersign) <- names(wgcna)
    layersign[grep("^mi", names(wgcna), ignore.case = TRUE)] <- -1
    ff <- list()
    for (k in names(wgcna)) {
      rho <- cor(ME, Y[, phenotype], use = "pairwise")[, 1]
      ff[[k]] <- layersign[k] * rho
    }
    names(ff) <- NULL
    ff <- unlist(ff)
    ff <- 0.5 * (1 + ff) ## signed...
    ff <- ff[match(rownames(R), names(ff))]
    names(ff) <- rownames(R)
    ff[is.na(ff)] <- 1 ## really??? NEED RETHINK
    ww <- outer(ff, ff)
    ## ww[is.na(ww)] <- 0
    ww <- ww / max(ww, na.rm = TRUE)
    R <- R * ww
  }

  if (nmax > 0) {
    if (!is.null(phenotype)) {
      y1 <- Y[, phenotype]
      y1 <- y1[match(rownames(ME), names(y1))]
      rho <- cor(ME, y1, use = "pairwise")[, 1]
      ii <- head(order(-abs(rho)), nmax)
    } else {
      ii <- head(order(-Matrix::rowMeans(R**2)), nmax)
    }
    R <- R[ii, ii]
  }

  if (justdata) {
    return(R)
  }

  # Plot the correlation heatmap matrix (note: this plot will overwrite
  # the dendrogram plot)
  if (is.null(main)) main <- "Eigengene Adjacency Heatmap"

  if (plotDendro && plotHeatmap) {
    layout.matrix <- matrix(1:2, nrow = 1, ncol = 2)
    layout(mat = layout.matrix, heights = 1, widths = c(dendro.width, 1))
    if (dendro.horiz && dendro.labels) {
      mar1[4] <- mar2[2] ## copy left margin
    }
  }
  if (plotDendro) par(mar = mar1)

  # fixclust=FALSE
  R0 <- R
  R0[is.na(R0)] <- 0

  if (fixclust) {
    ii <- rownames(R)
    hc <- fastcluster::hclust(as.dist(1 - R0[ii, ii]), method = "average")
  } else {
    hc <- fastcluster::hclust(as.dist(1 - R0), method = "average")
  }
  if (plotDendro) {
    par(cex = cex.lab)
    plot(as.dendrogram(hc),
      horiz = TRUE,
      ylab = "Eigengene dendrogram"
    )
    par(cex = 1)
  }

  if (plotHeatmap) {
    ii <- hc$labels[hc$order]
    ii <- intersect(ii, rownames(R))
    R1 <- R[rev(ii), ii]
    # nsamples <- nrow(Y)
    nsamples <- nSamples[rownames(R1), colnames(R1)]
    par(mar = mar2)
    wgcna.plotLabeledCorrelationHeatmap(
      R1,
      nSamples = nsamples,
      text = text,
      pstar = pstar,
      colorlabel = colorlabel,
      cluster = FALSE,
      setpar = FALSE,
      main = main,
      cex.lab = cex.lab,
      cex.text = cex.text
    )
  }
  invisible(R)
}

#' @export
wgcna.plotMultiEigengeneCorrelation <- function(wgcna, addtraits = TRUE,
                                                phenotype = NULL, nmax = -1, main = NULL,
                                                showvalues = FALSE, showsig = TRUE,
                                                cex.text = 0.6, cex.lab = 0.8,
                                                fixcluster = TRUE, setpar = TRUE) {
  ## Show inter-correlation of modules
  me <- lapply(wgcna, function(w) w$net$MEs)
  if (length(me) == 1) {
    me <- list(me[[1]], me[[1]])
  }

  comb <- combn(length(me), 2)
  ncomb <- ncol(comb)
  nsamples <- nrow(wgcna[[1]]$datExpr)
  Y <- wgcna[[1]]$datTraits

  ## for miRNA we need to flip sign
  msign <- c(1, -1)[1 + 1 * (names(wgcna) %in% c("mi", "mir"))]

  if (setpar) {
    nc <- ceiling(sqrt(ncomb))
    nr <- ceiling(ncomb / nc)
    par(mfrow = c(nr, nc), mar = c(8, 10, 3, 1))
  }

  k <- 1
  for (k in 1:ncol(comb)) {
    i <- comb[1, k]
    j <- comb[2, k]
    M1 <- me[[i]]
    M2 <- me[[j]]

    if (addtraits) {
      M1 <- cbind(M1, Y)
      M2 <- cbind(M2, Y)
    }
    if (FALSE && !addtraits && !is.null(phenotype)) {
      y <- Y[, phenotype, drop = FALSE]
      M1 <- cbind(M1, y)
      M2 <- cbind(M2, y)
    }


    R1 <- cor(M1, M2, use = "pairwise.complete")

    if (nmax > 0) {
      ii <- head(order(-apply(abs(R1), 1, max)), nmax)
      jj <- head(order(-apply(abs(R1), 2, max)), nmax)
      R1 <- R1[ii, jj]
    }

    ## cluster unweighted matrix
    ii <- fastcluster::hclust(dist(R1), method = "average")$order
    jj <- fastcluster::hclust(dist(t(R1)), method = "average")$order
    R1 <- R1[ii, jj]

    ## This conditions the correlation on phenotype. Important.
    do.condition <- !is.null(phenotype)
    if (do.condition) {
      y <- Y[, phenotype]
      w1 <- cor(M1[, rownames(R1)], y, use = "pairwise")[, 1]
      w2 <- cor(M2[, colnames(R1)], y, use = "pairwise")[, 1]
      if (msign[i] != 0) w1 <- msign[i] * w1
      if (msign[j] != 0) w2 <- msign[j] * w2
      w1 <- pmax(w1, 0)
      w2 <- pmax(w2, 0)
      ww <- outer(w1, w2)
      ww <- ww / max(ww, na.rm = TRUE)
      R1 <- R1 * ww
    }

    main <- paste(names(me)[i], "vs.", names(me)[j])
    if (do.condition) main <- paste(main, "(conditioned)")

    wgcna.plotLabeledCorrelationHeatmap(
      R1,
      nsamples,
      text = showvalues,
      pstar = showsig,
      is.dist = FALSE,
      cluster = !fixcluster,
      setpar = FALSE,
      main = main,
      cex.text = cex.text,
      cex.lab = cex.lab
    )
  }
}


#' @export
wgcna.plotEigenGeneGraph <- function(wgcna, add_traits = TRUE, main = NULL,
                                     multi = FALSE, vcex = 1, labcex = 1,
                                     mincor = 0.5, as.phylo = TRUE) {
  ## require(igraph)
  if (multi) {
    ME <- lapply(wgcna, function(w) as.matrix(w$net$MEs))
    ME <- wgcna.mergeME(ME)
    if (add_traits) ME <- cbind(ME, wgcna[[1]]$datTraits)
  } else {
    ME <- wgcna$net$MEs
    if (add_traits) ME <- cbind(ME, wgcna$datTraits)
  }
  if (NCOL(ME) <= 2) ME <- cbind(ME, ME) ## error if ncol(ME)<=2 !!!!

  sdx <- matrixStats::colSds(as.matrix(ME * 1), na.rm = TRUE)
  if (any(sdx == 0)) ME <- ME + runif(length(ME), 0, 1e-5)

  ## Recalculate MEs with color as labels
  corx <- cor(ME, use = "pairwise")
  corx[is.na(corx)] <- 0

  layout <- NULL
  if(as.phylo) {
    clust <- fastcluster::hclust(as.dist(1 - corx))
    phylo <- ape::as.phylo(clust)
    gr <- igraph::as.igraph(phylo, directed = FALSE)
    layout <- igraph::layout_with_kk
  } else {
    adj <- (corx > mincor) * corx
    gr <- igraph::graph_from_adjacency_matrix(adj, diag=FALSE,
      mode="undirected", weighted=TRUE)
    layout <- igraph::layout_with_fr
  }
  
  is.node <- grepl("Node", igraph::V(gr)$name)
  module.name <- igraph::V(gr)$name
  if (multi) {
    module.size <- lapply(wgcna, function(w) table(w$net$labels))
    names(module.size) <- NULL
    module.size <- unlist(module.size)
    module.colors <- sapply(wgcna, function(w) w$me.colors)
    names(module.colors) <- NULL
    module.colors <- unlist(module.colors)
  } else {
    module.size <- table(wgcna$net$labels)
    module.colors <- wgcna$me.colors
  }
  module.size <- module.size / mean(module.size)
  module.colors <- adjustcolor(module.colors[module.name],alpha.f=0.3)
  igraph::V(gr)$label <- igraph::V(gr)$name
  igraph::V(gr)$label[is.node] <- NA
  igraph::V(gr)$color <- module.colors
  igraph::V(gr)$size <- vcex * 18 * (module.size[module.name])**0.4
  igraph::V(gr)$size[is.na(igraph::V(gr)$size)] <- 0

  ## par(mfrow = c(1, 1), mar = c(1, 1, 1, 1) * 0)
  igraph::plot.igraph(
    gr,
    layout = layout,
    vertex.label.cex = 0.85 * labcex,
    edge.width = 3
  )
  if (!is.null(main)) title(main, line = -1.5)
}


#' Plot Multi-dimensional scaling (MDS) of centered data matrix.
#'
#' @export
wgcna.plotMDS <- function(wgcna, main = NULL, scale = FALSE) {
  cc <- wgcna.labels2colors(wgcna$net$color)
  pc <- svd(t(scale(wgcna$datExpr, scale = scale)), nv = 2)$u[, 1:2]
  # pc <- svd(t(scale(wgcna$datExpr)),nv=1)$u[,1:2]
  colnames(pc) <- c("MDS-x", "MDS-y")
  if (is.null(main)) main <- "MDS of features"
  plot(pc, col = cc, main = main)
}

#'
#'
#' @export
wgcna.plotFeatureUMAP <- function(wgcna, nhub = 3, method = "clust",
                                  scale = FALSE, main = NULL,
                                  plotlib = "base", annot = NULL) {
  if (method == "clust" && "clust" %in% names(wgcna)) {
    pos <- wgcna$clust[["umap2d"]]
  } else if (method == "umap") {
    cX <- t(scale(wgcna$datExpr, scale = scale)) ## WGCNA uses correlation
    pos <- uwot::umap2(cX)
    colnames(pos) <- c("UMAP-x", "UMAP-y")
    rownames(pos) <- colnames(wgcna$datExpr)
    ##  } else if(method=="mds") {
  } else {
    pos <- svd(t(scale(wgcna$datExpr, scale = scale)), nv = 2)$u[, 1:2]
    colnames(pos) <- c("MDS-x", "MDS-y")
    rownames(pos) <- colnames(wgcna$datExpr)
  }

  if (is.null(main)) main <- "Feature UMAP colored by module"

  hubgenes <- NULL
  if (nhub > 0) {
    ## get top hub genes
    mm <- wgcna$stats$moduleMembership
    hubgenes <- apply(mm, 2, function(x) head(names(sort(-x)), nhub), simplify = FALSE)
    sel <- which(names(hubgenes) != "MEgrey")
    hubgenes <- unlist(hubgenes[sel])
  }

  col1 <- wgcna$net$colors
  genes1 <- names(which(col1 != "grey"))
  if (!is.null(annot)) {
    rownames(pos) <- probe2symbol(rownames(pos), annot, "gene_name", fill_na = TRUE)
    names(col1) <- probe2symbol(names(col1), annot, "gene_name", fill_na = TRUE)
    genes1 <- probe2symbol(genes1, annot, "gene_name", fill_na = TRUE)
    if (nhub > 0) {
      hubgenes <- setNames(probe2symbol(hubgenes, annot, "gene_name", fill_na = TRUE), names(hubgenes))
    }
  }
  pgx.scatterPlotXY(
    pos,
    var = col1,
    col = sort(unique(col1)),
    hilight = genes1,
    hilight2 = hubgenes,
    cex.lab = 1.2,
    label.clusters = FALSE,
    title = main,
    plotlib = plotlib
  )
}


#' Plot module significance.
#'
#' @export
wgcna.plotModuleSignificance <- function(wgcna, trait, main = NULL, abs = FALSE) {
  ## cc <- paste0("ME",wgcna$net$color)
  cc <- wgcna.labels2colors(wgcna$net$color)
  if ("stats" %in% names(wgcna)) {
    traitSignificance <- wgcna$stats$traitSignificance
  } else {
    traitSignificance <- as.data.frame(cor(wgcna$datExpr, wgcna$datTraits, use = "p"))
    names(traitSignificance) <- names(wgcna$datTraits)
    rownames(traitSignificance) <- colnames(wgcna$datExpr)
  }
  geneSig <- traitSignificance[, trait]
  if (is.null(main)) main <- paste("Module significance with", trait)
  if (abs) geneSig <- abs(geneSig)
  WGCNA::plotModuleSignificance(
    geneSig,
    colors = cc, main = main, boxplot = FALSE
  )
}


#'
#'
#' @export
wgcna.plotConsensusSampleDendroAndColors <- function(cons, i,
                                                     what = c("both", "me", "traits")[1],
                                                     show.me = TRUE, show.traits = TRUE,
                                                     show.contrasts = TRUE,
                                                     clust.expr = TRUE,
                                                     setLayout = TRUE,
                                                     marAll = c(0.2, 7, 1.5, 0.5),
                                                     colorHeightMax = 0.6,
                                                     main = NULL) {
  wgcna.plotSampleDendroAndColors(
    wgcna = cons$layers[[i]],
    main = main,
    datExpr = cons$datExpr[[i]],
    datTraits = cons$datTraits,
    datME = cons$net$multiME[[i]]$data,
    what = what,
    show.me = show.me,
    show.traits = show.traits,
    show.contrasts = show.contrasts,
    marAll = marAll,
    clust.expr = clust.expr,
    setLayout = setLayout,
    colorHeightMax = colorHeightMax
  )
}

#'
#'
#' @export
wgcna.plotSampleDendroAndColors <- function(wgcna, input.type = "wgcna",
                                            what = c("me", "traits", "both")[3],
                                            show.me = TRUE, show.traits = TRUE,
                                            show.contrasts = TRUE,
                                            datTraits = NULL, datExpr = NULL, datME = NULL,
                                            clust.expr = TRUE, setLayout = TRUE,
                                            marAll = c(0.2, 7, 1.5, 0.5),
                                            colorHeightMax = 0.6,
                                            main = NULL, justdata = FALSE) {
  if (input.type == "net") {
    ME0 <- wgcna$MEs
    if (is.null(datExpr)) stop("must supply datExpr")
    if (is.null(datTraits)) stop("must supply datTraits")
  } else {
    ME0 <- wgcna$net$MEs
    datTraits <- 1 * wgcna$datTraits
    datExpr <- wgcna$datExpr
  }

  if (!is.null(datME)) {
    ME0 <- datME
  }

  ME <- ME0[, 0]
  samples <- rownames(ME)
  if (show.me) {
    ME <- cbind(ME, ME0)
  }
  if (show.traits) {
    sel <- grep("_vs_", colnames(datTraits), invert = TRUE)
    ME <- cbind(ME, datTraits[samples, sel, drop = FALSE])
  }
  if (show.contrasts) {
    sel <- grep("_vs_", colnames(datTraits))
    ME <- cbind(ME, datTraits[samples, sel, drop = FALSE])
  }

  if (NCOL(ME) <= 2) ME <- cbind(ME, ME) ## error if ncol(ME)<=2 !!!!
  sdx <- matrixStats::colSds(as.matrix(ME * 1), na.rm = TRUE)
  ME <- ME[, which(sdx > 0), drop = FALSE]

  ## Recalculate MEs with color as labels
  if (clust.expr) {
    corx <- cor(t(datExpr), use = "pairwise")
  } else {
    corx <- cor(t(ME0), use = "pairwise")
  }
  corx[is.na(corx)] <- 0
  sampleTree <- fastcluster::hclust(as.dist(1 - corx), method = "average")

  corx <- cor(ME, use = "pairwise")
  corx[is.na(corx)] <- 0
  jj <- fastcluster::hclust(as.dist(1 - corx))$order
  colors <- WGCNA::numbers2colors(ME[, jj])

  if (justdata) {
    return(ME)
  }

  if (is.null(main)) {
    if (what == "me") main <- "Sample dendrogram and module heatmap"
    if (what == "traits") main <- "Sample dendrogram and trait heatmap"
    if (what == "both") main <- "Sample dendrogram and module+traits heatmap"
  }

  ## Plot the dendrogram and the module colors underneath
  WGCNA::plotDendroAndColors(
    dendro = sampleTree,
    colors = colors,
    groupLabels = colnames(ME)[jj],
    dendroLabels = rownames(ME),
    hang = 0.03,
    addGuide = FALSE,
    guideHang = 0.05,
    setLayout = setLayout,
    marAll = marAll,
    main = main,
    colorHeightMax = colorHeightMax
  )
}


#' @seealso WGCNAplus::plotLabeledCorrelationHeatmap
#' @export
wgcna.plotLabeledCorrelationHeatmap <- function(R, nSamples,
                                                cluster = TRUE, text = TRUE,
                                                main = NULL, justdata = FALSE,
                                                colorlabel = TRUE, pstar = TRUE,
                                                zlim = NULL, colorpal = NULL,
                                                cex.text = 0.6, cex.lab = NULL,
                                                setpar = TRUE, is.dist = FALSE) {
  WGCNAplus::plotLabeledCorrelationHeatmap(
    R = R, nSamples = nSamples, cluster = cluster, text = text, main = main,
    justdata = justdata, colorlabel = colorlabel, pstar = pstar, zlim = zlim,
    colorpal = colorpal, cex.text = cex.text, cex.lab = cex.lab,
    setpar = setpar, is.dist = is.dist
  )
}


#' Plot module hub genes as graph in circle layout. This shows the
#' main connectivity structure for the modules. It can be used for
#' (edge) preservation analysis for comparing different group of
#' samples..
#'
#' @export
wgcna.plotModuleHubGenes <- function(wgcna, modules = NULL,
                                     alpha = 0.5, setpar = TRUE) {
  if (is.null(modules)) {
    modules <- colnames(wgcna$stats$moduleMembership)
  }
  modules <- intersect(modules, colnames(wgcna$stats$moduleMembership))
  if (length(modules) == 0) {
    message("ERROR. no valid modules")
    return(NULL)
  }

  if (setpar) {
    nr <- floor(length(modules)**0.5)
    nc <- ceiling(length(modules) / nr)
    nr
    nc
    par(mfrow = c(nr, nc), mar = c(0, 1, 2.5, 1))
  }
  for (k in modules) {
    mm <- wgcna$stats$moduleMembership[, k]
    mm.score <- head(sort(mm, decreasing = TRUE), 30)
    topgenes <- names(mm.score)
    A <- cor(wgcna$datExpr[, topgenes])
    diag(A) <- 0
    A <- (A - min(A, na.rm = TRUE)) / (max(A, na.rm = TRUE) - min(A, na.rm = TRUE))
    A[is.na(A)] <- 0
    gr <- igraph::graph_from_adjacency_matrix(
      A,
      mode = "undirected", weighted = TRUE, diag = FALSE
    )
    norm.mm.score <- (mm.score - min(mm.score)) / (max(mm.score) - min(mm.score))
    clr <- sub("ME", "", k)
    if (!is.na(as.integer(clr))) clr <- as.integer(clr)
    if (clr == "black") clr <- "grey40"
    plot(gr,
      layout = igraph::layout_in_circle,
      edge.width = 6 * igraph::E(gr)$weight**8,
      vertex.size = 5 + 15 * norm.mm.score,
      vertex.color = adjustcolor(clr, alpha.f = alpha),
      vertex.frame.color = clr
    )
    title(paste("Module", k), line = 0.33)
  }
}

#' @export
wgcna.plotGeneNetwork <- function(wgcna, genes, col = NULL,
                                  edge.alpha = 0.3,
                                  rgamma = 4,
                                  edge.width = 6,
                                  alpha = 0.5,
                                  min.rho = 0.5,
                                  setpar = TRUE) {
  A <- cor(wgcna$datExpr[, genes])
  A <- A * (abs(A) > min.rho)
  A[is.na(A)] <- 0
  gr <- igraph::graph_from_adjacency_matrix(
    A,
    mode = "undirected", weighted = TRUE, diag = FALSE
  )
  vcex <- matrixStats::colVars(wgcna$datExpr[, genes])
  vcex <- vcex / max(abs(vcex))
  if (is.null(col)) {
    col <- wgcna$net$color[genes]
  }
  col <- sub("black", "grey40", col)
  ecol <- c("darkred", "darkgreen")[1 + 1 * (igraph::E(gr)$weight > 0)]
  table(ecol)
  ecol <- adjustcolor(ecol, alpha.f = edge.alpha)
  ewt <- abs(igraph::E(gr)$weight)
  ewt <- (ewt / max(abs(ewt)))**rgamma
  plot(gr,
    layout = igraph::layout_in_circle,
    edge.width = edge.width * ewt,
    edge.color = ecol,
    vertex.size = 5 + 20 * vcex,
    vertex.color = adjustcolor(col, alpha.f = alpha),
    vertex.frame.color = col
  )
}

#' @seealso WGCNAplus::plotModuleHeatmap
#' @export
wgcna.plotModuleHeatmap <- function(wgcna,
                                    module,
                                    genes = NULL,
                                    rgamma = 4,
                                    min.rho = 0,
                                    cex = 0.8,
                                    nmax = -1,
                                    cluster = TRUE,
                                    type = c("expression", "correlation")[1],
                                    heatmap.mar = c(7, 7),
                                    main = NULL) {
  WGCNAplus::plotModuleHeatmap(
    wgcna = wgcna, module = module, genes = genes, rgamma = rgamma,
    min.rho = min.rho, cex = cex, nmax = nmax, cluster = cluster,
    type = type, heatmap.mar = heatmap.mar, main = main
  )
}

#'
#'
#' @export
wgcna.plotPowerAnalysis <- function(datExpr, networktype = "signed",
                                    cex = 1, maxpower = 20, nmax = 2000,
                                    plots = c(
                                      "sft.modelfit", "mean.k",
                                      "dendro.IQR"
                                    ), main = NULL,
                                    RsquaredCut = 0.85, setPar = TRUE) {
  RsquaredCut <- RsquaredCut[1]

  ## Choose a set of soft-thresholding powers
  powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
  if (maxpower > 20) {
    powers <- c(powers, seq(from = 20, to = maxpower, by = 5))
  }

  ## subsample for speed
  if (ncol(datExpr) > nmax && nmax > 0) {
    ii <- sample(1:ncol(datExpr), nmax)
    datExpr <- datExpr[, ii]
  }

  ## Call the network topology analysis function
  sft <- WGCNA::pickSoftThreshold(
    datExpr,
    powerVector = powers,
    RsquaredCut = RsquaredCut,
    networkType = networktype,
    verbose = 0
  )

  ## This is more robust

  if (setPar) {
    np <- length(plots)
    nc <- ceiling(sqrt(np))
    par(mfrow = c(nc, nc), mar = c(3.3, 3.5, 1, 1), mgp = c(2, 0.9, 0))
    par(mfrow = c(1, np), mar = c(3.8, 3.8, 1, 1), mgp = c(2.4, 0.95, 0))
  }

  ## Plot the results:
  if ("sft.modelfit" %in% plots) {
    ## Scale-free topology fit index as a function of the soft-thresholding power
    y <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
    base::plot(
      x = sft$fitIndices[, 1],
      y = y,
      ylim = c(min(y), 1),
      type = "n",
      xlab = "Soft threshold (power)",
      ylab = "SFT model fit (signed R^2)",
      main = main
    )
    abline(h = 0, col = "black", lty = 3)
    text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
      labels = powers, cex = cex, col = "red"
    )
    ## this line corresponds to using an R^2 cut-off of h
    abline(h = RsquaredCut, col = "red", lty = 2)
    ## if(legend) legend("bottomright", legend=paste("opt. power =",optPower))
  }

  ## Mean connectivity as a function of the soft-thresholding power
  if ("mean.k" %in% plots) {
    base::plot(sft$fitIndices[, "Power"], sft$fitIndices[, "mean.k."],
      type = "n",
      xlab = "Soft threshold (power)",
      ylab = "Mean connectivity",
      main = main
    )
    text(sft$fitIndices[, "Power"], sft$fitIndices[, "mean.k."],
      labels = powers,
      cex = cex, col = "red"
    )
  }

  ht <- NULL
  if ("dendro.IQR" %in% plots) {
    ht <- wgcna.checkDendroHeights(datExpr, n = 200, powers = powers)
    base::plot(
      sft$fitIndices[, 1], ht$IQR,
      type = "n",
      xlab = "Soft threshold (power)",
      ylab = "Dendrogram height IQR",
      main = main
    )
    text(sft$fitIndices[, 1], ht$IQR,
      labels = powers,
      cex = cex, col = "red"
    )
  }
}

#'
#'
#' @export
wgcna.plotPowerAnalysis_multi <- function(exprList,
                                          cex = 1, maxpower = 20,
                                          nmax = 2000,
                                          networktype = "signed",
                                          plots = c(
                                            "sft.modelfit", "mean.k",
                                            "dendro.IQR"
                                          ),
                                          main = NULL,
                                          cex.legend = 1,
                                          RsquaredCut = 0.85,
                                          setPar = TRUE) {
  if (0) {
    networktype <- "signed"
    cex <- 1
    maxpower <- 20
    nmax <- 2000
    plots <- c("sft.modelfit", "mean.k", "dendro.IQR")
    main <- NULL
    RsquaredCut <- 0.85
  }

  RsquaredCut <- RsquaredCut[1]

  ## Choose a set of soft-thresholding powers
  powers <- c(c(1:10), seq(from = 12, to = 20, by = 2))
  if (maxpower > 20) {
    powers <- c(powers, seq(from = 20, to = maxpower, by = 5))
  }

  ## process each list
  sft <- list()
  for (i in 1:length(exprList)) {
    datExpr <- Matrix::t(exprList[[i]])

    ## subsample for speed
    if (nmax > 0 && nrow(exprList[[i]]) > nmax) {
      ii <- sample(1:ncol(datExpr), nmax)
      datExpr <- datExpr[, ii]
    }

    ## Call the network topology analysis function
    k <- names(exprList)[i]
    sft[[k]] <- WGCNA::pickSoftThreshold(
      datExpr,
      powerVector = powers,
      RsquaredCut = RsquaredCut,
      networkType = networktype,
      verbose = 0
    )
  }

  if (setPar) {
    np <- length(plots)
    par(mfrow = c(1, np), mar = c(3.8, 4.5, 3, 1), mgp = c(2.6, 0.95, 0))
  }

  ## Plot the results:
  if ("sft.modelfit" %in% plots) {
    ## Scale-free topology fit index as a function of the soft-thresholding power
    Y <- c()
    for (i in 1:length(sft)) {
      y1 <- -sign(sft[[i]]$fitIndices[, "slope"]) * sft[[i]]$fitIndices[, "SFT.R.sq"]
      Y <- cbind(Y, y1)
    }
    colnames(Y) <- names(sft)
    x <- sft[[1]]$fitIndices[, "Power"]
    Y <- pmax(Y, 0)
    matplot(
      x = x,
      y = Y,
      ylim = c(0, 1),
      type = "l",
      col = 2:99,
      lty = 1,
      lwd = 0.6,
      xlab = "Soft threshold (power)",
      ylab = "SFT model fit (signed R^2)",
      main = main
    )
    # abline(h = 0, col = "black", lty=1)
    for (i in 1:ncol(Y)) {
      text(powers, Y[, i], labels = "█", cex = cex, col = "white")
      text(powers, Y[, i], labels = powers, cex = cex, col = 1 + i)
    }
    ## this line corresponds to using an R^2 cut-off of h
    abline(h = RsquaredCut, col = "grey10", lty = 2)
    legend("bottomright",
      legend = colnames(Y), fill = 2:10,
      cex = cex.legend, y.intersp = 0.9
    )
    title("SFT model fit")
  }

  ## Mean connectivity as a function of the soft-thresholding power
  if ("mean.k" %in% plots) {
    Y <- sapply(sft, function(s) s$fitIndices[, "mean.k."])
    matplot(
      powers,
      Y,
      type = "l",
      col = 2:99,
      lty = 1,
      lwd = 0.6,
      xlab = "Soft threshold (power)",
      ylab = "Mean connectivity",
      main = main
    )
    for (i in 1:ncol(Y)) {
      text(powers, Y[, i], labels = "█", cex = cex, col = "white")
      text(powers, Y[, i], labels = powers, cex = cex, col = 1 + i)
    }
    legend("topright",
      legend = colnames(Y), fill = 2:10,
      cex = cex.legend, y.intersp = 0.9
    )
    title("Mean connectivity")
  }

  ht <- NULL
  if ("dendro.IQR" %in% plots) {
    ht <- list()
    for (i in 1:length(exprList)) {
      ht[[i]] <- wgcna.checkDendroHeights(
        Matrix::t(exprList[[i]]),
        n = 200, powers = powers
      )
    }
    Y <- sapply(ht, function(h) h$IQR)
    matplot(
      powers,
      Y,
      type = "l",
      col = 2:99,
      lty = 1,
      lwd = 0.6,
      xlab = "Soft threshold (power)",
      ylab = "Dendrogram height IQR",
      main = main
    )
    for (i in 1:ncol(Y)) {
      text(powers, Y[, i], labels = "█", cex = cex, col = "white")
      text(powers, Y[, i], labels = powers, cex = cex, col = 1 + i)
    }
    legend("bottomright",
      legend = names(exprList), fill = 2:10,
      cex = cex.legend, y.intersp = 0.9
    )
    title("Dendrogram IQR")
  }
}
