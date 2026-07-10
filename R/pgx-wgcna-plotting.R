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

#' @seealso WGCNAplus::plotPreservationSummaries
#' @export
wgcna.plotPreservationSummaries <- function(pres, setpar = TRUE) {
  WGCNAplus::plotPreservationSummaries(pres = pres, setpar = setpar)
}

#' @seealso WGCNAplus::plotPreservationModuleTraits
#' @export
wgcna.plotPreservationModuleTraits <- function(pres,
                                               subplots = c("zsummary", "consmt", "wt.consmt"),
                                               order.by = "name",
                                               setpar = TRUE, rm.na = FALSE) {
  WGCNAplus::plotPreservationModuleTraits(
    pres = pres, subplots = subplots, order.by = order.by,
    setpar = setpar, rm.na = rm.na
  )
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

#' @seealso WGCNAplus::plotTraitCorrelationBarPlots
#' @export
wgcna.plotTraitCorrelationBarPlots <- function(res, trait, multi = FALSE,
                                               colored = TRUE, beside = TRUE,
                                               main = NULL, cex.main = 1.3,
                                               setpar = TRUE) {
  WGCNAplus::plotTraitCorrelationBarPlots(
    res = res, trait = trait, multi = multi, colored = colored,
    beside = beside, main = main, cex.main = cex.main, setpar = setpar
  )
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
#' @seealso WGCNAplus::plotEigenGeneClusterDendrogram
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
  WGCNAplus::plotEigenGeneClusterDendrogram(
    wgcna = wgcna, ME = ME, add_traits = add_traits, horiz = horiz,
    setMargins = setMargins, method = method, showlabels = showlabels,
    plot = plot, multi = multi, main = main
  )
}


#' Plot the adjacency correlation heatmap matrix of eigengenes with or
#' without traits. This can show how traits cluster together with the
#' eigengenes.
#' @seealso WGCNAplus::plotEigenGeneAdjacencyHeatmap
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
  WGCNAplus::plotEigenGeneAdjacencyHeatmap(
    wgcna = wgcna, add_traits = add_traits, traits = traits, add_me = add_me,
    marx = marx, main = main, multi = multi, phenotype = phenotype,
    colorlabel = colorlabel, text = text, pstar = pstar, power = power,
    setMargins = setMargins, mar1 = mar1, mar2 = mar2, cex.lab = cex.lab,
    cex.text = cex.text, plotDendro = plotDendro, plotHeatmap = plotHeatmap,
    dendro.horiz = dendro.horiz, dendro.width = dendro.width,
    dendro.labels = dendro.labels, nmax = nmax, fixclust = fixclust,
    mask.intra = mask.intra, justdata = justdata
  )
}

#' @seealso WGCNAplus::plotMultiEigengeneCorrelation
#' @export
wgcna.plotMultiEigengeneCorrelation <- function(wgcna, addtraits = TRUE,
                                                phenotype = NULL, nmax = -1, main = NULL,
                                                showvalues = FALSE, showsig = TRUE,
                                                cex.text = 0.6, cex.lab = 0.8,
                                                fixcluster = TRUE, setpar = TRUE) {
  WGCNAplus::plotMultiEigengeneCorrelation(
    wgcna = wgcna, addtraits = addtraits, phenotype = phenotype, nmax = nmax,
    main = main, showvalues = showvalues, showsig = showsig,
    cex.text = cex.text, cex.lab = cex.lab, fixcluster = fixcluster,
    setpar = setpar
  )
}


#' WGCNAplus::plotEigenGeneGraph() renders module-color node fills at full
#' opacity rather than this function's alpha.f=0.3 transparency - a
#' cosmetic rendering difference only, no data/value change.
#' @seealso WGCNAplus::plotEigenGeneGraph
#' @export
wgcna.plotEigenGeneGraph <- function(wgcna, add_traits = TRUE, main = NULL,
                                     multi = FALSE, vcex = 1, labcex = 1,
                                     mincor = 0.5, as.phylo = TRUE) {
  WGCNAplus::plotEigenGeneGraph(
    wgcna = wgcna, add_traits = add_traits, main = main, multi = multi,
    vcex = vcex, labcex = labcex, mincor = mincor, as.phylo = as.phylo
  )
}


#' Plot Multi-dimensional scaling (MDS) of centered data matrix.
#' @seealso WGCNAplus::plotMDS
#' @export
wgcna.plotMDS <- function(wgcna, main = NULL, scale = FALSE) {
  WGCNAplus::plotMDS(wgcna = wgcna, main = main, scale = scale)
}

#' @seealso WGCNAplus::plotFeatureUMAP
#' @export
wgcna.plotFeatureUMAP <- function(wgcna, nhub = 3, method = "clust",
                                  scale = FALSE, main = NULL,
                                  plotlib = "base", annot = NULL) {
  WGCNAplus::plotFeatureUMAP(
    wgcna = wgcna, nhub = nhub, method = method, scale = scale,
    main = main, plotlib = plotlib, annot = annot
  )
}


#' Plot module significance.
#' @seealso WGCNAplus::plotModuleSignificance
#' @export
wgcna.plotModuleSignificance <- function(wgcna, trait, main = NULL, abs = FALSE) {
  WGCNAplus::plotModuleSignificance(wgcna = wgcna, trait = trait, main = main, abs = abs)
}


#' @seealso WGCNAplus::plotConsensusSampleDendroAndColors
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
  WGCNAplus::plotConsensusSampleDendroAndColors(
    cons = cons, i = i, what = what, show.me = show.me,
    show.traits = show.traits, show.contrasts = show.contrasts,
    clust.expr = clust.expr, setLayout = setLayout, marAll = marAll,
    colorHeightMax = colorHeightMax, main = main
  )
}

#' @seealso WGCNAplus::plotSampleDendroAndColors
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
  WGCNAplus::plotSampleDendroAndColors(
    wgcna = wgcna, input.type = input.type, what = what, show.me = show.me,
    show.traits = show.traits, show.contrasts = show.contrasts,
    datTraits = datTraits, datExpr = datExpr, datME = datME,
    clust.expr = clust.expr, setLayout = setLayout, marAll = marAll,
    colorHeightMax = colorHeightMax, main = main, justdata = justdata
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
#' @seealso WGCNAplus::plotModuleHubGenes
#' @export
wgcna.plotModuleHubGenes <- function(wgcna, modules = NULL,
                                     alpha = 0.5, setpar = TRUE) {
  WGCNAplus::plotModuleHubGenes(wgcna = wgcna, modules = modules, alpha = alpha, setpar = setpar)
}

#' @seealso WGCNAplus::plotGeneNetwork
#' @export
wgcna.plotGeneNetwork <- function(wgcna, genes, col = NULL,
                                  edge.alpha = 0.3,
                                  rgamma = 4,
                                  edge.width = 6,
                                  alpha = 0.5,
                                  min.rho = 0.5,
                                  setpar = TRUE) {
  WGCNAplus::plotGeneNetwork(
    wgcna = wgcna, genes = genes, col = col, edge.alpha = edge.alpha,
    rgamma = rgamma, edge.width = edge.width, alpha = alpha,
    min.rho = min.rho, setpar = setpar
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

#' @seealso WGCNAplus::plotPowerAnalysis
#' @export
wgcna.plotPowerAnalysis <- function(datExpr, networktype = "signed",
                                    cex = 1, maxpower = 20, nmax = 2000,
                                    plots = c(
                                      "sft.modelfit", "mean.k",
                                      "dendro.IQR"
                                    ), main = NULL,
                                    RsquaredCut = 0.85, setPar = TRUE) {
  WGCNAplus::plotPowerAnalysis(
    datExpr = datExpr, networktype = networktype, cex = cex,
    maxpower = maxpower, nmax = nmax, plots = plots, main = main,
    RsquaredCut = RsquaredCut, setPar = setPar
  )
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
