##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## ----------------------------------------------------------------------
## mixOmics related functions
## ----------------------------------------------------------------------


#' @title Create Hive Plot for mixed network
#'
#' @param res Input data
#' @param ngs NGS object containing expression data
#' @param ct Contrast name or column in \code{ngs$samples} to use for node coloring
#' @param showloops Whether to show looping edges
#' @param numlab Number of node labels to show
#' @param cex Label size scaling factor
#'
#' @return Hive plot object
#'
#' @description Creates a Hive Plot to visualize a mixed gene-gene set network.
#'
#' @details  It extracts the network graph and formats it into a Hive Plot layout.
#'
#' Node importance scores are extracted from the NGS object based on the specified \code{ct} contrast.
#' These scores are used to determine node sizes in the plot.
#'
#' The number of node labels can be reduced by setting \code{numlab} to avoid overplotting.
#'
#' @export
mixHivePlot <- function(res, ngs, ct, showloops = FALSE, numlab = 6, cex = 1) {
  cat("<mixHivePlot> called\n")
  if (is.null(showloops)) showloops <- FALSE

  gr <- res$graph

  ## -------------------------------------------------------------
  ## Prepare the HivePlot data
  ## -------------------------------------------------------------
  df <- data.frame(res$edges)
  hpd <- edge2HPD(df, axis.cols = rep("grey", 3))

  hpd$edges$looping <- res$edges$looping
  hpd$edges$weight <- res$edges$importance
  loop.nodes <- unique(
    c(
      as.character(res$edges$from)[res$edges$looping],
      as.character(res$edges$to)[res$edges$looping]
    )
  )
  hpd$nodes$looping <- hpd$nodes$lab %in% loop.nodes
  hpd <- mineHPD(hpd, option = "rad <- tot.edge.count")
  hpd$nodes$degree <- hpd$nodes$radius
  hpd$edges$from <- hpd$nodes$lab[hpd$edges$id1]
  hpd$edges$to <- hpd$nodes$lab[hpd$edges$id2]

  if (is.null(ct) || is.null(ngs)) {
    fx <- tapply(igraph::V(gr)$importance, sub(".*:", "", igraph::V(gr)$name), max)
  } else if (ct %in% names(ngs$gx.meta$meta)) {
    ## use fold change as radial layout
    fc <- ngs$gx.meta$meta[[ct]]$meta.fx
    names(fc) <- rownames(ngs$gx.meta$meta[[ct]])
    gs <- ngs$gset.meta$meta[[ct]]$meta.fx
    names(gs) <- rownames(ngs$gset.meta$meta[[ct]])
    fc <- fc / max(abs(fc), na.rm = TRUE)
    gs <- gs / max(abs(gs), na.rm = TRUE)
    fx <- c(fc, gs)
  } else if (ct %in% colnames(ngs$samples)) {
    group <- ngs$samples[, ct]

    design <- stats::model.matrix(~ 0 + group)
    colnames(design) <- sub("group", "", colnames(design))
    fit <- limma::eBayes(limma::lmFit(ngs$X, design))
    stat <- limma::topTable(fit, number = Inf)
    fx <- stat$F
    names(fx) <- rownames(ngs$X)
  } else {
    stop("FATAL:: mixHivePlot: unknown contrast/conditions=", ct, "\n")
  }

  fx <- fx / max(abs(fx), na.rm = TRUE)
  g <- sub("[1-9]:", "", hpd$nodes$lab)
  hpd$nodes$radius <- rank(fx[g], na.last = "keep")
  hpd$nodes$radius <- 100 * hpd$nodes$radius / max(hpd$nodes$radius, na.rm = TRUE)

  maxgrp <- unlist(lapply(res$W, function(w) max.col(w)))

  names(maxgrp) <- as.vector(sapply(res$W, rownames))

  ## use importance as node size
  importance <- igraph::V(gr)$importance
  names(importance) <- igraph::V(gr)$name

  hpd$nodes$size <- abs(importance[hpd$nodes$lab])
  hpd$nodes$size <- 1.6 * (hpd$nodes$size / max(hpd$nodes$size))**0.5
  hpd$nodes$axis <- as.integer(sub(":.*", "", hpd$nodes$lab))
  hpd$nodes$color <- c("red3", "blue2")[1 + 1 * (fx[g] > 0)]
  wt1 <- hpd$edges$weight ## edge.importance
  wt1 <- rank(abs(wt1), na.last = "keep") * sign(wt1)
  hpd$edges$weight <- 3 * abs(wt1 / max(abs(wt1)))**2
  hpd$edges$color <- psych::alpha("grey70", 0.3)

  jj <- which(hpd$edges$looping)
  if (showloops && length(jj)) {
    hpd$edges[jj, ]
    hpd$edges$color <- psych::alpha("grey70", 0.2)
    hpd$edges$color[jj] <- psych::alpha("red3", 0.3)
  }

  axis.names <- names(res$X)
  makeAcronym <- function(x) {
    x <- gsub("[)(]", "", x)
    sapply(strsplit(x, split = "[_ -]"), function(s) {
      if (length(s) == 1) {
        return(substring(s, 1, 2))
      }
      toupper(paste(substring(s, 1, 1), collapse = ""))
    })
  }
  axis.names <- sapply(names(res$X), makeAcronym) ## see pgx-functions

  ## -------------------------------------------------------------
  ## Finally do the plotting
  ## -------------------------------------------------------------
  hpd$nodes$size <- cex * hpd$nodes$size
  hpd$edges$weight <- cex * hpd$edges$weight
  mr <- max(hpd$nodes$radius)
  plotHive(hpd,
    ch = 5, bkgnd = "white",
    axLabs = axis.names,
    axLab.pos = c(1, 1.2, 1.2) * 0.15 * mr, #

    axLab.gpar = grid::gpar(
      col = "black", fontsize = 18 * cex,
      lwd = 4, fontface = "bold"
    )
  )

  tt <- paste("edge.width = edge.importance",
    "node.size = variable importance",
    "axis = fold-change or F-stat",
    sep = "\n"
  )
  grid::grid.text(tt,
    x = 0.2 * mr, y = -1.0 * mr, default.units = "native",
    just = "left", gp = grid::gpar(fontsize = 9, col = "black")
  )

  ## axis 1
  rot.xy <- function(x, y, deg) {
    a <- deg * pi / 180
    rx <- cos(a) * x - sin(a) * y
    ry <- sin(a) * x + cos(a) * y
    cbind(rx, ry)
  }
  rot <- c(0, 120, 240)
  mr <- max(hpd$nodes$radius)
  yoff <- c(0, -0, +0)

  k <- 1
  for (k in 1:3) {
    kk <- which(hpd$nodes$axis == k)
    if (showloops) {
      kk <- which(hpd$nodes$axis == k & hpd$nodes$looping)
    }
    jj <- Matrix::head(kk[order(-hpd$nodes$size[kk])], numlab) ## number of labels
    rr <- hpd$nodes$radius
    rx <- rot.xy(0, rr[jj] + 5, rot[k])


    lab <- sub(".*:", "", hpd$nodes$lab[jj])
    ##    pt <- maptools::pointLabel(rx[, 1], rx[, 2], labels = lab, cex = cex * 2, doPlot = FALSE)
    ##    px <- cbind(pt$x, pt$y)
    px <- cbind(rx[, 1], rx[, 2])
    px[, 1] <- px[, 1] + 4
    grid::grid.text(lab,
      ## x = 10 + rx[,1], y = rx[,2],
      x = px[, 1], y = px[, 2],
      default.units = "native", just = "left",
      gp = grid::gpar(fontsize = 12 * cex, col = "black")
    )

    grid::grid.segments(rx[, 1], rx[, 2], px[, 1] - 1, px[, 2],
      default.units = "native"
    )
  }
}

#' @describeIn mixHivePlot function generates a hive plot visualization of variable loadings
#' from a lmer model result object.
#' @export
mixPlotLoadings <- function(res, showloops = FALSE, cex = 1) {
  cat("<mixPlotLoadings> called\n")
  levels <- levels(res$Y)
  ny <- length(levels)
  klrpal <- c("blue2", "orange2")
  klrpal <- rep(RColorBrewer::brewer.pal(n = 8, "Set2"), 10)[1:ny]

  names(klrpal) <- levels

  plotly::layout(matrix(1:6, 1, 6), widths = c(1, 0.5, 1, 0.5, 1, 0.5))
  k <- 1
  for (k in 1:3) {
    W <- res$W[[k]]
    graphics::par(mar = c(5, 8 * cex, 4, 0), mgp = c(2.2, 0.8, 0))
    graphics::barplot(t(W),
      horiz = TRUE, las = 1,
      border = NA, col = klrpal,
      names.arg = sub(".*:", "", rownames(W)),
      xlim = c(0, 1.1) * max(rowSums(W, na.rm = TRUE)),
      cex.axis = 1 * cex, cex.names = 1.1 * cex,
      cex.lab = 1 * cex, xlab = "importance"
    )
    graphics::title(names(res$loadings)[k],
      cex.main = 1.3 * cex,
      adj = 0.33, xpd = NA
    )
    graphics::legend("topright",
      legend = names(klrpal),
      cex = 1.1 * cex, pch = 15, col = klrpal, #
      y.intersp = 0.85, inset = c(0.15, 0.03)
    )

    if (k < 99) {
      ## add correlation lines
      graphics::par(mar = c(5, 0, 4, 0))

      plot(0,
        type = "n", xlim = c(0, 1), ylim = c(0, nrow(W)),
        xaxt = "n", yaxt = "n", bty = "n", xlab = ""
      )

      g1 <- rownames(res$W[[k]])
      g2 <- rownames(res$W[[ifelse(k < 3, k + 1, 1)]])

      sel <- which(res$edges[, "from"] %in% c(g1, g2) &
        res$edges[, "to"] %in% c(g1, g2))
      sel
      ee <- res$edges[sel, ]
      ii <- apply(ee[, 1:2], 1, function(e) which(e %in% g1))
      jj <- apply(ee[, 1:2], 1, function(e) which(e %in% g2))
      ee$from <- res$edges[sel, ][cbind(1:nrow(ee), ii)]
      ee$to <- res$edges[sel, ][cbind(1:nrow(ee), jj)]


      lwd <- ee$importance

      lwd <- rank(abs(lwd), na.last = "keep")**1.5
      lwd <- 3.0 * cex * (lwd / max(lwd))
      lty <- 1 + 1 * (sign(ee$rho) < 0)
      xy <- cbind(match(ee$from, g1), match(ee$to, g2))
      xy[, 2] <- (xy[, 2] - 0.5) / length(g2) * length(g1)
      klr <- rep(psych::alpha("grey70", 0.3), nrow(ee))
      if (showloops) {
        klr <- rep(psych::alpha("grey70", 0.2), nrow(ee))
        klr[which(ee$looping)] <- psych::alpha("red3", 0.3)
      }
      graphics::segments(0, xy[, 1] - 0.5, 1, xy[, 2], lwd = lwd, col = klr, lty = lty)
      rr <- paste(round(range(abs(ee$rho)), 2), collapse = ",")

      graphics::title(sub = paste0("[", rr, "]"), line = -1.2, cex.sub = cex)
    }
  }
}


## ----------------------------------------------------------------------
## Variable importance functions
## ----------------------------------------------------------------------


#' @title Variable importance for survival models
#'
#' @param X Numeric matrix of predictor variables, rows are variables
#' @param time Vector of event/censoring times
#' @param status Vector indicating event (1) or censoring (0)
#' @param methods Methods for computing variable importance. Options are "glmnet", "randomforest", "boruta", "xgboost", "pls".
#'
#' @return Named list with variable importance scores by method.
#'
#' @description Compute variable importance scores for predictors in a survival analysis model using different methods.
#'
#' @details This function calculates variable importance scores using a variety of methods suitable for survival analysis.
#' The input data consists of a predictor matrix \code{X}, a vector of event/censoring times \code{time} and a status indicator vector \code{status}.
#'
#' The following methods can be selected via the \code{methods} parameter:
#' \itemize{
#' \item glmnet: Absolute value of coefficients from elastic net Cox model
#' \item randomForest: Variable importance from random survival forest
#' \item boruta: Boruta variable selection algorithm
#' \item xgboost: Importance scores from XGBoost survival model
#' \item pls: Absolute coefficients from partial least squares Cox model
#' }
#'
#' Variable importance scores are returned for each method in a named list.
#' These scores can be used to select important predictors for survival modeling.
#'
#' @export
pgx.survivalVariableImportance <- function(X, time, status,
                                           methods = c("glmnet", "randomforest", "boruta", "xgboost", "pls")) {
  ## ----------------------------------------------------------------------
  ## multi-class version
  ## ----------------------------------------------------------------------

  imp <- list()
  xnames <- rownames(X)
  sdx <- apply(X, 1, stats::sd, na.rm = TRUE)
  if (nrow(X) == 1) X <- rbind(X, X)

  if (!inherits(status, "logical") && all(status %in% c(0, 1, NA))) {
    stop("status must be logical or 0/1")
  }

  y <- survival::Surv(time + 0.0001, status)

  if ("glmnet" %in% methods) {
    fam <- "cox"

    NFOLD <- 5
    out0 <- glmnet::cv.glmnet(t(X), y, alpha = 0, family = fam, standardize = TRUE, nfolds = NFOLD)
    cf0 <- glmnet::coef.glmnet(out0, s = "lambda.min")[, 1]
    out1 <- glmnet::cv.glmnet(t(X), y, alpha = 1, family = fam, standardize = TRUE, nfolds = NFOLD)
    cf1 <- glmnet::coef.glmnet(out1, s = "lambda.min")[, 1]
    out0a <- glmnet::cv.glmnet(t(X), y, alpha = 0, family = fam, standardize = FALSE, nfolds = NFOLD)
    cf0a <- glmnet::coef.glmnet(out0a, s = "lambda.min")[, 1]
    out1a <- glmnet::cv.glmnet(t(X), y, alpha = 1, family = fam, standardize = FALSE, nfolds = NFOLD)
    cf1a <- glmnet::coef.glmnet(out1a, s = "lambda.min")[, 1]
    imp[["coxnet.a0"]] <- (cf0 / max(abs(cf0)) + cf0a / max(abs(cf0a))) * sdx[names(cf0)]
    imp[["coxnet.a1"]] <- (cf1 / max(abs(cf1)) + cf1a / max(abs(cf1a))) * sdx[names(cf1)]
  }

  if ("randomforest" %in% methods) {
    df <- data.frame(time = time, status = status, t(X))
    fit_rf <- randomForestSRC::rfsrc(survival::Surv(time, status) ~ ., data = df)
    vimp <- randomForestSRC::vimp(fit_rf)$importance
    imp[["randomForest"]] <- vimp
  }

  if ("boruta" %in% methods) {
    imp4 <- rep(0, nrow(X))
    niter <- 4
    for (k in 1:niter) {
      jj <- sample(ncol(X), ncol(X) * 0.9)
      out3 <- Boruta(t(X[, jj, drop = FALSE]), y[jj])
      fd <- factor(out3$finalDecision, levels = c("Rejected", "Tentative", "Confirmed"))
      fd <- (as.integer(fd) - 1) / 2
      imp4 <- imp4 + fd / niter
    }

    names(imp4) <- rownames(X)
    imp[["Boruta"]] <- imp4
  }

  if ("xgboost" %in% methods) {
    yy <- ifelse(!status, -time, time)
    bst <- xgboost::xgboost(
      data = t(X), label = yy, booster = "gbtree",
      max_depth = 2, eta = 1, nthread = 2, nrounds = 2,
      verbose = 0, objective = "survival:cox"
    )

    xgmat <- xgboost::xgb.importance(model = bst)
    imp5 <- xgmat$Gain**0.33
    imp5 <- imp5[match(rownames(X), xgmat$Feature)]
    names(imp5) <- rownames(X)
    imp5[which(is.na(imp5))] <- 0
    Matrix::tail(sort(imp5))
    imp[["xgboost"]] <- imp5

    ## linear model
    bst2 <- xgboost::xgboost(
      data = t(X), label = yy, booster = "gblinear",
      max_depth = 2, eta = 1, nthread = 2, nrounds = 2,
      verbose = 0, objective = "survival:cox"
    )
    xgmat <- xgboost::xgb.importance(model = bst2)
    imp6 <- xgmat$Weight
    names(imp6) <- xgmat$Feature
    imp6 <- imp6[match(rownames(X), names(imp6))]
    Matrix::tail(sort(imp6))
    imp[["xgboost.lin"]] <- imp6
  }

  if ("pls" %in% methods) {
    res <- plsRcox::plsRcox(t(X), time = time, event = status, nt = 5)
    cf <- res$Coeffs[, 1]
    cf[is.na(cf)] <- 0
    cf <- cf * sdx[names(cf)] ## really?
    imp[["pls.cox"]] <- abs(cf) / max(abs(cf), na.rm = TRUE)
  }

  P <- do.call(cbind, imp)
  P <- abs(P) ## always positive??
  P[is.na(P)] <- 0
  P <- P[xnames, , drop = FALSE]
  return(P)
}


CARET.METHODS <- c("svmLinear", "rpart", "glm", "glmStepAIC", "pda", "sparseLDA", "pls", "plsda")
CARET.METHODS <- c("pls", "plsda")

## for compatibility before v3.3.0
#' @export
pgx.multiclassVariableImportance <- function(X, y, methods) {
  res <- pgx.variableImportance(X = X, y = y, methods = methods)
  res$importance
}

#' @describeIn pgx.survivalVariableImportance Calculates variable importance scores for predictors of a multiclass response using various methods.
#' @param y Multiclass factor response variable. Contains the class labels for each sample
#' @export
pgx.variableImportance <- function(X, y,
                                   methods = c(
                                     "glmnet", "randomforest",
                                     "xgboost", "splsda", "correlation", "ftest",
                                     "boruta", CARET.METHODS
                                   )[1:6],
                                   scale = TRUE, reduce = 1000, resample = 0, add.noise = 0) {
  ## variables
  imp <- list()
  runtime <- list()
  xnames <- rownames(X)

  if (nrow(X) == 1) X <- rbind(X, X)

  ## drop missing??
  sel <- which(!is.na(y) & y != "")
  X <- X[, sel, drop = FALSE]
  y <- y[sel]

  ## convert to factor
  y <- factor(as.character(y))

  ## resample to minimum size to balance groups
  if (resample > 0) {
    NSIZE <- resample
    jj <- tapply(1:length(y), y, function(ii) {
      if (length(ii) < NSIZE) ii <- sample(c(ii, ii), NSIZE, replace = TRUE)
      return(ii)
    })
    jj <- unlist(jj)
    dbg("[pgx.variableImportance] resampling levels from: ", table(y))
    dbg("[pgx.variableImportance] resampling levels to: ", table(y[jj]))
    X <- X[, jj]
    y <- y[jj]
  }


  ## add noise
  if (add.noise > 0) {
    sdx0 <- matrixStats::rowSds(X, na.rm = TRUE)
    sdx1 <- add.noise * (0.5 * sdx0 + 0.5 * mean(sdx0, na.rm = TRUE))
    X <- X + sdx1 * matrix(rnorm(length(X)), nrow(X), ncol(X)) ## add some noise
  }

  if (scale) {
    sdx <- matrixStats::rowSds(X, na.rm = TRUE)
    X <- (X - rowMeans(X, na.rm = TRUE)) / (sdx + 1e-8)
  }

  ## reduce features using fast correlation/cov
  if (reduce > 0 && reduce < nrow(X)) {
    dy <- model.matrix(~ 0 + y)
    rho <- cov(t(X), dy) ## note: using covariance
    imp1 <- sqrt(rowMeans(rho**2))
    sel <- head(names(sort(-imp1)), reduce)
    X <- X[sel, ]
  }

  if ("correlation" %in% methods) {
    runtime[["correlation"]] <- system.time({
      dy <- model.matrix(~ 0 + y)
      rho <- cov(t(X), dy)
    })
    imp[["correlation"]] <- sqrt(rowMeans(rho**2))
  }

  if ("ftest" %in% methods) {
    runtime[["ftest"]] <- system.time({
      res <- playbase::gx.limmaF(X, y, lfc = 0, fdr = 1)
    })
    imp1 <- -log10(res[rownames(X), "P.Value"])
    names(imp1) <- rownames(X)
    imp[["ftest"]] <- imp1
  }

  if ("glmnet" %in% methods) {
    fam <- "multinomial"
    NFOLD <- 5
    runtime[["glmnet.a0"]] <- system.time({
      out0 <- glmnet::cv.glmnet(t(X), y, alpha = 0, family = fam, standardize = TRUE, nfold = NFOLD)
      cf0 <- Matrix::rowMeans(do.call(cbind, coef(out0, s = "lambda.min"))[-1, ]**2)
      out0a <- glmnet::cv.glmnet(t(X), y, alpha = 0, family = fam, standardize = FALSE, nfold = NFOLD)
      cf0a <- Matrix::rowMeans(do.call(cbind, coef(out0a, s = "lambda.min"))[-1, ]**2)
    })

    runtime[["glmnet.a1"]] <- system.time({
      out1 <- glmnet::cv.glmnet(t(X), y, alpha = 0.99, family = fam, standardize = TRUE, nfold = NFOLD)
      cf1 <- Matrix::rowMeans(do.call(cbind, coef(out1, s = "lambda.min"))[-1, ]**2)
      out1a <- glmnet::cv.glmnet(t(X), y, alpha = 0.99, family = fam, standardize = FALSE, nfold = NFOLD)
      cf1a <- Matrix::rowMeans(do.call(cbind, coef(out1a, s = "lambda.min"))[-1, ]**2)
    })

    sdx <- matrixStats::rowSds(X, na.rm = TRUE)
    imp[["glmnet.a0"]] <- (cf0 / max(abs(cf0)) + cf0a / max(abs(cf0a))) * sdx[names(cf0)]
    imp[["glmnet.a1"]] <- (cf1 / max(abs(cf1)) + cf1a / max(abs(cf1a))) * sdx[names(cf1)]
  }

  if ("randomforest" %in% methods) {
    runtime[["randomForest"]] <- system.time({
      fit_rf <- randomForest::randomForest(t(X), factor(y))
    })
    imp[["randomForest"]] <- fit_rf$importance[, 1]
  }

  if ("boruta" %in% methods) {
    imp4 <- rep(0, nrow(X))
    niter <- 4
    runtime[["Boruta"]] <- system.time({
      for (k in 1:niter) {
        jj <- sample(ncol(X), ncol(X) * 0.9)
        out3 <- Boruta::Boruta(t(X[, jj, drop = FALSE]), y[jj])
        fd <- factor(out3$finalDecision,
          levels = c("Rejected", "Tentative", "Confirmed")
        )
        fd <- (as.integer(fd) - 1) / 2
        imp4 <- imp4 + fd / niter
      }
    })
    names(imp4) <- rownames(X)
    imp[["Boruta"]] <- imp4
  }

  if ("xgboost" %in% methods) {
    ny <- length(table(y))
    yy <- as.integer(factor(y)) - 1
    runtime[["xgboost"]] <- system.time({
      bst <- xgboost::xgboost(
        data = t(X), label = yy, booster = "gbtree",
        max_depth = 2, eta = 1, nthread = 2, nrounds = 2,
        num_class = ny,
        verbose = 0, objective = "multi:softmax"
      )
    })

    xgmat <- xgboost::xgb.importance(model = bst)
    imp5 <- xgmat$Gain**0.2
    imp5 <- imp5[match(rownames(X), xgmat$Feature)]
    names(imp5) <- rownames(X)
    imp[["xgboost"]] <- imp5

    ## linear model
    runtime[["xgboost.lin"]] <- system.time({
      bst2 <- xgboost::xgboost(
        data = t(X), label = yy, booster = "gblinear",
        max_depth = 2, eta = 1, nthread = 2, nrounds = 2,
        num_class = ny,
        verbose = 0, objective = "multi:softmax"
      )
    })
    xgmat <- xgboost::xgb.importance(model = bst2)
    imp6 <- xgmat$Weight
    names(imp6) <- xgmat$Feature
    imp6 <- imp6[match(rownames(X), names(imp6))]
    imp[["xgboost.lin"]] <- imp6
  }

  if ("splsda" %in% methods) {
    n <- min(25, nrow(X))
    colnames(X) <- names(y) <- paste0("sample", 1:length(y))
    runtime[["splsda"]] <- system.time({
      res <- mixOmics::splsda(t(X), y, keepX = c(n, n))
    })
    impx <- rowMeans(res$loadings$X[rownames(X), ]**2)
    imp[["splsda"]] <- impx
  }

  #  CARET.METHODS = c("glmnet",'svmLinear','rpart','xgbLinear','glm','glmStepAIC',
  #  'pda','sparseLDA','pls','plsda')
  caret.methods <- intersect(methods, CARET.METHODS)
  if (length(caret.methods)) {
    repeatedcv <- caret::trainControl(method = "repeatedcv", number = 25, repeats = 1)
    m <- "pls"
    for (m in caret.methods) {
      message("computing variable importance for ", m)
      t0 <- NULL
      if (m == "plsda") {
        runtime[["plsda"]] <- system.time({
          fit <- try(caret::plsda(x = t(X), y = y, ncomp = min(20, nrow(X) / 3)))
        })
        if (!"try-error" %in% class(fit)) {
          imp[[m]] <- sqrt(rowMeans(caret::varImp(fit)**2))
        }
      } else {
        runtime[[m]] <- system.time({
          fit <- try(caret::train(
            x = t(X), y = y, method = m,
            trControl = repeatedcv, trace = FALSE
          ))
        })
        if (!"try-error" %in% class(fit)) {
          imp[[m]] <- sqrt(rowMeans(caret::varImp(fit)$importance**2))
        }
      }
    } ## for caret.methods
  }
  names(imp)

  imp <- imp[!sapply(imp, is.null)]
  imp <- lapply(imp, function(v) v[match(xnames, names(v))])
  P <- do.call(cbind, imp)
  colnames(P) <- names(imp)
  rownames(P) <- xnames
  P <- abs(P) ## always positive??
  P[is.na(P)] <- 0
  ## P <- P[xnames, , drop = FALSE]

  ## normalize values by ranking
  ii <- which(P == 0)
  P <- 1 / apply(-P, 2, rank)
  P[ii] <- 0

  ## parse runtimes
  tt <- unlist(sapply(runtime[colnames(P)], function(t) t[1]))
  names(tt) <- sub(".user.self", "", names(tt))
  tt <- tt[match(colnames(P), names(tt))]
  names(tt) <- colnames(P)
  tt
  sort(tt)

  res <- list(importance = P, runtime = tt)
  return(res)
}


#' @export
plotImportance <- function(P, p.sign = NULL, top = 50, runtime = NULL) {
  P[is.na(P)] <- 0
  if (is.null(p.sign)) p.sign <- rep(1, nrow(P))
  p.sign <- sign(p.sign)
  P1 <- 1 / apply(-abs(P), 2, rank) * p.sign
  P1 <- P1[order(-rowMeans(P1, na.rm = TRUE)), ]
  if (any(p.sign < 0)) {
    sel <- unique(head(rownames(P1), top / 2), tail(rownames(P1), top / 2))
  } else {
    sel <- head(rownames(P1), top)
  }
  P1 <- P1[sel, ]
  P1 <- P1[order(rowMeans(P1, na.rm = TRUE)), ]

  par(mfrow = c(2, 2), mar = c(5, 15, 4, 4))
  if (is.null(runtime)) par(mfrow = c(1, 2))
  barplot(t(P1),
    las = 1, horiz = TRUE,
    xlab = "cumulative importance", main = "variable importance"
  )
  legend("bottomright", legend = colnames(P1), fill = grey.colors(ncol(P1)))

  par(mar = c(8, 4, 4, 15))
  playbase::gx.imagemap(scale(P1, center = FALSE))

  if (!is.null(runtime)) {
    par(mar = c(5, 15, 4, 4))
    runtime <- sort(runtime, decreasing = TRUE)
    runtime <- pmax(runtime, 0.1 * min(runtime[runtime > 0], na.rm = TRUE))
    log.runtime <- log10(runtime)
    log.runtime <- log.runtime - min(log10(runtime))
    barplot(log.runtime,
      las = 1, horiz = TRUE,
      xlab = "time (log10)", main = "algorithm runtime"
    )
  }
}






## =====================================================================================
## =========================== END OF FILE =============================================
## =====================================================================================
