##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

## ---------------------------------------------------------------
## ------------- Functions for WordCloud ------------------------
## ---------------------------------------------------------------


#' @title Calculate word frequencies for word cloud
#'
#' @param pgx A PGX object
#' @param progress A progress object to track progress.
#' @param pg.unit Increment value for progress tracking.
#'
#' @return A matrix with word frequencies for word cloud visualization.
#'
#' @description Calculates word frequencies from gene set names for generating a word cloud.
#'
#' @details This function takes a NGStest object containing gene set analysis results.
#' It extracts the gene set names, filters out common words, and calculates word frequencies.
#'
#' Gene set names are split into words and filtered to remove common stopwords. The top 1000 most frequent
#' words are kept. Word frequencies are calculated based on occurrence in the gene set names.
#'
#' The output is a matrix of word frequencies suitable for generating a word cloud visualization.
#'
#' @export
pgx.calculateWordCloud <- function(pgx, progress = NULL, pg.unit = 1) {
  if (is.null(pgx$gset.meta)) {
    message("[pgx.calculateWordCloud] ERROR: no gset.meta in object")
    return(NULL)
  }

  if (!is.null(progress)) progress$set(message = "WordCloud", value = 0)

  ## get gset meta foldchange-matrix
  S <- sapply(pgx$gset.meta$meta, function(x) x$meta.fx)
  rownames(S) <- rownames(pgx$gset.meta$meta[[1]])
  S[is.na(S)] <- 0
  S <- S[order(-rowMeans(S**2)), , drop = FALSE]
  dim(S)

  ## exclude down, GSE gene sets??????
  S <- S[grep("dn|down|^gse", rownames(S), ignore.case = TRUE, invert = TRUE), , drop = FALSE]

  if (nrow(S) <= 10 || NCOL(S) == 0) {
    message("[pgx.calculateWordCloud] WARNING:: not enough valid genesets left")
    return(NULL)
  }

  if (!is.null(progress)) progress$inc(0.2 * pg.unit, detail = "calculating word frequencies")

  ## Determine top most frequent terms
  sname <- gsub(".*:", "", tolower(rownames(S)))
  sname <- gsub("b.cell", "bcell", sname)
  sname <- gsub("t.cell", "tcell", sname)
  words <- strsplit(sname, split = "[-_ ]")
  names(words) <- rownames(S)
  terms <- names(sort(-table(unlist(words))))
  stopwords <- strsplit(
    "with int like strand tconv pid lee and or mouse dn small big human homo sapiens mus musculus drug early late hsa gse culture in the a g for line r up down events anti large targets tissue vitro process cells ctrl regulation processing common pathway months days white pre post process all mice from",
    split = " "
  )[[1]]
  terms <- terms[which(!terms %in% stopwords)]
  terms <- terms[sapply(terms, nchar) > 2]
  terms <- grep("[0-9]|^\\(", terms, invert = TRUE, value = TRUE)
  terms <- Matrix::head(terms, 1000)

  ## Calculate incidence matrix
  words2 <- lapply(words, function(w) intersect(w, terms))
  words2 <- words2[sapply(words2, length) > 0]
  idx <- lapply(1:length(words2), function(i) cbind(i, match(words2[[i]], terms)))
  idx <- do.call(rbind, idx)

  W <- Matrix::sparseMatrix(idx[, 1], idx[, 2], x = 1)
  rownames(W) <- names(words2)
  colnames(W) <- terms
  dim(W)

  ## filter on minimal size and maximum ratio???
  nn <- Matrix::colSums(W, na.rm = TRUE)
  nr <- nn / nrow(W)
  W <- W[, which(nn >= 3 & nr <= 0.5), drop = FALSE]
  dim(W)

  if (ncol(W) < 1) {
    message("[pgx.calculateWordCloud] WARNING:: no valid words left")
    return(NULL)
  }

  ## align geneset expression matrix
  S <- S[rownames(W), , drop = FALSE]

  if (!is.null(progress)) progress$inc(0.3 * pg.unit, detail = "computing GSEA")

  ## compute for average contrast
  rms.FC <- Matrix::rowMeans(S**2)**0.5
  rms.FC <- rms.FC + 0.01 * stats::rnorm(length(rms.FC))
  gmt <- apply(W, 2, function(x) names(which(x != 0)))
  suppressWarnings(res <- fgsea::fgseaSimple(gmt, rms.FC, nperm = 1000))
  res$leadingEdge <- sapply(res$leadingEdge, paste, collapse = "//")
  colnames(res)[1] <- "word"

  ## --------- only significant and positive
  res <- res[(res$padj < 1 & res$NES > 0), ]
  res <- res[order(-abs(res$NES)), ]

  ## now compute significant terms for all contrasts
  all.gsea <- list()
  for (i in 1:ncol(S)) {
    fc <- as.vector(S[, i])
    names(fc) <- rownames(W)
    fc <- fc + 0.01 * stats::rnorm(length(fc))
    gmt1 <- gmt[as.character(res$word)]
    res1 <- fgsea::fgseaSimple(gmt1, fc, nperm = 1000)
    res1$leadingEdge <- sapply(res1$leadingEdge, paste, collapse = "//")

    colnames(res1)[1] <- "word"
    all.gsea[[colnames(S)[i]]] <- res1
  }
  all.gsea[["rms.FC"]] <- res

  if (!is.null(progress)) progress$inc(0.25 * pg.unit, detail = "clustering")

  if (NCOL(W) <= 3) {
    ## t-SNE doesn't like 1-2 columns...
    W <- cbind(W, W, W, W, W)
    W <- W + 1e-2 * matrix(stats::rnorm(length(W)), nrow(W), ncol(W))
  }
  nb <- floor(min(max(ncol(W) / 4, 1), 10))
  nb
  message("[pgx.calculateWordCloud] dim(W) = ", paste(dim(W), collapse = "x"))
  message("[pgx.calculateWordCloud] setting perplexity = ", nb)

  pos1 <- try(Rtsne::Rtsne(t(as.matrix(W)), perplexity = nb, check_duplicates = FALSE)$Y)
  class(pos1)
  if ("try-error" %in% class(pos1)) {
    pos1 <- try(Rtsne::Rtsne(t(as.matrix(W)),
      perplexity = ceiling(nb / 2),
      check_duplicates = FALSE
    )$Y)
  }
  pos2 <- uwot::umap(t(as.matrix(W)), n_neighbors = max(nb, 2))
  rownames(pos1) <- rownames(pos2) <- colnames(W)
  colnames(pos1) <- colnames(pos2) <- c("x", "y")
  pos1 <- pos1[match(res$word, rownames(pos1)), ]
  pos2 <- pos2[match(res$word, rownames(pos2)), ]

  # sometimes we have words that NA is tsne, make sure we remove them
  # (likely special characters) in windows or wsl
  pos1 <- pos1[!is.na(rownames(pos1)), ]
  pos2 <- pos2[!is.na(rownames(pos2)), ]

  ordered_words <- all.gsea[[1]]$word
  res$tsne <- res$tsne[ordered_words, ]
  res$umap <- res$umap[ordered_words, ]

  all.res <- list(gsea = all.gsea, S = S, W = W, tsne = pos1, umap = pos2)
  return(all.res)
}



#' @title Plot wordcloud graphic from pgx
#'
#' @param pgx A NGStest object containing gene set results.
#' @param progress A progress object to track progress.
#' @param pg.unit Increment value for progress tracking.
#'
#' @return A wordcloud plot
#'
#' @description Plot wordcloud graphic from pgx
#'
#' @export
pgx.plotWordCloud <- function(pgx, contrast) {
  if (!"wordcloud" %in% names(pgx)) {
    message("[pgx.plotWordCloud] ERROR: pgx object has no wordcloud results")
    return()
  }

  res <- pgx$wordcloud
  gsea1 <- res$gsea[[contrast]]

  ## sometimes we have words that NA is tsne, make sure we remove
  ## them (likely special characters) in windows or wsl
  res$tsne <- res$tsne[!is.na(rownames(res$tsne)), ]
  res$umap <- res$umap[!is.na(rownames(res$umap)), ]

  ordered_words <- gsea1$word
  res$tsne <- res$tsne[ordered_words, ]
  res$umap <- res$umap[ordered_words, ]

  ## end sometimes we have words that NA is tsne, make sure we
  ## remove them (likely special characters) in windows or wsl
  df <- data.frame(gsea1, tsne = res$tsne, umap = res$umap)
  df <- df[order(-df$NES), ]

  cex1 <- 1 + round((5 * rank(abs(df$NES)) / nrow(df))**2)
  cex2 <- (-log10(df$padj))**1.0
  size <- 10 * abs(cex1 * cex2)**1
  size[is.na(size)] <- 0
  minsize <- tail(sort(size), 250)[1]

  color.pal <- c("Blues", "Greys", "Accent", "Dark2")[1]

  par(mar = c(1, 1, 1, 1) * 0)
  suppressWarnings(suppressMessages(
    wordcloud::wordcloud(
      words = df$word, freq = size,
      colors = RColorBrewer::brewer.pal(8, color.pal),
      scale = c(2, 0.1) * 0.9, min.freq = minsize
    )
  ))
}
