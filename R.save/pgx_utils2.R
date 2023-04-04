# pgx.getMetaMatrix
get_meta_matrix <- function(pgx, methods = "meta", level = "gene") {
  fc0 <- NULL
  qv0 <- NULL
  if (level == "gene") {
    all.methods <- colnames(unclass(pgx$gx.meta$meta[[1]]$fc))
    if (is.null(methods)) methods <- all.methods
    if (any(methods %in% all.methods)) {
      methods <- intersect(methods, all.methods)
      fc0 <- sapply(pgx$gx.meta$meta, function(x) {
        rowMeans(unclass(x$fc)[, methods, drop = FALSE], na.rm = TRUE)
      })
      qv0 <- sapply(pgx$gx.meta$meta, function(x) {
        apply(unclass(x$q)[, methods, drop = FALSE], 1, max)
      }) ## maxQ
      rownames(fc0) <- rownames(qv0) <- rownames(pgx$gx.meta$meta[[1]])
    } else if (methods[1] == "meta") {
      fc0 <- sapply(pgx$gx.meta$meta, function(x) x$meta.fx)
      qv0 <- sapply(pgx$gx.meta$meta, function(x) x$meta.q)
      rownames(fc0) <- rownames(qv0) <- rownames(pgx$gx.meta$meta[[1]])
    } else {
      warning("WARNING:: pgx.getMetaFoldChangeMatrix: unknown method")
      return(NULL)
    }
  }
  if (level == "geneset") {
    all.methods <- colnames(unclass(pgx$gset.meta$meta[[1]]$fc))
    if (is.null(methods)) methods <- all.methods
    if (any(methods %in% all.methods)) {
      fc0 <- sapply(pgx$gset.meta$meta, function(x) {
        rowMeans(unclass(x$fc)[, methods, drop = FALSE], na.rm = TRUE)
      })
      qv0 <- sapply(pgx$gset.meta$meta, function(x) {
        apply(unclass(x$q)[, methods, drop = FALSE], 1, max)
      })
      rownames(fc0) <- rownames(qv0) <- rownames(pgx$gset.meta$meta[[1]])
    } else if (methods[1] == "meta") {
      fc0 <- sapply(pgx$gset.meta$meta, function(x) x$meta.fx)
      qv0 <- sapply(pgx$gset.meta$meta, function(x) x$meta.q)
      rownames(fc0) <- rownames(qv0) <- rownames(pgx$gset.meta$meta[[1]])
    } else {
      warning("WARNING:: pgx.getMetaFoldChangeMatrix: unknown method")
      return(NULL)
    }
  }
  res <- list(fc = fc0, qv = qv0)
  return(res)
}

# pgx.getMetaFoldChangeMatrix
get_meta_fold_change_matrix <- function(pgx, what = "meta", level = "gene") {
  get_meta_matrix(pgx, methods = what, level = level)
}
