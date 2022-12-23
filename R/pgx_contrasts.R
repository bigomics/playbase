# makeFullContrasts
make_full_contrasts <- function(labels, by.sample = FALSE) {
  levels <- sort(unique(as.character(labels)))
  cc <- t(combn(levels, 2))
  contr.matrix <- c()
  for (i in nrow(cc):1) {
    ctr <- 1 * (levels == cc[i, 1]) - 1 * (levels == cc[i, 2])
    contr.matrix <- cbind(ctr, contr.matrix)
    colnames(contr.matrix)[1] <- paste(cc[i, ], collapse = "_vs_")
  }
  rownames(contr.matrix) <- levels
  if (by.sample) {
    design <- model.matrix(~ 0 + labels)
    colnames(design) <- sub("^labels", "", colnames(design))
    rownames(design) <- names(labels)
    design <- design[, rownames(contr.matrix)]
    contr.matrix <- design %*% contr.matrix
  }
  return(contr.matrix)
}

# makeContrastsFromLabelMatrix
make_contrasts_from_label_matrix <- function(lab.matrix) {
  ct.names <- colnames(lab.matrix)
  main.grp <- sapply(strsplit(ct.names, split = "_vs_"), "[", 1)
  ctrl.grp <- sapply(strsplit(ct.names, split = "_vs_"), "[", 2)
  main.grp <- sub(".*:", "", main.grp)
  ctrl.grp <- sub("@.*", "", ctrl.grp)
  contr.mat <- matrix(0, nrow(lab.matrix), ncol(lab.matrix))
  rownames(contr.mat) <- rownames(lab.matrix)
  colnames(contr.mat) <- colnames(lab.matrix)
  for (i in 1:ncol(lab.matrix)) {
    lab1 <- trimws(lab.matrix[, i])
    j1 <- which(lab1 == main.grp[i])
    j0 <- which(lab1 == ctrl.grp[i])
    contr.mat[j1, i] <- +1 / length(j1)
    contr.mat[j0, i] <- -1 / length(j0)
  }
  contr.mat
}
