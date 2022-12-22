#makeFullContrasts
make_full_contrasts <- function(labels, by.sample=FALSE) {
  levels <- sort(unique(as.character(labels)))
  cc <- t(combn(levels,2))
  contr.matrix <- c()
  for(i in nrow(cc):1) {
    ctr <- 1*(levels==cc[i,1]) - 1*(levels==cc[i,2])
    contr.matrix <- cbind(ctr,contr.matrix)
    colnames(contr.matrix)[1] <- paste(cc[i,],collapse="_vs_")
  }
  rownames(contr.matrix) <- levels
  if(by.sample) {
    design <- model.matrix( ~ 0 + labels )
    colnames(design) <- sub("^labels","",colnames(design))
    rownames(design) <- names(labels)
    design <- design[,rownames(contr.matrix)]
    contr.matrix <- design %*% contr.matrix
  }
  return(contr.matrix)
}