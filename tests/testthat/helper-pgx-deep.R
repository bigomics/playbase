## Helper for test-pgx-deep.R: makes MultiOmicsSAE available (from the
## playbase namespace if present, else by sourcing the needed R files) and
## builds a tiny synthetic multi-omics fixture.

if (requireNamespace("torch", quietly = TRUE) && requireNamespace("R6", quietly = TRUE)) {
  library(R6)
  library(torch)
  library(matrixStats)
  if (!"playbase" %in% loadedNamespaces() || !exists("MultiOmicsSAE", asNamespace("playbase"))) {
    for (f in c("pgx-feature-utils.R", "pgx-deep.R")) {
      sys.source(file.path("..", "..", "R", f), envir = globalenv())
    }
  }
}

## 3 views x 30 features x 24 samples, 3 balanced classes with a separable
## signal in the first 10 features of every view.
make_sae_fixture <- function(nfeat = 30, nsamp = 24, nclass = 3) {
  set.seed(1)
  torch::torch_manual_seed(1)
  cls <- factor(rep(paste0("class", seq_len(nclass)), each = nsamp / nclass))
  samples <- paste0("s", seq_len(nsamp))
  X <- lapply(c(gx = "g", px = "p", mx = "m"), function(prefix) {
    x <- matrix(rnorm(nfeat * nsamp), nfeat, nsamp)
    for (j in seq_len(nclass)) {
      x[1:10, cls == levels(cls)[j]] <- x[1:10, cls == levels(cls)[j]] + 3 * (j - 2)
    }
    dimnames(x) <- list(paste0(prefix, seq_len(nfeat)), samples)
    x
  })
  names(cls) <- samples
  list(X = X, Y = list(class = cls))
}
