#' Test for the globaltest method of gset.fitContrastsWithAllMethods
set.seed(1746)
X <- matrix(rnorm(60 * 12), 60, 12,
  dimnames = list(paste0("m", 1:60), paste0("s", 1:12))
)
X[1:6, 7:12] <- X[1:6, 7:12] + 2 ## "up" is up in the second group
gmt <- list(up = paste0("m", 1:6), null = paste0("m", 31:36), none = c("x1", "x2"))
G <- Matrix::sparseMatrix(
  i = match(unlist(gmt[1:2]), rownames(X)), j = rep(1:2, each = 6), x = 1,
  dims = c(60, 3), dimnames = list(rownames(X), names(gmt))
)
contr.matrix <- matrix(rep(c(-1, 1), each = 6), ncol = 1,
  dimnames = list(colnames(X), "B_vs_A")
)

## max() over the all-NA row of "none" warns; that is the meta step's, not ours
res <- suppressWarnings(playbase::gset.fitContrastsWithAllMethods(
  gmt, X, G, design = NULL, contr.matrix = contr.matrix, methods = "globaltest"
))
out <- res$meta[["B_vs_A"]]

test_that("globaltest is run and reported", {
  expect_equal(colnames(unclass(out$p)), "globaltest")
  expect_equal(rownames(out), names(gmt))
  expect_true(is.na(unclass(out$p)[3, 1])) ## "none" has no measured member
})

test_that("globaltest finds the shifted set, signed up, and not the null one", {
  expect_lt(out["up", "meta.p"], 1e-3)
  expect_gt(out["up", "meta.fx"], 0)
  expect_gt(out["null", "meta.p"], 0.05)
})
