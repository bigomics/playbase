#' Test for the globaltest method of gset.fitContrastsWithAllMethods
set.seed(1746)
X <- matrix(rnorm(60 * 12), 60, 12,
  dimnames = list(paste0("m", 1:60), paste0("s", 1:12))
)
X[c(1:6, 13:15), 7:12] <- X[c(1:6, 13:15), 7:12] + 2 ## up in the second group
X[c(7:12, 16:18), 7:12] <- X[c(7:12, 16:18), 7:12] - 2 ## down
gmt <- list(
  up = paste0("m", 1:6),
  "DB:down y (z)" = paste0("m", 7:12),
  mixed = paste0("m", 13:18),
  null = paste0("m", 31:36),
  small = paste0("m", 1:4),
  small.null = paste0("m", 41:44),
  none = c("x1", "x2")
)
G <- playbase::gmt2mat(gmt, bg = rownames(X))
contr.matrix <- matrix(rep(c(-1, 1), each = 6),
  ncol = 1,
  dimnames = list(colnames(X), "B_vs_A")
)
fit <- function(methods, ...) {
  res <- playbase::gset.fitContrastsWithAllMethods(
    gmt, X, G,
    design = NULL, contr.matrix = contr.matrix, methods = methods, ...
  )
  out <- res$meta[["B_vs_A"]]
  rownames(out$p) <- rownames(out) ## dropped to save memory
  out
}

out <- fit("globaltest")

test_that("globaltest is run and reported", {
  expect_equal(colnames(unclass(out$p)), "globaltest")
  expect_equal(rownames(out), names(gmt))
  expect_true(is.na(out$p["none", 1])) ## no measured member
  expect_true(is.na(out["none", "meta.p"])) ## NA, not max() = -Inf
})

test_that("globaltest finds the shifted sets, signed, and not the null one", {
  expect_lt(out["up", "meta.p"], 1e-3)
  expect_gt(out["up", "meta.fx"], 0)
  expect_lt(out["DB:down y (z)", "meta.p"], 1e-3)
  expect_lt(out["DB:down y (z)", "meta.fx"], 0)
  expect_lt(out["mixed", "meta.p"], 1e-3) ## non-directional: both ways count
  expect_gt(out["null", "meta.p"], 0.05)
})

test_that("without padding fisher tests small sets and vetoes the hitless", {
  ## padded, all 60 features are hits both ways; unpadded with the
  ## both-directions rule, "up" and "small" (< 6 members) are never tested
  out <- fit(c("fisher", "globaltest"), fisher.pad = FALSE)
  expect_equal(names(which(is.na(out$p[, "fisher"]))), "none")
  expect_lt(out["up", "meta.p"], 0.05)
  expect_lt(out["small", "meta.p"], 0.05)
  expect_equal(out["small.null", "meta.p"], 1) ## not globaltest's alone
  expect_true(is.na(out["none", "meta.p"]))
})
