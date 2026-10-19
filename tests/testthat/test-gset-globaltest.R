#' Test gset.globaltest, the in-house port of globaltest::gt
set.seed(537)
X <- matrix(rnorm(40 * 10, sd = rep(runif(40, 0.2, 3), 10)), 40, 10,
  dimnames = list(paste0("m", 1:40), paste0("s", 1:10))
)
y <- rep(0:1, each = 5)
X[1:5, y == 1] <- X[1:5, y == 1] + 2
X[6, ] <- 1 ## constant
X[7, 2] <- NA
gmt <- list(
  up = paste0("m", 1:5), mixed = paste0("m", c(1:3, 20:30)),
  big = paste0("m", 1:30), const = "m6", one = "m7", na = c("m7", "m8")
)

test_that("p-values are identical to globaltest::gt", {
  skip_if_not_installed("globaltest")
  for (st in c(FALSE, TRUE)) {
    gt <- globaltest::gt(y, t(X),
      subsets = gmt, model = "logistic", standardize = st
    )
    expect_identical(
      playbase::gset.globaltest(X, y, gmt, standardize = st)$p.value,
      unname(globaltest::result(gt)[, "p-value"])
    )
  }
})

test_that("gset.globaltest output and edge cases", {
  gmt2 <- c(gmt, list(none = "x1", "DB:x y (z)" = "m1", up = paste0("m", 1:5)))
  out <- playbase::gset.globaltest(X, factor(y), gmt2)
  expect_equal(rownames(out), make.unique(names(gmt2)))
  expect_lt(out["up", "p.value"], 1e-3)
  expect_gt(out["up", "score"], 0)
  expect_equal(out["const", "p.value"], 1)
  expect_true(all(is.na(out["none", ])))
  expect_equal(out["up.1", ], out["up", ], ignore_attr = TRUE)
  expect_equal(out$SIZE, c(5, 14, 30, 1, 1, 2, NA, 1, 5))
  expect_equal(nrow(playbase::gset.globaltest(X, y, list())), 0)
})
