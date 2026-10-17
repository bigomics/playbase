make_gmt <- function() {
  gmt <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 2, 3, 4, 5, 6),
    j = c(1, 1, 1, 2, 2, 3, 3, 3),
    x = 1, dims = c(6, 3),
    dimnames = list(paste0("G", 1:6), c("A:set1", "B:set2", "C:set3"))
  )
  gmt
}

test_that("sizes come from the catalog in column order", {
  gmt <- make_gmt()
  catalog <- data.frame(
    name = colnames(gmt), set_uid = c("u1", "u2", "u3"),
    n_members = c(10L, 20L, 30L)
  )
  pgx <- list(GMT = gmt, gset.catalog = catalog)
  gset.size <- Matrix::colSums(pgx$GMT != 0)
  out <- gset_size_raw(pgx, list(), gset.size)
  expect_equal(out, c("A:set1" = 10L, "B:set2" = 20L, "C:set3" = 30L))
})

test_that("catalog rows in another order are aligned by name", {
  gmt <- make_gmt()
  catalog <- data.frame(
    name = c("C:set3", "A:set1", "B:set2"), set_uid = c("u3", "u1", "u2"),
    n_members = c(30L, 10L, 20L)
  )
  pgx <- list(GMT = gmt, gset.catalog = catalog)
  gset.size <- Matrix::colSums(pgx$GMT != 0)
  out <- gset_size_raw(pgx, list(), gset.size)
  expect_equal(out, c("A:set1" = 10L, "B:set2" = 20L, "C:set3" = 30L))
})

test_that("without a catalog the legacy size vector is used", {
  gmt <- make_gmt()
  pgx <- list(GMT = gmt)
  gset.size <- Matrix::colSums(pgx$GMT != 0)
  testthat::local_mocked_bindings(
    .pb_libx = function(name) c("A:set1" = 5, "B:set2" = 6, "X:other" = 7)
  )
  custom <- list(gmt = list(), info = list(GSET_SIZE = c("C:set3" = 9)))
  out <- gset_size_raw(pgx, custom, gset.size)
  expect_equal(out, c("A:set1" = 5, "B:set2" = 6, "C:set3" = 9))
})

test_that("same-name sets take the first set's size, as their GMT columns do", {
  gmt <- make_gmt()
  colnames(gmt)[2] <- "A:set1"
  catalog <- data.frame(
    name = colnames(gmt), set_uid = c("u1", "u2", "u3"),
    n_members = c(10L, 20L, 30L)
  )
  ## compute_testGenesets subsets the GMT by set name
  pgx <- list(GMT = gmt[, colnames(gmt)], gset.catalog = catalog)
  gset.size <- Matrix::colSums(pgx$GMT != 0)
  expect_equal(unname(gset.size), c(3, 3, 3))
  out <- gset_size_raw(pgx, list(), gset.size)
  expect_equal(unname(out), c(10L, 10L, 30L))
})
