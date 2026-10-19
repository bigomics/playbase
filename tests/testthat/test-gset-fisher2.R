#' Test gset.fisher2 direction handling
test_that("gset.fisher2 matches set names exactly, per direction", {
  ## each set has hits in one list only; "s21" is a prefix of "s214"
  gmt <- list(s21 = paste0("g", 11:15), s214 = paste0("g", 1:8))
  ft <- function(...) {
    playbase::gset.fisher2(paste0("g", 1:10), paste0("g", 11:20), gmt,
      background = paste0("g", 1:100), fdr = 1, min.genes = 0, verbose = 0,
      ...
    )
  }
  expect_equal(nrow(ft()), 0) ## neither has hits in both directions
  out <- ft(both.directions = FALSE)
  expect_equal(sign(out[c("s21", "s214"), "sign"]), c(-1, 1))
  expect_equal(out[c("s21", "s214"), "overlap"], c("5/5", "8/8"))
})
