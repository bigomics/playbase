## DESeq2 cannot estimate size factors when every gene holds a zero, as in
## single-cell counts. The fit must skip DESeq2, not fail the other methods.

zero_heavy_counts <- function() {
  set.seed(7)
  counts <- matrix(
    stats::rpois(60 * 8, 20),
    nrow = 60,
    dimnames = list(paste0("G", 1:60), paste0("S", 1:8))
  )
  counts[cbind(1:60, rep(1:8, length.out = 60))] <- 0
  counts
}

b_vs_a <- function() {
  matrix(
    rep(c(-1, 1), each = 4),
    ncol = 1,
    dimnames = list(paste0("S", 1:8), "b_vs_a")
  )
}

test_that("the no-design DESeq2 fit names the contrast it could not fit", {
  expect_error(
    suppressMessages(playbase:::.ngs.fitContrastsWithDESEQ2.nodesign(
      zero_heavy_counts(), b_vs_a()
    )),
    "DESeq2 could not fit contrast b_vs_a"
  )
})

test_that("a DESeq2 failure skips that method and keeps the others", {
  counts <- zero_heavy_counts()
  samples <- data.frame(
    group = rep(c("a", "b"), each = 4),
    row.names = colnames(counts)
  )
  res <- suppressMessages(playbase::ngs.fitContrastsWithAllMethods(
    counts = counts,
    samples = samples,
    contr.matrix = b_vs_a(),
    methods = c("trend.limma", "deseq2.wald")
  ))
  expect_named(res$outputs, "trend.limma")
})
