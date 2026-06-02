test_that("DE report tables render integrated cross-contrast evidence", {
  skip_if_not_installed("omicsai")

  features <- c("g1", "g2", "g3")
  make_meta <- function(fx, q) {
    data.frame(
      meta.fx = fx,
      meta.p = q / 2,
      meta.q = q,
      avg.0 = c(4, 5, 6),
      avg.1 = c(6, 4, 7),
      fc = I(cbind(method_a = fx, method_b = fx * 0.9)),
      p = I(cbind(method_a = q / 2, method_b = q / 3)),
      q = I(cbind(method_a = q, method_b = q)),
      row.names = features
    )
  }

  pgx <- list(
    name = "synthetic",
    description = "two related comparisons",
    organism = "human",
    datatype = "RNA-seq",
    samples = data.frame(
      group = c("A", "A", "B", "B"),
      row.names = paste0("s", 1:4)
    ),
    contrasts = cbind(
      A_vs_B = c("A", "A", "B", "B"),
      A2_vs_B = c("A2", "A2", "B", "B")
    ),
    genes = data.frame(
      symbol = c("GENE1", "GENE2", "GENE3"),
      gene_title = c("first gene", "second gene", "third gene"),
      row.names = features
    ),
    X = matrix(0, nrow = 3, ncol = 4)
  )
  slice <- list(meta = list(
    A_vs_B = make_meta(c(2, -1, 0.5), c(0.001, 0.01, 0.8)),
    A2_vs_B = make_meta(c(1.5, -0.7, 0.2), c(0.002, 0.04, 0.9))
  ))

  out <- de_build_report_tables(slice, pgx, ntop = 2L, cross_n = 3L)

  expect_match(out$text, "synthetic")
  expect_match(out$text, "Cross-contrast marker footprint")
  expect_match(out$text, "GENE1")
  expect_equal(nrow(out$data$contrast_summary), 2L)
  expect_equal(nrow(out$data$cross_contrast), 3L)
  expect_false(grepl("\\{\\{[^}]+\\}\\}", out$text, perl = TRUE))
})
