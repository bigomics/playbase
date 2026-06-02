test_that("pathway report tables render integrated cross-contrast evidence", {
  skip_if_not_installed("omicsai")

  terms <- c(
    "GOBP:rRNA processing",
    "PATHWAY_REACTOME:DNA Replication_R-HSA-69306",
    "CELLTYPE:irrelevant"
  )
  make_meta <- function(fx, q) {
    data.frame(
      meta.fx = fx,
      meta.p = q / 2,
      meta.q = q,
      q = I(cbind(method_a = q, method_b = q)),
      row.names = terms
    )
  }

  pgx <- list(
    name = "synthetic pathways",
    description = "high light response in two genotypes",
    organism = "arabidopsis",
    datatype = "RNA-seq",
    samples = data.frame(
      genotype = c("wild type", "wild type", "knockdown", "knockdown"),
      treatment = c("high light", "none", "high light", "none"),
      row.names = paste0("s", 1:4)
    ),
    contrasts = cbind(
      highlight_vs_none = c("highlight", "none", NA, NA),
      highlight_vs_none.1 = c(NA, NA, "highlight", "none.1")
    )
  )
  slice <- list(meta = list(
    highlight_vs_none = make_meta(c(1.2, -0.8, 5), c(0.001, 0.01, 0.001)),
    highlight_vs_none.1 = make_meta(c(0.9, -1.1, 6), c(0.02, 0.002, 0.001))
  ))

  out <- pathways_build_report_tables(slice, pgx, ntop = 2L, cross_n = 3L)

  expect_match(out$text, "synthetic pathways")
  expect_match(out$text, "Cross-contrast pathway footprint")
  expect_match(out$text, "genotype=wild type\\(1\\)")
  expect_match(out$text, "genotype=knockdown\\(1\\)")
  expect_match(out$text, "rRNA processing")
  expect_equal(.pathways_clean_term("PATHWAY:Ribosome_Homo sapiens_hsa03010"),
               "Ribosome [hsa03010]")
  expect_false(grepl("CELLTYPE", out$text))
  expect_equal(nrow(out$data$contrast_landscape), 2L)
  expect_equal(nrow(out$data$cross_contrast), 2L)
  expect_false(grepl("\\{\\{[^}]+\\}\\}", out$text, perl = TRUE))
})
