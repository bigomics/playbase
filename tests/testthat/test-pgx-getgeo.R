#' Test for title2pheno
test_that("title2pheno extracts phenotype terms", {
  titles <- c(
    "GSE1234-tissue_liver-disease-cancer",
    "GSE4321-tissue_lung-disease-cancer"
  )

  expected <- list(
    tissue = c("liver", "lung", "muscle"),
    disease = c("cancer", "IPF", "myopathy")
  )

  result <- title2pheno(titles, split = "-")

  # Check class, dim, and first column
  expect_equal(class(result), c("matrix", "array"))
  expect_equal(dim(result), c(length(titles), 4))
  expect_equal(result[, 1], c("GSE1234", "GSE4321"))
})
