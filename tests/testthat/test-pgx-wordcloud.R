##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

## Builds a pgx whose word-incidence matrix collapses to a single column:
## 10 gene sets share the word "class)", 10 contribute only unique words.
mock_pgx_one_word <- function() {
  shared <- sprintf("METABOLITE_CLASS:Alpha%d (main_class)", 1:10)
  unique <- paste0("PATHWAY:", c(
    "aardvark", "bumblebee", "crocodile", "dromedary", "elephant",
    "flamingo", "giraffe", "hedgehog", "jackal", "kangaroo"
  ))
  fx <- c(rep(3, 10), rep(0.1, 10))
  meta <- data.frame(meta.fx = fx, row.names = c(shared, unique))
  list(gset.meta = list(meta = list(contrast1 = meta)))
}

test_that("pgx.calculateWordCloud survives a single surviving word", {
  skip_if_not_installed("fgsea")
  skip_if_not_installed("Rtsne")
  skip_if_not_installed("uwot")

  res <- pgx.calculateWordCloud(mock_pgx_one_word())

  expect_false(is.null(res))
  expect_equal(nrow(res$gsea[["rms.FC"]]), 1L)
  expect_equal(dim(res$tsne), c(1L, 2L))
  expect_equal(dim(res$umap), c(1L, 2L))
})

test_that("pgx.calculateWordCloud returns NULL when no gene set yields a word", {
  ## every word is either too short or contains a digit -> no terms at all
  meta <- data.frame(meta.fx = seq_len(20), row.names = sprintf("X:g%d h%d", 1:20, 1:20))
  expect_null(pgx.calculateWordCloud(list(gset.meta = list(meta = list(c1 = meta)))))
})
