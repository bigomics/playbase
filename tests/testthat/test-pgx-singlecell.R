## Gene renaming ahead of pgx.predictCellTypes(): reference atlases use
## human symbols, so count rows are renamed through the annotation table
## and rows that land on the same symbol are summed.

test_that(".sum_rows_by sums duplicate ids and keeps a sparse integer matrix", {
  counts <- Matrix::sparseMatrix(
    i = c(1, 2, 3, 3), j = c(1, 1, 2, 3), x = c(1, 2, 5, 7),
    dims = c(3, 3), dimnames = list(c("a", "b", "c"), paste0("cell", 1:3))
  )
  out <- .sum_rows_by(counts, c("X", "X", "Y"))
  expect_s4_class(out, "CsparseMatrix")
  expect_equal(rownames(out), c("X", "Y"))
  expect_equal(colnames(out), colnames(counts))
  expect_equal(as.matrix(out)["X", ], c(cell1 = 3, cell2 = 0, cell3 = 0))
  expect_equal(as.matrix(out)["Y", ], c(cell1 = 0, cell2 = 5, cell3 = 7))
})

test_that(".sum_rows_by accepts a dense matrix", {
  counts <- matrix(1:6, nrow = 3, dimnames = list(c("a", "b", "c"), c("c1", "c2")))
  out <- .sum_rows_by(counts, c("X", "Y", "X"))
  expect_equal(as.matrix(out), rbind(X = c(c1 = 4, c2 = 10), Y = c(c1 = 2, c2 = 5)))
})

test_that(".counts_to_human_symbols maps Ensembl ids and mouse genes to human symbols", {
  counts <- Matrix::sparseMatrix(
    i = 1:4, j = c(1, 1, 2, 2), x = c(1, 2, 3, 4), dims = c(4, 2),
    dimnames = list(
      c("ENSMUSG00000025902", "ENSMUSG00000033845", "ENSMUSG00000025903", "novel1"),
      c("c1", "c2")
    )
  )
  annot <- data.frame(
    symbol = c("Sox17", "Mrpl15", "Lypla1", NA),
    human_ortholog = c("SOX17", "MRPL15", "MRPL15", NA),
    row.names = rownames(counts)
  )
  out <- .counts_to_human_symbols(counts, annot)
  ## unmapped rows keep their own name; duplicates are summed
  expect_equal(rownames(out), c("SOX17", "MRPL15", "novel1"))
  expect_equal(as.matrix(out)["MRPL15", ], c(c1 = 2, c2 = 3))
  expect_equal(sum(out), sum(counts))
})

test_that(".counts_to_human_symbols falls back on symbol and leaves case alone", {
  counts <- Matrix::sparseMatrix(
    i = 1:2, j = 1:2, x = c(1, 1), dims = c(2, 2),
    dimnames = list(c("ENSG00000116922", "ENSG00000187634"), c("c1", "c2"))
  )
  annot <- data.frame(symbol = c("C1orf112", "SAMD11"), row.names = rownames(counts))
  out <- .counts_to_human_symbols(counts, annot)
  expect_equal(rownames(out), c("C1orf112", "SAMD11"))
})

test_that(".counts_to_human_symbols errors without a symbol column", {
  counts <- Matrix::sparseMatrix(i = 1, j = 1, x = 1, dims = c(1, 1),
                                 dimnames = list("g1", "c1"))
  annot <- data.frame(feature = "g1", row.names = "g1")
  expect_error(.counts_to_human_symbols(counts, annot), "human symbols")
})
