## Gene renaming ahead of pgx.predictCellTypes(): reference atlases use
## human symbols, so count rows go through rename_by2() on the
## human_ortholog column. Like Azimuth's ConvertGeneNames, it keeps the
## first gene per human symbol, does not add counts up, and drops genes
## without a symbol.

test_that("rename_by2 maps sparse mouse counts to human symbols without summing paralogs", {
  counts <- Matrix::sparseMatrix(
    i = 1:5, j = c(1, 1, 2, 2, 2), x = c(1, 2, 3, 4, 5), dims = c(5, 2),
    dimnames = list(
      c("ENSMUSG00000025902", "ENSMUSG00000000001", "ENSMUSG00000000002",
        "ENSMUSG00000000003", "novel1"),
      c("c1", "c2")
    )
  )
  annot <- data.frame(
    feature = rownames(counts),
    symbol = c("Sox17", "Ifi1", "Ifi2", "Ifi3", NA),
    human_ortholog = c("SOX17", "IFI", "IFI", "IFI", NA),
    row.names = rownames(counts)
  )
  out <- rename_by2(counts, annot, new_id = "human_ortholog")

  expect_s4_class(out, "dgCMatrix")
  ## first paralog keeps the human symbol with its own counts
  expect_equal(as.matrix(out)["IFI", ], c(c1 = 2, c2 = 0))
  ## later paralogs are set aside under names no reference gene has
  expect_setequal(rownames(out), c("SOX17", "IFI", "IFI.1", "IFI.2"))
  ## unmapped gene dropped
  expect_false("novel1" %in% rownames(out))
})

test_that("numeric single-cell QC thresholds filter cells", {
  counts <- Matrix::Matrix(
    rbind(
      G1 = c(10, 10, 10, 10),
      G2 = c(10, 0, 0, 0),
      G3 = c(10, 0, 0, 0),
      `MT-ND1` = c(0, 10, 0, 0),
      HBA1 = c(0, 0, 10, 0)
    ),
    sparse = TRUE
  )
  colnames(counts) <- paste0("cell", 1:4)
  samples <- data.frame(row.names = colnames(counts))
  make_object <- function(settings) {
    pgx.createSeuratObject(
      counts, samples, batch = NULL, cellcyclescores = FALSE,
      preprocess = FALSE, sc_compute_settings = settings
    )
  }

  expect_equal(ncol(make_object(list())), 4L)
  expect_equal(ncol(make_object(list(regress_mt = TRUE))), 4L)
  expect_equal(ncol(make_object(list(mt_threshold = FALSE))), 4L)
  expect_setequal(colnames(make_object(list(nfeature_threshold = c(1, 4)))),
                  paste0("cell", 1:3))
  expect_setequal(colnames(make_object(list(mt_threshold = 10))),
                  c("cell1", "cell3", "cell4"))
  expect_setequal(colnames(make_object(list(hb_threshold = 10))),
                  c("cell1", "cell2", "cell4"))
})
