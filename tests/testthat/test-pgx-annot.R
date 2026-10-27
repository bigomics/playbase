test_that("getProbeAnnotation annotates unmapped probes on vendor gene column", {
  ## mock: official symbols and one UniProt ID map; TF resolves as alias
  db <- data.frame(
    row.names = c("SCD", "SNX1", "P04637", "TF"),
    symbol = c("SCD", "SNX1", "TP53", "F3"),
    gene_title = c("stearoyl-CoA desaturase", "sorting nexin 1", "tumor protein p53", "tissue factor")
  )
  calls <- list()
  testthat::local_mocked_bindings(
    getGeneAnnotation = function(organism, probes, ...) {
      calls[[length(calls) + 1]] <<- probes
      g <- data.frame(feature = probes, db[probes, ], row.names = make.unique(probes))
      g$human_ortholog <- g$symbol
      playbase:::cleanupAnnotation(g)
    }
  )
  probes <- c("Stearoyl-CoA desaturase_SCD", "Sorting nexin-1_SNX1", "P04637", "Serotransferrin_TF")
  annot <- data.frame(
    row.names = probes,
    PG.Genes = c("SCD", "SNX1;SNX1B", "TP53", "TF"),
    PG.UniProtIds = c("O00767", "Q13596", "P04637", "P02787")
  )
  g <- getProbeAnnotation("Human", probes, datatype = "proteomics", annot_table = annot)

  expect_equal(g$feature, probes)
  ## TF -> F3 is not the user's gene: row left as it was
  expect_equal(g$symbol, c("SCD", "SNX1", "TP53", "{Serotransferrin_TF}"))
  expect_equal(g$human_ortholog, c("SCD", "SNX1", "TP53", ""))
  expect_equal(g$gene_title[1:3], c("stearoyl-CoA desaturase", "sorting nexin 1", "tumor protein p53"))
  expect_equal(g$PG.Genes, annot$PG.Genes)
  ## only the rows that did not map are re-annotated
  expect_equal(calls[[2]], c("SCD", "SNX1", "TF"))
})
