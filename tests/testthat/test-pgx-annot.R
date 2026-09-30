## getProbeAnnotation() for methylation arrays: every probe carries its own
## manifest chromosome, a probe outside any gene included, while the gene
## columns stay as playbase.epigenetics annotated them.
for (arr in list(
  list(meth_type = "450K array", pkg = "IlluminaHumanMethylation450kanno.ilmn12.hg19"),
  list(meth_type = "EPIC array", pkg = "IlluminaHumanMethylationEPICanno.ilm10b4.hg19"),
  list(meth_type = "EPIC v2 array", pkg = "IlluminaHumanMethylationEPICv2anno.20a1.hg38")
)) {
  test_that(paste(arr$meth_type, "probes carry their manifest chromosome, gene or not"), {
    skip_if_not_installed("playbase.epigenetics")
    skip_if_not_installed(arr$pkg)
    manifest <- get("methyl_annotation", envir = asNamespace("playbase.epigenetics"))(arr$meth_type)
    probes <- rownames(manifest)
    genic <- nzchar(manifest$UCSC_RefGene_Name)
    ## annotate_methylomics() as it is: the gene's cytoband in chr, nothing
    ## for a probe without a gene, the manifest's pos for every probe.
    local_mocked_bindings(
      annotate_methylomics = function(organism, probes, meth_type) {
        g <- data.frame(
          feature = probes, symbol = ifelse(genic, manifest$UCSC_RefGene_Name, NA),
          chr = ifelse(genic, "1p36.33", NA), pos = manifest$pos, row.names = probes
        )
        attr(g, "genome") <- "genome-attr"
        g
      },
      .package = "playbase.epigenetics"
    )

    genes <- getProbeAnnotation("Human", probes, "methylomics", meth_type = arr$meth_type)

    expect_identical(genes$chr, as.character(manifest$chr))
    expect_false(anyNA(genes$chr) || any(!nzchar(genes$chr)))
    expect_true(any(!genic))
    expect_true(all(is.na(genes$symbol[!genic])))
    expect_identical(genes$pos, manifest$pos)
    expect_identical(attr(genes, "genome"), "genome-attr")
  })
}

## An ALIAS key that is also another gene's official symbol must resolve to
## that gene, not to whichever gene sorts first: "TPO" is thyroid peroxidase
## (and an alias of THPO), "EBF3" is EBF3 (and an alias of MAPRE3).
test_that("an official symbol wins over another gene's alias", {
  annot <- data.frame(
    ALIAS = c("TPO", "TPO", "EBF3", "EBF3", "GAPDH", "OLDNAME", "AMBIG", "AMBIG"),
    SYMBOL = c("THPO", "TPO", "MAPRE3", "EBF3", "GAPDH", "NEWNAME", "GENEA", "GENEB"),
    GENENAME = c("thrombopoietin", "thyroid peroxidase", "MAPRE3 title", "EBF3 title",
                 "gapdh", "renamed gene", "a", "b")
  )
  out <- resolve_alias_rows(annot, "ALIAS")
  expect_identical(out$SYMBOL[out$ALIAS == "TPO"], "TPO")
  expect_identical(out$SYMBOL[out$ALIAS == "EBF3"], "EBF3")
  expect_identical(out$SYMBOL[out$ALIAS == "GAPDH"], "GAPDH")
  ## An alias of exactly one gene still resolves; an ambiguous one does not.
  expect_identical(out$SYMBOL[out$ALIAS == "OLDNAME"], "NEWNAME")
  expect_false("AMBIG" %in% out$ALIAS)
  ## Other keytypes pass through untouched.
  expect_identical(resolve_alias_rows(annot, "SYMBOL"), annot)
})

test_that("gene annotation by ALIAS keeps the official symbol's gene", {
  skip_if_not_installed("org.Hs.eg.db")
  orgdb <- org.Hs.eg.db::org.Hs.eg.db
  keys <- c("TPO", "EBF3", "TNRC18", "AHRR", "GAPDH")
  a <- suppressMessages(AnnotationDbi_select_2pass(
    orgdb, keys = keys, columns = c("SYMBOL", "GENENAME"), keytype = "ALIAS"
  ))
  expect_identical(a$SYMBOL, keys)
  expect_identical(a$GENENAME[1:2], c("thyroid peroxidase", "EBF transcription factor 3"))
})

## The loose match (case and punctuation folded) is a fallback only: an exact
## key must win, or "SF3B1" becomes SF3B2's alias "SF3b1" and "NKX6-1" the
## alias "NKX6.1".
test_that("probe names keep an exact key over a case- or punctuation-folded one", {
  skip_if_not_installed("org.Hs.eg.db")
  orgdb <- org.Hs.eg.db::org.Hs.eg.db
  keys <- c("SF3B1", "SLIT2", "TIFA", "NKX6-1", "VMA21")
  expect_identical(unname(match_probe_names(keys, orgdb, "ALIAS")), keys)
  ## Folding still rescues a key that has no exact match.
  expect_identical(unname(match_probe_names("gapdh", orgdb, "SYMBOL")), "GAPDH")
})
