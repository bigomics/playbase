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
