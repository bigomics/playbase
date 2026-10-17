##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

## Reference data from the lib-x release.
##
## Every reference object playbase needs is read here by a short name through
## omicspgx::libx_read, which caches each object per release and records the
## read in the open provenance context.

## lib-x path of each reference object, by short name.
.pb_LIBX_PATHS <- c(
  gsetxgene = "featuresets/legacy/gsetxgene",
  gset_size = "featuresets/legacy/gset_size.parquet",
  igsets = "featuresets/legacy/igsets.parquet",
  gset_genes = "featuresets/legacy/gset_genes.rds",
  families = "featuresets/legacy/families.parquet",
  msetxmetabolite = "featuresets/legacy/msetxmetabolite",
  metabolite_id = "annotation/metabolites/id.parquet",
  phosphosite = "annotation/proteins_special/phosphosite.parquet",
  graphite_ppi = "networks/graphite_ppi.parquet",
  signature_immunemeta = "refmat/SIGNATURE_IMMUNEMETA.rds",
  rtcga_survival = "tcga_survival.rds",
  geiger_pgx = "examples/GEIGER_PGX.rds",
  lm22 = "refmat/LM22.rds",
  immprot_signature1000 = "refmat/IMMPROT_SIGNATURE1000.rds",
  dice_signature1000 = "refmat/DICE_SIGNATURE1000.rds",
  immunostates_matrix = "refmat/IMMUNOSTATES_MATRIX.rds",
  rna_tissue_matrix = "refmat/RNA_TISSUE_MATRIX.rds",
  gtex_rna_tissue_tpm = "refmat/GTEX_RNA_TISSUE_TPM.rds",
  hpa_rna_celline = "refmat/HPA_RNA_CELLINE.rds",
  ccle_rna_celline = "refmat/CCLE_RNA_CELLINE.rds",
  ccle_rna_cancertype = "refmat/CCLE_RNA_CANCERTYPE.rds",
  reactome_dir = "annotation/pathways/reactome/svg/",
  reactome_map = "annotation/pathways/reactome/diagram_map.tsv",
  l1000_activity = "cmap/L1000-activity_n10d5263.rds",
  l1000_repurposing = "cmap/L1000_repurposing_drugs.txt",
  gmt_kea = "featuresets/raw/playdata/kinase_substrates_kea.gmt",
  gmt_chea = "featuresets/raw/playdata/tf_targets_chea.gmt",
  gmt_hgnc = "featuresets/raw/playdata/hgnc_genefamilies_EDITED.gmt"
)

## Reads the lib-x object `name` (a name of .pb_LIBX_PATHS) from the release
## `dir`. Tables and .rds files come back as R objects, .tsv and .gmt files as
## paths, `reactome_dir` as the svg directory and `l1000_repurposing` as the
## drug annotation data.frame. Errors when no lib-x release is set or the
## release lacks it.
.pb_libx <- function(name, dir = omicspgx::libx_dir()) {
  if (!is.character(name) || length(name) != 1L || !name %in% names(.pb_LIBX_PATHS)) {
    stop("[.pb_libx] unknown lib-x name: ", paste(name, collapse = ", "))
  }
  path <- .pb_LIBX_PATHS[[name]]
  switch(name,
    reactome_dir = omicspgx::libx_path(path, dir = dir),
    l1000_repurposing = utils::read.csv(omicspgx::libx_read(path, dir = dir), sep = "\t", comment.char = "#"),
    omicspgx::libx_read(path, dir = dir)
  )
}

## The lib-x release directory `dir`, or the configured release when `dir` is
## NULL (see omicspgx::libx_dir). Errors when there is none or `dir` is not a
## lib-x release.
.pb_libx_release <- function(dir = NULL) {
  dir <- omicspgx::libx_dir(dir %||% omicspgx::libx_dir())
  if (is.null(dir)) {
    stop("no lib-x release: pass libx.dir, or set LIBX_DIR or options(omicspgx.libx_dir)")
  }
  dir
}

## Reads the assets of the lib-x release `dir` whose paths match `pattern`,
## named by file name.
.pb_libx_assets <- function(pattern, dir) {
  paths <- grep(pattern, omicspgx::libx_manifest(dir)$path, value = TRUE)
  stats::setNames(lapply(paths, omicspgx::libx_read, dir = dir), basename(paths))
}

## Opens a lib-x provenance context labelled `op_id`. The returned `stamp(pgx)`
## closes it into `pgx$versions`; `close()` closes it unstamped and does
## nothing once stamped, so the caller can `on.exit(prov$close())` to close it
## on error.
.pb_libx_provenance <- function(op_id) {
  omicspgx::libx_provenance_begin(op_id)
  open <- TRUE
  end <- function() {
    open <<- FALSE
    omicspgx::libx_provenance_end()
  }
  list(
    stamp = function(pgx) .pb_libx_stamp(pgx, end()),
    close = function() if (open) end()
  )
}

## Stamps `pgx` with the provenance record `prov`. A pgx stamped before keeps
## its release and manifest, gains the new assets, and becomes "mixed" when
## the releases or reference sources differ.
.pb_libx_stamp <- function(pgx, prov) {
  old <- pgx$versions
  if (!is.null(old$libx_assets)) {
    known <- function(x) unique(x[!is.na(x)])
    releases <- known(c(old$libx_release, prov$release))
    sources <- known(c(old$reference_source, prov$reference_source))
    if (length(known(old$libx_release))) {
      prov$release <- old$libx_release
      prov$manifest_sha256 <- old$libx_manifest_sha256
    }
    prov$assets <- rbind(old$libx_assets, prov$assets)
    prov$reference_source <- if (length(releases) > 1L || length(sources) > 1L) "mixed"
                             else if (length(sources)) sources else NA_character_
  }
  omicspgx::pgx_stamp_libx(pgx, prov)
}

## Blue-grey-red palette of `n` colours.
BLUERED <- function(n) {
  grDevices::colorRampPalette(c("#3A5FCD", "#E5E5E5", "#CD5555"))(n)
}

#' @title Get the genes of gene sets
#'
#' @description Looks up the member genes of gene sets in the lib-x gene set
#' collection.
#'
#' @param gs Character vector of gene set names.
#'
#' @return A named list with, per gene set, the character vector of its genes.
#'
#' @export
getGSETS <- function(gs) {
  igsets <- .pb_libx("igsets")
  gset_genes <- .pb_libx("gset_genes")
  lapply(igsets[gs], function(i) gset_genes[i])
}
