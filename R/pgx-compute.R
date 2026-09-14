##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

#' Create a pgx object
#' This function creates a pgx object from files, which is the core object in the
#' OmicsPlayground. It then runs the specified differential expression methods.
#' @param counts.file Path to counts data file. Rows are genes, columns are samples.
#' @param samples.file Path to samples data file. Rows are samples, columns are sample info.
#' @param contrasts.file (optional) Path to contrasts file. Rows and columns define contrasts.
#' @param preprocess (optional) Named list of preprocessing settings, forwarded to
#'   [pgx.createPGX()] and from there to [pgx.preprocess()]. NULL (the default)
#'   keeps the historical behaviour of this entry point, where `X` is a plain
#'   `log2(counts + prior)` with no filtering, imputation or normalization.
#'   Supplying a list is what makes a script reproduce the app: the wizard sends
#'   the same list, so the same settings give the same `X`. See
#'   \code{\link[playbase.preprocess]{pgx.preprocess}} for the settings and
#'   their defaults.
#' @param gxmethods a string with the gene-level methods to use. The default value is \code{"trend.limma,edger.qlf,deseq2.wald"}
#' @param gsetmethods a string with the gene-set methods to use. The default value is \code{"fisher,gsva,fgsea"}
#' @param extra a string with the extra modules to use. The default value is \code{"meta.go,deconv,infer,drugs,wordcloud"}
#' @return list. represents a pgx object. It contains the data and analysis results.
#' @examples
#' \dontrun{
#' library(playbase)
#' counts <- system.file("extdata", "counts.csv", package = "playbase")
#' contrasts <- system.file("extdata", "contrasts.csv", package = "playbase")
#' samples <- system.file("extdata", "samples.csv", package = "playbase")
#' mypgx <- pgx.createFromFiles(counts, samples, contrasts)
#' }
#' @export
pgx.createFromFiles <- function(counts.file,
                                samples.file,
                                contrasts.file = NULL,
                                preprocess = NULL,
                                gxmethods = "trend.limma,edger.qlf,deseq2.wald",
                                gsetmethods = "fisher,gsva,fgsea",
                                extra = "meta.go,deconv,infer,drugs,wordcloud",
                                pgx.dir = "./data",
                                libx.dir = "./libx") {
  ## read counts table (allow dup rownames)
  counts <- read.as_matrix(counts.file)

  ## compile sample table
  samples <- read.as_matrix(samples.file)
  samples <- data.frame(samples, check.names = FALSE)

  ## parse requested phenotypes
  if (!is.null(contrasts.file) && file.exists(contrasts.file)) {
    message("reading contrasts file ", contrasts.file)
    contrasts <- read.as_matrix(contrasts.file)
  } else {
    ## take first (not-dotted) column in samples as phenotype vector
    group.col <- head(grep("group|condition", colnames(samples), ignore.case = TRUE), 1)
    if (length(group.col) == 0) {
      group.col <- head(grep("^.*", colnames(samples), invert = TRUE), 1)
    }
    if (length(group.col) == 0) {
      group.col <- colnames(samples)[1]
    }
    Y <- samples[, group.col, drop = FALSE]
    ## automatically guess contrasts
    contr <- pgx.makeAutoContrasts(Y, mingrp = 3, slen = 20, ref = NA)
    contrasts <- contrastAsLabels(contr$exp.matrix)
  }

  ## other params
  gx.methods <- strsplit(gxmethods, split = ",")[[1]]
  gset.methods <- strsplit(gsetmethods, split = ",")[[1]]
  extra.methods <- strsplit(extra, split = ",")[[1]]

  ## create initial PGX object
  pgx <- pgx.createPGX(
    counts,
    samples = samples,
    contrasts = contrasts,
    X = NULL,
    preprocess = preprocess,
    is.logx = NULL,
    dotimeseries = FALSE,
    batch.correct.method = "no_batch_correct",
    batch.pars = "<autodetect>",
    covariates = NULL,
    auto.scale = TRUE,
    filter.genes = TRUE,
    prune.samples = FALSE,
    only.known = TRUE,
    average.duplicated = FALSE,
    only.hugo = TRUE,
    convert.hugo = TRUE,
    only.proteincoding = TRUE,
    max.genesets = 10000
  )

  ## start computing PGX object
  pgx <- pgx.computePGX(
    pgx,
    max.genes = 40000,
    gx.methods = gx.methods,
    gset.methods = gset.methods,
    extra.methods = extra.methods,
    cluster.contrasts = FALSE,
    do.clustergenes = TRUE,
    do.clustergenesets = TRUE,
    do.cluster = TRUE,
    use.design = FALSE,
    prune.samples = TRUE,
    pgx.dir = pgx.dir,
    libx.dir = libx.dir,
    progress = NULL
  )

  ## save
  pgx
}

## The rows and samples of `counts` that `X` is made of, as `list(rows, cols)`.
##
## `pgx.preprocess()` returns `counts` at the shape it was handed and `X` at
## whatever the removals left it (D-24/D-39), so the two are no longer one row
## set nor one sample set and nothing downstream may assume they are. The index
## is derived, not looked up: the provenance record stays inside
## playbase.preprocess (D-37), and `pgx.alignXtoCounts()` recovers the relation
## from the data and the fixed pipeline order instead -- refusing, never
## guessing, when the names cannot answer it. While the two shapes still agree
## it is the identity, which is what the caller-supplied `X` and plain-log2
## paths get. `as.matrix()` is the coercion, not a copy: the single-cell path
## hands `pgx.createPGX()` a sparse `X`.
counts_index_of_X <- function(counts, X) {
  pgx.alignXtoCounts(list(counts = as.matrix(counts), X = as.matrix(X)))
}

## The rows only, and the identity when there is no `X` yet to align to.
counts_rows_of_X <- function(counts, X) {
  if (is.null(X)) {
    return(seq_len(nrow(counts)))
  }
  counts_index_of_X(counts, X)$rows
}

## The rows of `pgx$X` that stand for a set of `counts` rows, in the order that
## set names them. The inverse of the above, and what a feature filter needs: it
## decides on `counts`, and the index it produces cannot be applied to `X`. An
## `X` row whose counts row is not in the set is dropped, which is the filter
## doing its job.
x_rows_for_counts_rows <- function(pgx, keep) {
  i <- match(keep, counts_rows_of_X(pgx$counts, pgx$X))
  i[!is.na(i)]
}

#' @title This object's data on the count scale
#'
#' @description
#' The one back-transform rule (D-07). Everything that needs count-scale values
#' -- the edgeR/DESeq2 fitters, the deconvolution mixture, the cell-cycle and
#' gender signatures -- asks here rather than reaching for `pgx$counts` or
#' calling [pgx.recomputeCounts()] on its own.
#'
#' Batch correction moves `pgx$X` away from `pgx$counts`, and only then is the
#' upload the wrong answer: the reconstruction is what `X` is now made of. When
#' correction did not run, `X` is a transform of `counts` and `counts` is the
#' better count-scale matrix of the two, because the reconstruction is lossy
#' where normalization and imputation touched the data.
#'
#' [pgx.ranWithCorrection()] answers `NA` on an object that carries no
#' preprocessing record, and "we do not know" is not "it ran" -- so the upload
#' is used. playbase's own batch correction runs outside the record (D-37) and
#' writes nothing to it, so on every object playbase creates today the answer is
#' the upload. That is a gap in the record, not in the rule.
#'
#' The gap is not silent. `pgx.createPGX()` stamps `settings$batch.correct.method`
#' where it corrects, so this function can tell the two record-less cases apart:
#' an object that was never corrected, where the upload is simply the right
#' answer, and an object whose `X` was corrected by a step the record cannot
#' invert, where the upload is the only answer available but no longer the same
#' data `X` holds. The second **warns**, because it splits a result by method
#' rather than by contrast: `edgeR`/`DESeq2` and the deconvolution mixture come
#' here and get uncorrected values, while `trend.limma` and `ttest` read the
#' corrected `pgx$X` directly. Closing it needs the correction to reach the
#' record; until then the split is visible instead of silent.
#'
#' The two branches agree on scale, **not** on shape: the reconstruction has
#' `X`'s rows and samples, the upload has its own (D-24). A caller that needs
#' the shapes to line up subsets by name, or asks [pgx.alignXtoCounts()].
#'
#' @param pgx A pgx-shaped list carrying `counts` and `X`.
#'
#' @return A matrix on the count scale. Warns, and returns `pgx$counts`, when
#'   the object records a batch correction the preprocessing record cannot
#'   invert.
#'
#' @export
pgx.countScaleMatrix <- function(pgx) {
  if (isTRUE(pgx.ranWithCorrection(pgx))) {
    return(pgx.recomputeCounts(pgx))
  }
  mm <- pgx$settings$batch.correct.method
  if (!is.null(mm)) {
    warning(
      "[pgx.countScaleMatrix] this object was batch-corrected with '", mm,
      "', and that correction was applied to X outside the preprocessing ",
      "record, so there is no chain to invert and the count scale cannot ",
      "follow it. Returning the uploaded counts. Count-scale methods (edgeR, ",
      "DESeq2, deconvolution) therefore answer on UNCORRECTED data while ",
      "trend.limma and ttest read the corrected X: the same contrast can ",
      "differ by method.",
      call. = FALSE
    )
  }
  pgx$counts
}

#' Create a PGX object
#' This function creates a pgx object, which is the core object in the
#' OmicsPlayground.
#' @param counts Matrix of count data with genes as rows and samples as columns.
#' @param samples Data frame containing sample information.
#' @param organism Default "Human", it indicates the species used
#' for the gene annotation table and the probe to symbol conversion.
#' @param contrasts Data frame defining sample contrasts.
#' @param X (Optional) Matrix of normalized expression data. If NULL, will be calculated from counts.
#' @param preprocess (Optional) Named list of preprocessing options. If provided and
#'   `X` is NULL, `X` is built from `counts` via [pgx.preprocess()] (normalization,
#'   imputation, missingness filter, outlier removal) instead of a plain log2 transform.
#'   This is how the Shiny upload flow and the compute endpoint obtain identical `X`.
#' @param is.logx Logical indicating if count matrix is already log-transformed. If NULL, guessed automatically.
#' @param dotimeseries Logical indicating if timeseries analysis has been activated by the user at upload
#' @param batch.correct.method BC method. Default is "no_batch_correct" (meaning no batch correction).
#' @param batch.pars BC variable. Default "autodetect" as per QC/BC tab in upload.
#' @param covariates variables to regress out. Valid only for linear model-based tests.
#' @param dma Differential methylation analysis. If datatype=="methylomics", can be DMP (default) vs. DMR. Else NULL.
#' @param remove.xy.probes Logical. Only activated when datatype=="methylomics". Remove X- and Y-linked CpG probes.
#' @param meth_type Type of array: 450K array or EPIC array
#' @param auto.scale Logical indicating whether to automatically scale/center genes. Default is TRUE.
#' @param filter.genes Logical indicating whether to filter lowly expressed genes. Default is TRUE.
#' @param prune.samples Logical indicating whether to remove samples without contrasts. Default is FALSE.
#' @param only.known Logical indicating whether to keep only known genes. Default is TRUE.
#' @param average.duplicated Logical whether average duplicated features (if any). Default FALSE (thus keep all by making unique).
#' @param only.hugo Logical indicating whether to convert symbols to HUGO names. Default is TRUE.
#' @param convert.hugo Logical indicating whether to convert symbols to HUGO names. Default is TRUE.
#' @param only.proteincoding Logical indicating whether to keep only protein-coding genes. Default is TRUE.
#' @param custom.geneset Custom gene sets to test, as a named list with gmt and info elements.
#' @param max.genesets Maximum number of gene sets to test. Default is 5000.
#'
#' @details
#' pgx.createPGX creates a pgx object with the following slots:
#'
#' - `name`: Name of the dataset
#' - `organism`: Organism for the dataset
#' - `version`: Dataset version
#' - `date`: Date the dataset was created
#' - `creator`: Creator of the dataset
#' - `datatype`: Type of data (e.g. RNA-seq, microarray)
#' - `description`: Description of the dataset
#' - `metadata`: User-defined metadata (list with study_type, tissue_type, etc.)
#' - `samples`: Sample metadata
#' - `counts`: Raw count matrix
#' - `contrasts`: Contrast matrix
#' - `X`: Normalized expression matrix
#' - `total_counts`: Total counts per sample
#' - `counts_multiplier`: Counts multiplier for each sample
#' - `genes`: Gene annotation data.frame (initially NULL)
#' - `all_genes`: Full list of genes
#' - `probe_type`: Probe type according to biomaRt classification(e.g. ensemble_id)
#' - `GMT`: Gene set matrix
#' @import data.table
#' @return List. PGX object containing input data and parameters.
#' @export
pgx.createPGX <- function(counts,
                          samples,
                          contrasts,
                          organism,
                          custom.geneset = NULL,
                          annot_table = NULL,
                          max.genesets = 5000,
                          name = "Data set",
                          datatype = "RNA-seq",
                          datatype_subtype = NULL,
                          azimuth_ref = "pbmcref",
                          probe_type = NULL,
                          creator = "unknown",
                          description = "No description provided.",
                          metadata = NULL,
                          X = NULL,
                          preprocess = NULL,
                          norm_method = "CPM",
                          is.logx = NULL,
                          dotimeseries = FALSE,
                          batch.correct.method = "no_batch_correct",
                          batch.pars = "<autodetect>",
                          covariates = NULL,
                          dma = NULL, ## new
                          remove.xy.probes = FALSE, ## new
                          meth_type = NULL, ## new
                          auto.scale = TRUE,
                          filter.genes = TRUE,
                          exclude.genes = NULL,
                          prune.samples = FALSE,
                          only.known = TRUE,
                          average.duplicated = FALSE,
                          only.hugo = TRUE, ## DEPRECATED
                          convert.hugo = FALSE,
                          only.proteincoding = TRUE,
                          remove.xxl = TRUE, ## DEPRECATED
                          remove.outliers = TRUE, ## DEPRECATED
                          add.gmt = TRUE,
                          ortholog_species = "Human",
                          include_default_gmt = TRUE,
                          species_go = NULL,
                          settings = list(),
                          sc_compute_settings = list()) {
  message("[pgx.createPGX]===========================================")
  message("[pgx.createPGX]=========== pgx.createPGX =================")
  message("[pgx.createPGX]===========================================")
  message("\n")
  message("[pgx.createPGX] datatype = ", datatype, "\n")
  if (!is.null(datatype_subtype)) {
    message("[pgx.createPGX] datatype_subtype = ", datatype_subtype, "\n")
  }

  if (is.null(counts)) stop("[pgx.createPGX] FATAL: counts must be provided")
  if (is.null(samples)) stop("[pgx.createPGX] FATAL: samples must be provided")
  if (is.null(organism)) stop("[pgx.createPGX] FATAL: organism must be provided")

  message("[pgx.createPGX] dim.counts: ", dim(counts)[1], " x ", dim(counts)[2])
  message("[pgx.createPGX] class.counts: ", class(counts))
  message("[pgx.createPGX] counts has ", sum(is.na(counts)), " missing values")

  ## Opt-in: build X from raw counts via the shared preprocessing pipeline.
  ## Runs before de-duplication so counts and X are averaged/uniquified together,
  ## matching the Shiny upload flow (normalization module -> createPGX).
  if (is.null(X) && !is.null(preprocess) && datatype != "scRNA-seq") {
    message("[pgx.createPGX] building X via pgx.preprocess()")
    ## `contrasts` is never validated above and is legitimately NULL for an
    ## upload with no comparisons defined yet; pgx.preprocess() tolerates that
    ## and treats the samples as one single group.
    pp <- pgx.preprocess(counts,
      samples = samples, contrasts = contrasts,
      annot = annot_table, options = preprocess
    )
    ## Phase 0 Task 0.2: persist the caller's option list verbatim (D-18).
    ## Additive only - nothing reads this key yet. NB this is the caller's
    ## PARTIAL list; the resolved opt is built inside pgx.preprocess() and
    ## is not returned. See bead for that limitation.
    settings$options <- preprocess
    ## `counts` comes back at the shape it went in (D-39). The removals shrank
    ## `X` and recorded the index they kept instead of applying it to the
    ## upload, so `X` is the smaller thing derived from it and every step from
    ## here that renames or subsets features asks counts_rows_of_X() which of
    ## the two it is talking about.
    ##
    ## That holds THROUGH PREPROCESSING, which is as far as D-24/D-39 reach. It
    ## is not a claim about what finally lands in `pgx$counts`: the three
    ## ANNOTATION filters further down -- `filter.genes`, the
    ## `only.known`/`only.proteincoding`/`exclude.genes` block, and methylomics
    ## `remove.xy.probes` -- still cut it, because `pgx$genes` is counts-shaped
    ## and the gene-set universe is built from `pgx$genes$symbol`. See the note
    ## at the first of them.
    counts <- pp$counts
    X <- pp$X
    if (!is.null(annot_table)) annot_table <- pp$annot
  }

  ndup <- sum(duplicated(rownames(counts)))
  if (ndup > 0) {
    if (average.duplicated) {
      message("[pgx.createPGX] ", ndup, " duplicated feature(s) detected. Averaging....")
      ## Asked BEFORE the merge, like the make-unique branch below, and for a
      ## sharper reason: the merge is what would hide the answer. `counts` and
      ## `X` are averaged separately over row sets a preprocessing removal may
      ## have made different (D-24), and the merge leaves BOTH name vectors
      ## unique -- so afterwards pgx.alignXtoCounts() sees no repeated name to
      ## object to and matches `counts["A"]`, the mean of every copy, to
      ## `X["A"]`, the mean of only the copies `X` kept. It answers, wrongly,
      ## where its own contract is to refuse. Here the copies are still there
      ## and the split-duplicate refusal fires on its own.
      if (!is.null(X)) invisible(counts_rows_of_X(counts, X))
      counts <- playbase::counts.mergeDuplicateFeatures(counts, is.counts = TRUE)
      if (!is.null(X)) X <- playbase::counts.mergeDuplicateFeatures(X, is.counts = FALSE)
    } else {
      message("[pgx.createPGX] ", ndup, " duplicated feature(s) detected. Making unique to keep all...")
      ## One naming, taken from `counts` and handed down to the rows `X` kept.
      ## make_unique() numbers the copies it can see, so renaming `X` from its
      ## own rownames would give one feature a different suffix on each side
      ## wherever a removal split a duplicate.
      xrows <- counts_rows_of_X(counts, X)
      rownames(counts) <- playbase::make_unique(rownames(counts))
      if (!is.null(X)) rownames(X) <- rownames(counts)[xrows]
      if (!is.null(annot_table)) rownames(annot_table) <- rownames(counts)
    }
  }

  if (datatype == "scRNA-seq") {
    pgx <- pgx.createSingleCellPGX(
      counts = counts,
      samples = samples,
      contrasts = contrasts,
      organism = organism,
      azimuth_ref = azimuth_ref,
      batch = NULL,
      sc_compute_settings = sc_compute_settings
    )
    return(pgx)
  }

  if (is.null(X)) {
    min.nz <- min(counts[counts > 0], na.rm = TRUE)
    prior <- ifelse(grepl("CPM|TMM|TPM", norm_method), 1, min.nz)
    message("[pgx.createPGX] creating X as log2(counts+p) with p = ", prior)
    X <- log2(counts + prior)
  }

  if (!is.null(annot_table)) {
    message("[pgx.createPGX] dim(annot_table) = ", nrow(annot_table), " x ", ncol(annot_table))
    ndiff <- sum(rownames(annot_table) != rownames(counts))
    message("[pgx.createPGX] WARNING: annot_table has ", ndiff, " different rownames as counts")
    ndups <- sum(duplicated(rownames(annot_table)))
    message("[pgx.createPGX] annot_table has ", ndups, " duplicated rows")
    if (nrow(annot_table) != nrow(counts)) {
      message("[pgx.createPGX] WARNING: annot_table has different nrows. forcing dimensions.")
      ii <- match(rownames(counts), rownames(annot_table))
      annot_table <- annot_table[ii, ]
      rownames(annot_table) <- rownames(counts)
    }
  }

  if (sum(is.na(X)) > 0) {
    message("[pgx.createPGX] X has ", sum(is.na(X)), " missing values")
  }

  ## D-24 deletes the two guards that stood here -- `dim(counts) == dim(X)` and
  ## `rownames(counts) == rownames(X)`. They are what forced the trim: they made
  ## a shrunk `X` illegal unless `counts` was cut down with it, and a cut-down
  ## `counts` is what ratchets, because upload_server.R:1095 seeds the next
  ## Reanalyse from it. What stands in their place is the weaker statement that
  ## is still true -- every row and sample of `X` is one `counts` still has --
  ## and counts_rows_of_X() stops when it is not. The index it returns is
  ## discarded on purpose: this call is the assertion, and `invisible()` says so,
  ## because a bare expression whose only effect is its `stop()` reads as dead
  ## code to the next cleanup.
  if (!is.null(X)) invisible(counts_rows_of_X(counts, X))

  if (datatype == "multi-omics") {
    has.colons <- mean(grepl("[:]", rownames(counts)), na.rm = TRUE) > 0.9
    if (!has.colons) stop("[pgx.createPGX] FATAL: features must have multi-omics prefix\n")
  }

  ## -------------------------------------------------------------------
  ## clean up input files
  ## -------------------------------------------------------------------
  samples <- as.data.frame(samples, drop = FALSE)
  counts <- as.matrix(counts)
  X <- as.matrix(X)
  if (is.null(contrasts)) contrasts <- samples[, 0]

  ## convert old-style contrast matrix to sample-wise labeled contrasts
  contrasts <- contrasts.convertToLabelMatrix(contrasts, samples)
  contrasts <- fixContrastMatrix(contrasts)

  ## ---------------------------------------------------------------------
  ## Time series conducted if user checked the box during upload
  ## ---------------------------------------------------------------------
  if (dotimeseries) contrasts <- contrasts.addTimeInteraction(contrasts, samples)

  ## -------------------------------------------------------------------
  ## Auto-scaling (scale down huge values, often in proteomics)
  ## -------------------------------------------------------------------
  # res <- counts.autoScaling(counts)
  # counts <- res$counts
  # counts_multiplier <- res$counts_multiplier
  counts_multiplier <- Inf
  # remove(res)

  ## -------------------------------------------------------------------
  ## conform all matrices
  ## -------------------------------------------------------------------
  message("[createPGX] conforming matrices...")

  ## prune unused samples
  contrasts[contrasts %in% c("", " ", "NA")] <- NA
  used.samples <- names(which(rowSums(!is.na(contrasts)) > 0))
  if (prune.samples && length(used.samples) < ncol(counts)) {
    counts <- counts[, used.samples, drop = FALSE]
    samples <- samples[used.samples, , drop = FALSE]
    contrasts <- contrasts[used.samples, , drop = FALSE] ## sample-based!!!
  }

  ## align samples
  ## `counts`, `samples` and `contrasts` span every uploaded sample; `X` spans
  ## the ones outlier removal left (D-24, the column half). The set used to be
  ## intersected with `colnames(X)` first, which cut the dropped samples back
  ## out of `counts` and `samples` -- the same ratchet on the sample axis that
  ## the row trim was on the feature axis, seeded from the same place
  ## (upload_server.R:1093-1095 hands Reanalyse `samples`, `contrasts` and
  ## `counts` off the computed object). `X` is cut to the samples the object has,
  ## keeping whichever of them it still carries.
  kk <- intersect(colnames(counts), rownames(samples))
  counts <- counts[, kk, drop = FALSE]
  X <- X[, intersect(kk, colnames(X)), drop = FALSE]
  samples <- samples[kk, , drop = FALSE]
  samples <- utils::type.convert(samples, as.is = TRUE) ## automatic type conversion
  if (all(kk %in% rownames(contrasts))) {
    contrasts <- contrasts[kk, , drop = FALSE]
  }

  ## sanity checks
  if (ncol(X) == 0) {
    info("[createPGX] FATAL. ncol(X) == 0")
    return(NULL)
  }

  ## -------------------------------------------------------------------
  ## Special case for PTM phospho-proteomics.
  ## -------------------------------------------------------------------
  is.phospho <- annotate_phospho_residue(rownames(counts), detect.only = TRUE)
  if (datatype == "proteomics" && is.phospho) {
    info("[createPGX] annotating rownames with phospho residue...")
    xrows <- counts_rows_of_X(counts, X)
    newnames <- annotate_phospho_residue(rownames(counts))
    newnames <- make_unique(newnames)
    rownames(counts) <- newnames
    rownames(X) <- newnames[xrows]
    if (!is.null(annot_table)) {
      rownames(annot_table) <- newnames
      pos.col <- grep("site|position|phosho", colnames(annot_table), ignore.case = TRUE)
      phosphosite <- sub(".*_|[.].*", "", newnames)
      if (length(pos.col)) {
        i <- pos.col[1]
        annot_table[, i] <- phosphosite
      } else {
        annot_table$site_position <- phosphosite
      }
    }
  }

  ## -------------------------------------------------------------------
  ## create pgx object
  ## -------------------------------------------------------------------
  message("[createPGX] creating pgx object...")

  ## remove special characters from description (other columns too??)
  description <- gsub("[\"\']", " ", description) ## remove quotes (important!!)
  description <- gsub("[\n]", ". ", description) ## replace newline
  description <- trimws(gsub("[ ]+", " ", description)) ## remove ws

  ## add to setting info
  settings$filter.genes <- filter.genes
  settings$exclude.genes <- exclude.genes
  settings$only.known <- only.known
  settings$only.proteincoding <- only.proteincoding
  settings$convert.hugo <- convert.hugo
  settings$custom.geneset <- !is.null(custom.geneset)

  ## add versions info
  versions <- list()
  versions$playbase_version <- packageVersion("playbase")
  versions$playdata_version <- packageVersion("playdata")

  ## Reorder uniprots
  if (datatype == "proteomics") {
    message("[pgx.createPGX] Reordering uniprots in counts, X, annot_table")
    feature.lengths <- NULL
    if (!is.null(annot_table)) {
      kk <- grep("length|size", tolower(colnames(annot_table)))
      if (length(kk) > 0) feature.lengths <- annot_table[, kk[1]]
    }
    ## Taken before the loop, because the loop is what makes the two namings
    ## diverge; the rename itself is row by row, so row i keeps meaning row i.
    xrows <- counts_rows_of_X(counts, X)
    for (i in 1:nrow(counts)) {
      rownames(counts)[i] <- reorder_uniprots(rownames(counts)[i], feature.lengths[i])$feature
    }
    rownames(X) <- rownames(counts)[xrows]
  }
  if (!is.null(annot_table)) rownames(annot_table) <- rownames(counts)

  pgx <- list(
    name = name,
    organism = organism,
    version = packageVersion("playbase"), # useless, just keep for back compatibility
    date = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    creator = creator,
    datatype = datatype,
    datatype_subtype = datatype_subtype,
    description = description,
    metadata = metadata,
    ortholog_species = ortholog_species,
    samples = data.frame(samples, check.names = FALSE),
    counts = as.matrix(counts),
    contrasts = contrasts,
    X = X,
    norm_method = norm_method,
    total_counts = Matrix::colSums(counts, na.rm = TRUE),
    counts_multiplier = counts_multiplier,
    covariates = covariates,
    dma = dma,
    settings = settings,
    versions = versions,
    sc_compute_settings = sc_compute_settings
  )

  ## Create gene annotation table
  pgx$genes <- NULL
  pgx$probe_type <- probe_type

  message("[createPGX] annotating genes")
  pgx$genes <- getProbeAnnotation(
    organism = pgx$organism,
    probes = rownames(pgx$counts),
    datatype = pgx$datatype,
    meth_type = meth_type,
    probetype = pgx$probe_type,
    ortholog_species = pgx$ortholog_species,    
    annot_table = annot_table
  )

  ## Reorder uniprots in pgx$genes. Valid for all datatypes.
  message("[pgx.createPGX] Reordering uniprot column in pgx$genes")
  hh <- grep("uniprot", tolower(colnames(pgx$genes)))
  if (length(hh) > 0) {
    feature.lengths <- NULL
    kk <- grep("length|size", tolower(colnames(pgx$genes)))
    if (length(kk) > 0) feature.lengths <- as.character(pgx$genes[, kk[1]])
    for (i in 1:nrow(pgx$genes)) {
      pgx$genes[i, hh[1]] <- reorder_uniprots(pgx$genes[i, hh[1]], feature.lengths[i])$feature
    }
  }

  if (is.null(pgx$genes)) stop("[pgx.createPGX] FATAL: Could not build gene annotation")

  if (!"symbol" %in% colnames(pgx$genes) && "gene_name" %in% colnames(pgx$genes)) {
    dbg("[pgx.createPGX] WARNING! no symbol column. copying deprecated gene_name column as symbol")
    pgx$genes$symbol <- pgx$genes$gene_name
  }

  if (all(is.na(pgx$genes$symbol))) {
    dbg("[pgx.createPGX] WARNING! all symbol NA. copying rownames as symbol")
    pgx$genes$symbol <- gsub(".*:|[.].*", "", rownames(pgx$genes))
  }

  ## -------------------------------------------------------------------
  ## Filter out not-expressed
  ## -------------------------------------------------------------------
  ## THE ANNOTATION AXIS -- read once, for the three filters below.
  ##
  ## D-24/D-39 stop the PREPROCESSING removals from cutting `counts`: they
  ## record the index they kept and `X` alone shrinks. The three filters that
  ## follow -- here, the only.known/only.proteincoding/exclude.genes block, and
  ## methylomics remove.xy.probes -- do NOT work that way. They decide on
  ## `pgx$genes`, which is built from `rownames(pgx$counts)` just above and is
  ## counts-shaped, and they cut `counts` and `genes` together while `X` follows
  ## through x_rows_for_counts_rows(). So `pgx$counts` leaving pgx.createPGX()
  ## is the upload MINUS whatever these three removed, not the upload, and
  ## upload_server.R:1095 seeds the next Reanalyse from it.
  ##
  ## Measured (playbase.preprocess-5xx, review-c1-evidence/gates-on-ratchet.R):
  ## with the gates on, a Reanalyse round 2 is not round 1 -- max |dX| 0.2368
  ## log2 where only.proteincoding fires, because dropping features moves the
  ## library size CPM divides by. It does not compound: rounds 3 and 4 are
  ## byte-identical to round 2 and `rownames(X)` never moves, because all three
  ## filters are idempotent on a fixed feature set. A one-step shift, not a
  ## ratchet -- but round 1 is still not reproducible from a Reanalyse.
  ##
  ## Closing it is not a matter of deleting the counts cut: `pgx$genes` is what
  ## pgx.add_GMT() builds the gene-set universe from, so keeping `genes` whole
  ## would let genesets survive on members no longer analysed, and shrinking
  ## `genes` below `counts` would break the counts-parallel reading that
  ## pgx.filterZeroCounts() and 166 other `pgx$genes` sites depend on. That
  ## decision is the bead's, not this comment's.
  if (filter.genes) {
    nexpr <- sum(rowSums(pgx$counts, na.rm = TRUE) == 0)
    message("[pgx.createPGX] Filtering out ", nexpr, " not-expressed genes...")
    pgx <- pgx.filterZeroCounts(pgx)
    ii <- match(rownames(pgx$counts), rownames(pgx$genes))
    pgx$genes <- pgx$genes[ii, , drop = FALSE]
  }

  ## -------------------------------------------------------------------
  ## Filter genes
  ## -------------------------------------------------------------------
  do.filter <- (only.known || only.proteincoding || !is.null(exclude.genes))
  if (do.filter) {
    if (only.known) {
      message("[pgx.createPGX] Removing genes without symbol...")
      no.symbol <- (is.na(pgx$genes$symbol) | pgx$genes$symbol %in% c("", "-"))
      pgx$genes <- pgx$genes[which(!no.symbol), , drop = FALSE]
    }

    if (only.proteincoding) {
      message("[pgx.createPGX] Removing Rik/ORF/LOC genes...")
      is.unknown <- grepl("^rik|^loc|^orf", tolower(pgx$genes$symbol))
      is.unknown <- is.unknown & !is.na(pgx$genes$symbol)
      pgx$genes <- pgx$genes[which(!is.unknown), , drop = FALSE]
    }

    if (!is.null(exclude.genes)) {
      message("[pgx.createPGX] excluding genes: ", exclude.genes)
      exstr <- strsplit(tolower(exclude.genes), split = "[ ,]")[[1]]
      exexpr <- paste(c(paste0("^", exstr), paste0(exstr, "$")), collapse = "|")
      exgene <- grepl(exexpr, tolower(pgx$genes$symbol))
      if (sum(exgene)) pgx$genes <- pgx$genes[which(!exgene), , drop = FALSE]
    }

    ## `genes` is counts-shaped while `X` keeps only the rows the preprocessing
    ## removals left (D-24), so `X` cannot be indexed by `genes`' rownames. `X`
    ## goes first, because the translation is an index into the `counts` the
    ## alignment was derived against. The second of the three cuts on the
    ## annotation axis -- see the note above `filter.genes`.
    keep <- match(rownames(pgx$genes), rownames(pgx$counts))
    pgx$X <- pgx$X[x_rows_for_counts_rows(pgx, keep), , drop = FALSE]
    pgx$counts <- pgx$counts[rownames(pgx$genes), , drop = FALSE]
  }

  ## Methylomics arrays: if user-specified, remove X- & Y-linked CpG probes.
  if (pgx$datatype == "methylomics" & remove.xy.probes) {
    kk <- intersect(c("chr", "map"), colnames(pgx$genes))[1]
    if (length(kk) > 0) {
      jj <- grep("chrX|chrY|^X|^Y", pgx$genes[, kk], ignore.case = TRUE)
      if (length(jj) > 0) {
        message("[pgx.createPGX] Methylomics: removing ", length(jj), " X- & Y-linked CpG probes...")
        ## `jj` indexes `genes`, which is `counts`-shaped; `X` may be smaller
        ## (D-24), so it drops the rows whose counts row is going. The third and
        ## last cut on the annotation axis -- see the note above `filter.genes`.
        keep <- setdiff(seq_len(nrow(pgx$counts)), jj)
        pgx$X <- pgx$X[x_rows_for_counts_rows(pgx, keep), , drop = FALSE]
        pgx$counts <- pgx$counts[-jj, , drop = FALSE]
        pgx$genes <- pgx$genes[-jj, , drop = FALSE]
      }
    }
  }

  ## -------------------------------------------------------------------
  ## collapse probe-IDs to gene symbol and aggregate duplicates
  ## -------------------------------------------------------------------
  ## if feature/rownames are not symbol, we paste symbol to row name.
  pp <- sub("^[a-zA-Z]+:", "", rownames(pgx$genes))
  mean_feature_is_symbol <- mean(pp == pgx$genes$symbol, na.rm = TRUE)

  ## NOTE: this was old chunk to convert rownames to HUGO gene
  ## symbol. It now serves to append symbol to rownames/feature names.
  if (convert.hugo && mean_feature_is_symbol < 0.10) {
    symbol <- pgx$genes$symbol
    symbol[is.na(symbol)] <- ""
    feature_is_symbol <- (sub("^[a-zA-Z]+:", "", rownames(pgx$genes)) == symbol)

    new.names <- combine_feature_names(pgx$genes, target = c("rownames", "_", "symbol"))
    new.names <- ifelse(feature_is_symbol, rownames(pgx$genes), new.names)
    new.names <- make_unique(new.names)

    xrows <- counts_rows_of_X(pgx$counts, pgx$X)
    rownames(pgx$genes) <- new.names
    pgx$genes$gene_name <- new.names ## gene_name should also be renamed??
    pgx$genes$feature <- new.names ## feature should also be renamed??
    rownames(pgx$counts) <- new.names
    ## One naming again: `new.names` is as long as `genes`/`counts`, and `X`
    ## takes the entries belonging to the rows it kept (D-24).
    rownames(pgx$X) <- new.names[xrows]
  }

  ## -------------------------------------------------------------------
  ## Infer cell cycle/gender
  ## -------------------------------------------------------------------
  ## This runs above the correction gate, which is why it used to see a
  ## different matrix from the `compute_extra()` call that scores the same two
  ## signatures. That split was an accident of call placement, not a design
  ## (D-07), and it is gone: both sites take `pgx.countScaleMatrix()`, which is
  ## the function's own default, so neither has to name a matrix at all.
  info("[createPGX] infer cell cycle")
  pgx <- compute_cellcycle_gender(pgx)

  ## -------------------------------------------------------------------
  ## Add GMT
  ## -------------------------------------------------------------------
  ## If no organism, no custom annotation table and no custom geneset,
  ## then create empty GMT
  unknown.organism <- (pgx$organism %in% c("No organism", "custom", "unkown"))
  unknown.datatype <- (pgx$datatype %in% c("custom", "unkown"))
  no3 <- unknown.organism && is.null(annot_table) && is.null(custom.geneset)
  if (no3 || unknown.datatype || !add.gmt) {
    message("[pgx.createPGX] WARNING: empty GMT matrix. No gene sets. ")
    pgx$GMT <- Matrix::Matrix(0, nrow = 0, ncol = 0, sparse = TRUE)
  } else {
    pgx <- pgx.add_GMT(
      pgx = pgx,
      custom.geneset = custom.geneset,
      max.genesets = max.genesets,
      include_default_gmt = include_default_gmt,
      species_go = species_go
    )
  }

  ## --------------------------------
  ## rm NA contrasts
  ## --------------------------------
  if (ncol(pgx$samples) > 1) {
    pgx$samples <- pgx$samples[, colMeans(is.na(pgx$samples)) < 1, drop = FALSE]
  }

  ## -------------------------------------------------------------------
  ## Batch correction if user-selected
  ## -------------------------------------------------------------------
  if (batch.correct.method != "no_batch_correct" && ncol(pgx$X) > 2) {
    batch <- NULL
    mm <- batch.correct.method[1]
    if (length(batch.pars) == 0) batch.pars <- "<autodetect>"
    ## Correction operates on `X`, so its covariates are `X`'s samples, not the
    ## object's. Outlier removal leaves the two different (D-24) and every
    ## method underneath is positional in X's column axis, so a full-length
    ## batch factor would be silently paired with the wrong columns. The gate
    ## above counts fittable samples for the same reason.
    ss <- colnames(pgx$X)
    X <- pgx$X
    samples <- pgx$samples[ss, , drop = FALSE]
    contrasts <- pgx$contrasts[ss, , drop = FALSE]

    message("[pgx.createPGX] batch.correct.method=", batch.correct.method)
    message("[pgx.createPGX] batch.pars=", batch.pars)

    pars <- playbase::get_model_parameters(X, samples, pheno = NULL, contrasts)
    if (any(grepl("<autodetect>", batch.pars))) batch.pars <- pars$batch.pars
    if (any(grepl("<none>", batch.pars))) batch.pars <- ""
    batch.pars <- intersect(batch.pars, colnames(samples))
    if (length(batch.pars)) batch <- samples[, batch.pars, drop = FALSE]
    pheno <- pars$pheno

    message("[pgx.createPGX] Batch correction using ", mm)
    if (sum(is.na(X)) == 0) {
      xlist <- playbase::runBatchCorrectionMethods(X, batch, pheno, methods = mm, ntop = Inf)
      cX <- xlist[[mm]]
    } else {
      impute_method <- "SVD2"
      pgx$impute_method <- impute_method ## recorded for the AI-report methods block
      is.mox <- is.multiomics(rownames(X))
      if (is.mox) {
        impX <- imputeMissing.mox(X, method = impute_method)
      } else {
        impX <- imputeMissing(X, method = impute_method)
      }
      xlist <- playbase::runBatchCorrectionMethods(impX, batch, pheno, methods = mm, ntop = Inf)
      cX <- xlist[[mm]]
      jj <- which(is.na(X), arr.ind = TRUE)
      cX[jj] <- NA ## Batch corrected X; original NAs restored
    }

    message("[pgx.createPGX] Batch correction completed\n")

    ## Correction changes X and nothing else. pgx$counts stays the matrix the
    ## user uploaded (playbase-lh8); the count-scale matrix the negative
    ## binomial fitters need is derived from X at the fitter boundary in
    ## compute_testGenes(), and is never persisted.
    pgx$X <- cX

    ## Recorded because this correction happens OUTSIDE the preprocessing record
    ## (D-37): pgx.ranWithCorrection() reads that record and so cannot see it,
    ## which left pgx.countScaleMatrix() answering "the upload" on an object
    ## whose `X` had moved, with nothing to say so. It is set only where the
    ## correction actually ran, never from the argument alone -- the gate above
    ## declines on fewer than three fittable samples.
    pgx$settings$batch.correct.method <- mm

    rm(xlist, cX)
  }

  rm(counts, X, samples, contrasts)

  message("\n\n")
  message("[pgx.createPGX]======================================")
  message("[pgx.createPGX]======== pgx.createPGX: DONE! ========")
  message("[pgx.createPGX]======================================")
  message("\n\n")

  return(pgx)
}


#' @title Compute PGX
#' @description Main function to populate pgx with results. The function computes the analysis on a pgx object
#'
#' @param pgx A pgx object containing the input data
#' @param max.genes Maximum number of genes to test. Default is 19999.
#' @param gx.methods Methods for differential expression analysis at the gene level. Default is c("ttest.welch", "trend.limma", "edger.qlf").
#' @param gset.methods Methods for differential analysis at the gene set level. Default is c("fisher", "gsva", "fgsea").
#' @param do.cluster Logical indicating whether to run sample clustering. Default is TRUE.
#' @param do.clustergenesets Logical indicating whether to cluster gene sets.
#' @param do.clustergenes Logical indicating whether to cluster genes. Default is TRUE.
#' @param use.design Whether to use model design matrix for testing. Default is FALSE.
#' @param prune.samples Whether to remove samples without valid contrasts. Default is TRUE.
#' @param time Whether perform time series analysis or not. Default FALSE
#' @param extra.methods Additional analysis methods to run. Default is c("meta.go", "infer", "deconv", "drugs", "wordcloud", "wgcna")[c(1, 2)].
#' @param libx.dir Directory containing custom analysis modules.
#' @param progress A progress object for tracking status.
#' @param ai_features Optional list of AI features to run after core compute.
#'
#' @details
#' The slots created by pgx.computePGX are the following:
#'
#' - `tsne2d`: 2D tSNE coordinates matrix
#' - `tsne3d`: 3D tSNE coordinates matrix
#' - `cluster`: List containing sample clustering results
#' - `cluster.genes`: List containing gene clustering results
#' - `model.parameters`: Model parameters from normalization
#' - `timings`: Matrix of timings for computations
#' - `gx.meta`: Gene metadata data.frame
#' - `gset.meta`: Gene set metadata data.frame
#' - `gsetX`: Gene set scores matrix
#' - `cluster.gsets`: List of gene set clustering results
#' - `meta.go`: GO graph and metadata
#'
#' @return An updated pgx object containing analysis results.
#'
#' @export
pgx.computePGX <- function(pgx,
                           max.genes = 19999,
                           gx.methods = c("trend.limma", "edger.qlf", "deseq2.wald"),
                           gset.methods = c("fisher", "gsva", "fgsea"),
                           custom.geneset = list(gmt = NULL, info = NULL),
                           custom_fc = NULL,
                           do.cluster = TRUE,
                           cluster.contrasts = FALSE,
                           do.clustergenesets = TRUE,
                           do.clustergenes = TRUE,
                           use.design = FALSE,
                           prune.samples = TRUE,
                           extra.methods = c(
                             "meta.go", "infer", "deconv", "drugs",
                             "connectivity", "wordcloud", "wgcna",
                             "mofa"
                           )[c(1, 2)],
                           pgx.dir = NULL,
                           libx.dir = NULL,
                           progress = NULL,
                           ai_features = NULL,
                           user_input_dir = getwd()) {
  message("[pgx.computePGX]===========================================")
  message("[pgx.computePGX]========== pgx.computePGX =================")
  message("[pgx.computePGX]===========================================")
  message("\n")
  message("[pgx.computePGX] Starting pgx.computePGX")
  message("\n")

  if (!"contrasts" %in% names(pgx)) {
    stop("[pgx.computePGX] FATAL:: no contrasts in object")
  }

  if (!all(grepl("_vs_", colnames(pgx$contrasts)))) {
    stop("[pgx.computePGX] FATAL:: all contrast names must include _vs_")
  }

  ## -----------------------------------------------------------------------------
  ## Time series: check methods
  ## -----------------------------------------------------------------------------
  timeseries <- any(grepl("IA:*", colnames(pgx$contrasts)))
  if (timeseries) {
    ts.mm <- c("trend.limma", "deseq2.lrt", "deseq2.wald", "edger.lrt", "edger.qlf")
    cm <- intersect(gx.methods, ts.mm)
    if (length(cm) == 0) {
      message(
        "[pgx.computePGX] For time series analysis, gx.methods must be among ",
        paste0(ts.mm, collapse = "; "), " Skipping time series analysis."
      )
      hh <- grep("IA:*", colnames(pgx$contrasts))
      pgx$contrasts <- pgx$contrasts[, -hh, drop = FALSE]
    } else {
      gx.methods <- cm
    }
  }

  contr.matrix <- contrasts.convertToLabelMatrix(pgx$contrasts, pgx$samples)
  contr.matrix <- makeContrastsFromLabelMatrix(contr.matrix)
  contr.matrix <- sign(contr.matrix) ## sign is fine

  ## The design is built from `pgx$samples`, which names every uploaded sample;
  ## `X` names the ones outlier removal left (D-24, the column half). A sample
  ## with no expression data cannot be fitted, so it leaves the design here --
  ## once, where the design is built, rather than at each of its readers. The
  ## contrast pruning below then drops any comparison that emptied out.
  ##
  ## `ss` is ordered by `X`, and the subset is taken whether or not anything was
  ## dropped, so the design's rows ARE `X`'s columns, in `X`'s order. Readers
  ## that walk the design and the expression matrix together need that -- they
  ## had it by coincidence, because the sample table happened to be in the same
  ## order as `X`; `compute_testGenes()` copies these rows into
  ## `model.parameters$exp.matrix` untouched, and that is what they read.
  ss <- intersect(colnames(pgx$X), rownames(contr.matrix))
  if (length(ss) < nrow(contr.matrix)) {
    message(
      "[pgx.computePGX] ", nrow(contr.matrix) - length(ss),
      " sample(s) are not in X and cannot be tested; excluded from the design"
    )
  }
  contr.matrix <- contr.matrix[ss, , drop = FALSE]
  if (!identical(rownames(contr.matrix), colnames(pgx$X))) {
    stop(
      "[pgx.computePGX] the design and X disagree on samples: ",
      "every column of X must be in pgx$samples for a design row to exist for it"
    )
  }

  ## sanity check
  if (NCOL(contr.matrix) == 0) {
    message("[pgx.computePGX] WARNING: FATAL ERROR. zero contrasts")
    return(pgx)
  }

  ## select valid contrasts
  sel <- Matrix::colSums(contr.matrix == -1) > 0 & Matrix::colSums(contr.matrix == 1) > 0
  contr.matrix <- contr.matrix[, sel, drop = FALSE]

  ## -------------------------------------------------------------------
  ## Clustering
  ## -------------------------------------------------------------------
  ## Cluster by sample
  if (do.cluster || cluster.contrasts) {
    message("[pgx.computePGX] clustering samples...")
    mm <- c("pca", "tsne", "umap")
    pgx <- pgx.clusterSamples(pgx, dims = c(2, 3), perplexity = NULL, X = NULL, methods = mm)
  }

  ## Make contrasts by cluster
  if (cluster.contrasts) {
    ## NEED RETHINK: for the moment we use combination of t-SNE/UMAP
    posx <- cbind(pgx$cluster$pos[["umap2d"]], pgx$cluster$pos[["tsne2d"]])
    posx <- scale(posx)
    idx <- pgx.findLouvainClusters(posx, level = 1, prefix = "c", small.zero = 0.0)
    if (length(unique(idx)) == 1) {
      ## try again with finer settings if single cluster...
      idx <- pgx.findLouvainClusters(posx, level = 2, prefix = "c", small.zero = 0.01)
    }
    pgx$samples$.cluster <- idx ## really add??

    ## Add cluster contrasts
    message("[pgx.computePGX] adding cluster contrasts...")
    Y <- pgx$samples[, ".cluster", drop = FALSE]
    if (length(unique(Y[, 1])) < 2) {
      message("[pgx.computePGX] warning: only one cluster.")
    } else {
      ct <- makeDirectContrasts(Y, ref = "others")
      ctx <- contrastAsLabels(ct$exp.matrix)
      if (ncol(pgx$contrasts) == 0) {
        pgx$contrasts <- ctx
      } else {
        pgx$contrasts <- cbind(pgx$contrasts, ctx)
      }
    }
  }

  ## Cluster by genes
  if (do.clustergenes) {
    message("[pgx.computePGX] clustering genes...")
    mm <- "umap"
    if (pgx$datatype == "scRNAseq") mm <- c("pca", "tsne", "umap")
    pgx <- pgx.clusterGenes(pgx, methods = mm, level = "gene")
  }

  ## Shrink number of genes (highest SD/var)
  ##
  ## The shrink keeps the top `max.genes` features of `X` by the standard
  ## deviation of their log-CPM, and the matrix it ranks is conditional -- it
  ## always was (D-06). Legacy read `pgx$counts`, which by the time it got here
  ## had been cut to `X`'s features and samples, and had been OVERWRITTEN by the
  ## batch-corrected reconstruction whenever correction ran. The cut is gone
  ## (D-39) and the overwrite is gone (playbase-lh8), so both halves have to be
  ## said out loud instead of arriving by side effect:
  ##
  ##   universe -- what `X` is made of, not what was uploaded, so `counts` is
  ##     aligned down to `X` first (D-24).
  ##   values   -- the reconstruction when correction ran, the pristine upload
  ##     when it did not. This is `pgx.countScaleMatrix()`'s rule (D-07), spelt
  ##     out rather than called: the universe clause above applies to one branch
  ##     only, because the reconstruction already has `X`'s shape. A ranking
  ##     universe is the shrink's own question, and the one place the shared
  ##     rule is deliberately not the whole answer.
  ##
  ## `pgx.removeLowVariance()` carries this same gate and returns before it
  ## aligns anything. It is repeated here because `rank_on` is built eagerly,
  ## and building it on an object the gate would have skipped can refuse where
  ## legacy did nothing at all.
  if (max.genes > 0 && nrow(pgx$X) > max.genes) {
    rank_on <- if (isTRUE(pgx.ranWithCorrection(pgx))) {
      pgx.recomputeCounts(pgx)
    } else {
      a <- counts_index_of_X(pgx$counts, pgx$X)
      pgx$counts[a$rows, a$cols, drop = FALSE]
    }
    ## The verb appends a provenance record, and that record does not live on
    ## the objects playbase holds (D-37). A one-step chain anchored at `X`'s
    ## shape would make every later pgx.alignXtoCounts() refuse, and would turn
    ## pgx.ranWithCorrection()'s NA -- "unknown" -- into FALSE.
    pp <- pgx$settings$preprocessing
    pgx <- pgx.removeLowVariance(pgx, n = max.genes, rank_on = rank_on)
    pgx$settings$preprocessing <- pp
  }

  ## `counts` is NOT re-cut to `X` here any more (D-42). That re-cut undid D-39
  ## on every compute run, and `upload_server.R:1095` seeds Reanalyse from a
  ## computed object, so it put the ratchet back by a second route.

  pgx$timings <- c()
  GENETEST.METHODS <- c(
    "ttest", "ttest.welch", "ttest.rank",
    "voom.limma", "trend.limma", "notrend.limma",
    "edger.qlf", "edger.lrt", "deseq2.wald", "deseq2.lrt"
  )
  GENESETTEST.METHODS <- c(
    "fisher", "gsva", "ssgsea", "spearman",
    "camera", "fry", "fgsea"
  ) ## no GSEA, too slow...

  ## ------------------ gene level tests ---------------------
  if (!is.null(progress)) progress$inc(0.1, detail = "testing genes")

  timeseries <- any(grepl("^IA:*", colnames(pgx$contrasts)))

  message("[pgx.computePGX] testing genes...")
  pgx <- compute_testGenes(
    pgx,
    contr.matrix,
    max.features = max.genes,
    test.methods = gx.methods,
    custom_fc = custom_fc,
    ## use.design = use.design,
    prune.samples = prune.samples,
    timeseries = timeseries,
    remove.outputs = TRUE
  )

  ## ------------------ gene set tests -----------------------
  if (!is.null(progress)) progress$inc(0.2, detail = "testing gene sets")

  if ((pgx$organism != "No organism" && !is.null(pgx$GMT) && nrow(pgx$GMT) > 0) ||
    (pgx$organism == "No organism" && !is.null(custom.geneset$gmt))) {
    message("[pgx.computePGX] testing genesets...")

    pgx <- compute_testGenesets(
      pgx = pgx,
      custom.geneset = custom.geneset,
      test.methods = gset.methods,
      use.replaid = TRUE
    )

    ## Cluster by genes
    if (do.clustergenesets) {
      message("[pgx.computePGX] clustering genesets...")
      pgx <- pgx.clusterGenes(pgx, methods = "umap", X = NULL, level = "geneset")
    }
  } else {
    message("[pgx.computePGX] Skipping genesets test")
  }

  ## ------------------ extra analyses ---------------------
  if (!is.null(progress)) progress$inc(0.3, detail = "extra modules")
  message("[pgx.computePGX] computing extra modules: ", paste0(extra.methods, collapse = "; "))
  pgx <- compute_extra(
    pgx,
    extra = extra.methods,
    pgx.dir = pgx.dir,
    libx.dir = libx.dir,
    user_input_dir = user_input_dir
  )

  ## methylomics: ensure all OPG graphics & tables use beta.
  if (pgx$datatype == "methylomics") pgx$X <- playbase::mToBeta(pgx$X)

  if (!is.null(ai_features)) {
    if (!is.list(ai_features)) {
      stop("[pgx.computePGX] ai_features must be a list", call. = FALSE)
    }
    if (!is.null(ai_features$reports)) {
      message("[pgx.computePGX] generating AI reports...")
      pgx <- pgx.update_reports(pgx, ai = ai_features$reports)
    }
    if (!is.null(ai_features$wgcna_summaries)) {
      info("[pgx.computePGX] generating WGCNA module summaries...")
      pgx <- pgx.update_wgcna_summaries(pgx, ai = ai_features$wgcna_summaries)
    }
    if (!is.null(ai_features$infographics)) {
      info("[pgx.computePGX] generating AI infographics...")
      pgx <- pgx.update_infographics(pgx, ai = ai_features$infographics)
    }
  }

  info("[pgx.computePGX] DONE")
  return(pgx)
}


## ===================================================================
## =================== UTILITY FUNCTIONS =============================
## ===================================================================

#' @export
counts.autoScaling <- function(counts) {
  message("[createPGX] scaling counts...")
  counts_multiplier <- 1

  ## If the difference in total counts is too large, we need to
  ## euqalize them because the thresholds can become strange. Here
  ## we decide if normalizing is necessary (WARNING changes total
  ## counts!!!)
  totcounts <- Matrix::colSums(counts, na.rm = TRUE)
  totratio <- log10(max(1 + totcounts, na.rm = TRUE) / min(1 + totcounts, na.rm = TRUE))
  totratio

  if (totratio > 6) {
    message("[createPGX:autoscale] WARNING: too large total counts ratio. forcing normalization.")
    meancounts <- exp(mean(log(1 + totcounts), na.rm = TRUE))
    counts <- t(t(counts) / totcounts) * meancounts
  }

  ## Check if too big (more than billion reads). This is important
  ## for some proteomics intensity signals that are in billions of
  ## units.
  mean.counts <- mean(Matrix::colSums(counts, na.rm = TRUE))
  is.toobig <- log10(mean.counts) > 9
  if (is.toobig) {
    ## scale to about 10 million reads
    message("[createPGX:autoscale] WARNING: too large total counts. Scaling down to 10e6 reads.")
    unit <- 10**(round(log10(mean.counts)) - 7)
    unit
    counts <- counts / unit
    counts_multiplier <- unit
  }
  counts_multiplier
  message("[createPGX:autoscale] count_multiplier= ", counts_multiplier)

  list(counts = counts, counts_multiplier = counts_multiplier)
}

#' @export
counts.mergeDuplicateFeatures <- function(counts, is.counts = TRUE) {
  counts <- counts[rownames(counts) != "", ]
  counts[which(is.nan(counts))] <- NA
  ndup <- sum(duplicated(rownames(counts)))
  if (ndup > 0) {
    if (!is.counts) counts <- 2**counts
    message("[mergeDuplicateFeatures] ", ndup, " duplicated rownames: averaging rows (in counts).")
    counts <- playbase::rowmean(counts, group = rownames(counts), reorder = TRUE)
    counts[which(is.nan(counts))] <- NA
    if (!is.counts) counts <- log2(counts)
  }
  counts
}

#' @export
pgx.filterZeroCounts <- function(pgx) {
  ## There is second filter in the statistics computation. This
  ## first filter is primarily to reduce the counts table.
  ## AZ: added na.rm=TRUE to avoid introducing NAs and edit to keep NAs.
  keep <- (Matrix::rowMeans(pgx$counts > 0, na.rm = TRUE) > 0) ## at least in one...

  ## Positional, and the all-NA rows still first: `counts` is counts-shaped
  ## while `X` keeps only the rows the preprocessing removals left (D-24), so
  ## one set of rownames no longer indexes both. `genes` follows `counts`; `X`
  ## drops the rows whose counts row is going, through the alignment. This is
  ## the first of the three cuts on the annotation axis -- see the note above
  ## the `filter.genes` gate in pgx.createPGX() for what that costs.
  keep <- c(which(is.na(keep)), which(keep))
  pgx$X <- pgx$X[x_rows_for_counts_rows(pgx, keep), , drop = FALSE]
  pgx$counts <- pgx$counts[keep, , drop = FALSE]
  pgx$genes <- pgx$genes[keep, , drop = FALSE]

  pgx
}

#' @export
pgx.filterLowExpressed <- function(pgx, prior.cpm = 1) {
  AT.LEAST <- ceiling(pmax(2, 0.01 * ncol(pgx$counts)))
  message("filtering for low-expressed genes: > ", prior.cpm, " CPM in >= ", AT.LEAST, " samples")
  keep <- (rowSums(edgeR::cpm(pgx$counts) > prior.cpm, na.rm = TRUE) >= AT.LEAST)
  pgx$filtered <- NULL
  pgx$filtered[["low.expressed"]] <- paste(rownames(pgx$counts)[which(!keep)], collapse = ";")
  if (!is.null(pgx$X)) {
    ## Before `counts` is cut, because the alignment is an index into the
    ## matrix the chain ran against. The warning that stood here -- "counts and
    ## X should match dimensions" -- was the assumption, not a check: `keep` is
    ## positional in `counts`, and under D-24 `counts` is the whole upload while
    ## `X` is what the removals left.
    pgx$X <- pgx$X[x_rows_for_counts_rows(pgx, which(keep)), , drop = FALSE]
  }
  pgx$counts <- pgx$counts[keep, , drop = FALSE]
  ## `genes` follows `counts`, positionally, exactly as it does in
  ## pgx.filterZeroCounts(). Leaving it whole while `counts` shrank put the two
  ## out of step on the one axis the rest of the object indexes them by.
  if (!is.null(pgx$genes)) pgx$genes <- pgx$genes[keep, , drop = FALSE]
  message("filtering out ", sum(!keep), " low-expressed genes")
  message("keeping ", sum(keep), " expressed genes")
  pgx
}


#' Internal use: append gmt (list) to a sparse gene set matrix
#'
.append_gmt_to_matrix <- function(gmt, G, all_genes, minsize, maxsize) {
  if (is.null(all_genes)) {
    all_genes <- unique(unlist(gmt))
    if (!is.null(G)) all_genes <- unique(c(rownames(G), all_genes))
  }

  ## Only admmitted if valid symbol
  gmt <- lapply(gmt, function(s) intersect(s, all_genes))
  gmt.size <- sapply(gmt, length)
  if (sum(gmt.size >= minsize & gmt.size <= maxsize) == 0) {
    message("[.append_gmt_to_matrix] warning no valid gmt to add")
    return(G)
  }
  
  ## check if we have the new prefixed GMT format. If so we need to
  ## match to stripped names and prepend matching datatype prefix
  if( !is.null(G) && mean(grepl("[:]",rownames(G))) > 0.5) {
    gmt_genes <- unique(unlist(gmt))
    names(gmt_genes) <- gmt_genes
    stripped_G_names <- sub("^[a-zA-Z]+:","",rownames(G))
    jj <- match(gmt_genes, stripped_G_names) 
    ii <- which(!is.na(jj))
    if(length(ii)) {
      gmt_genes[ii] <- rownames(G)[jj[ii]]
      gmt <- lapply(gmt, function(m) unname(gmt_genes[m]))
    }
  }

  add_gmt <- createSparseGenesetMatrix(
    gmt.all = gmt,
    min.geneset.size = minsize,
    max.geneset.size = maxsize,
    min_gene_frequency = 1,
    all_genes = all_genes,
    annot = NULL,
    filter_genes = FALSE
  )

  # G and custom_gmt have to be SYMBOL alligned
  if (!is.null(add_gmt) && ncol(add_gmt) > 0) {
    G <- merge_sparse_matrix(G, Matrix::t(add_gmt))
    remove(add_gmt)
  }
  return(G)
}

pgx.add_GMT <- function(pgx,
                        custom.geneset = NULL,
                        max.genesets = 20000,
                        include_default_gmt = TRUE,
                        include_iea = TRUE, species_go = NULL) {
  ## An explicit NULL (e.g. params$include_default_gmt from a params.RData
  ## written before this parameter existed) bypasses the TRUE default
  ## above, since R only applies argument defaults when the argument is
  ## missing, not when it's passed as NULL.
  if (is.null(include_default_gmt)) include_default_gmt <- TRUE

  if (!"symbol" %in% colnames(pgx$genes)) {
    message(paste(
      "[pgx.add_GMT] ERROR: could not find 'symbol' column.",
      "Is this an old gene annotation?"
    ))
    return(pgx)
  }

  ## -----------------------------------------------------------
  ## Load Geneset matrix and filter genes by gene or homologous
  ## -----------------------------------------------------------
  message("[pgx.add_GMT] Creating GMT matrix... ")

  # Load geneset matrix from playdata. add metabolomics if data.type
  # is metabolomics. GSETxGENE is keyed on human gene symbols, so
  # prefer the guaranteed-human ortholog column over the
  # species-configurable one.
  target <- c("human_ortholog", "ortholog", "symbol", "gene_name", "rownames")
  ortho.col <- intersect(target, colnames(pgx$genes))
  if (length(ortho.col) == 0) {
    symbol <- toupper(pgx$genes$symbol)
  } else {
    symbol <- pgx$genes[, ortho.col[1]] ## human symbol!
  }

  ## check if we have genes/proteins
  symbol <- sub(".*:", "", symbol) ## strip prefix
  sum.px <- sum(symbol %in% colnames(playdata::GSETxGENE), na.rm = TRUE)
  has.px <- sum.px >= 10

  ## check if we have metabolites/lipids
  has.mx1 <- grepl("metabolomics|lipidomics", pgx$datatype, ignore.case = TRUE)
  has.mx2 <- pgx$datatype == "multi-omics" && any(grepl("mx|metabolomics|lipidomics", pgx$genes$data_type))
  has.mx3 <- pgx$datatype == "multi-omics" && !all(grepl("mx|metabolomics|lipidomics", pgx$genes$data_type))

  has.mx <- has.mx1 || has.mx2
  has.px2 <- !has.mx1 || has.mx3

  dbg("[pgx.add_GMT] 1: has.px = ", has.px)
  dbg("[pgx.add_GMT] 1: has.px2 = ", has.px2)
  dbg("[pgx.add_GMT] 1: has.mx = ", has.mx)

  ## Note!!!: Rownames of G must be in species symbol (not anymore
  ## human ortholog).
  G <- NULL

  ## add metabolomic gene sets
  if (has.mx) {
    info("[pgx.add_GMT] Retrieving metabolomics genesets")
    G <- mx.create_metabolite_sets(
      annot = pgx$genes,
      gmin = 0,
      metmin = 5,
      as_matrix = TRUE
    )
  }

  ## add SYMBOL (classic) gene sets
  if (has.px && include_default_gmt) {
    info("[pgx.add_GMT] Retrieving transcriptomics/proteomics genesets")
    G1 <- Matrix::t(playdata::GSETxGENE)
    G1 <- rename_by2(G1, pgx$genes, new_id = "symbol") ## symbol!
    G <- merge_sparse_matrix(G, G1)
  }

  # create a feature list that will be used to filter and reduce dimensions of G
  ##full_feature_list <- c(pgx$genes$symbol, pgx$genes$ortholog, rownames(pgx$genes)) ## why ?? 
  full_feature_list <- c(pgx$genes$symbol)
  full_feature_list <- setdiff(full_feature_list, c(NA,""))
  full_feature_list <- unique(full_feature_list)

  if (!is.null(G)) {
    G <- G[rownames(G) %in% full_feature_list, , drop = FALSE]
    G <- G[, Matrix::colSums(G != 0) > 0, drop = FALSE]
    if (nrow(G) == 0 || ncol(G) == 0) G <- NULL
  }

  ## Add organism specific GO gene sets. This is species gene
  ## symbol. Skip if the GMT has enough (>1000) terms.
  num_goterms <- sum(grepl("^GO",colnames(G)))
  info("[pgx.add_GMT] number of GO gene sets in GMT =",num_goterms)
  if(is.null(species_go)) {
    ##species_go <- (num_goterms < 1000)
    species_go <- !(pgx$organism %in% c("Human","Mouse","Rat"))
  }

  if (has.px2 && species_go) {
    ## add species GO genesets from AnnotationHub
    info("[pgx.add_GMT] Retrieving species GO for organism", pgx$organism,"...")
    go.main = go.ortho = NULL

    ## Lookup GO for main species
    go.main <- tryCatch({
      getOrganismGO(
        organism = pgx$organism,
        symbol.annot = pgx$genes,
        features = full_feature_list,
        db = c("annothub","gprofiler"),
        include_iea = include_iea)
    }, error = function(e) {
      message("Error in getOrganismsGO:", e)
    })

    ## Lookup GO for ortholog species
    if(!is.null(pgx$ortholog_species)) {
      go.ortho <- tryCatch({
        getOrganismGO(
          organism = pgx$ortholog_species,
          symbol.annot = pgx$genes,
          features = full_feature_list,
          db = c("annothub","gprofiler"),
          include_iea = include_iea)
      }, error = function(e) {
        message("Error in getOrganismsGO:", e)
      })
    }

    ## merge both
    go.genesets <- c(go.main, go.ortho)
    
    if (!is.null(go.genesets)) {
      go.genesets <- go.merge_duplicates(go.genesets)
      dbg("[pgx.add_GMT] Adding", length(go.genesets), "species GO genesets")
      all_genes <- unique(pgx$genes$symbol)      
      G <- .append_gmt_to_matrix(go.genesets, G, all_genes, minsize = 15, maxsize = 400)
    } ## end-if go.genesets
  } ## end-if !metabolics

  ## Add custom gene sets if provided
  if (!is.null(custom.geneset$gmt)) {
    ## convert gmt standard to SPARSE matrix: gset in rows, genes in columns.
    custom_gmt <- custom.geneset$gmt
    custom_gmt <- custom_gmt[sapply(custom_gmt,length)>1]    
    message(paste("[pgx.add_GMT] Adding",length(custom_gmt),"custom genesets"))
    ## Map feature id always to species specific symbols. This uses
    ## the feature annotation table pgx$genes so it also uses the
    ## ortholog columns for matching.
    custom_gmt <- gmt.map2symbol(custom_gmt, annot=pgx$genes, target="symbol") 
    all_genes <- unique(pgx$genes$symbol)
    G <- .append_gmt_to_matrix(custom_gmt, G, all_genes, minsize = 3, maxsize = 9999)
  }

  num_goterms <- sum(grepl("^GO",colnames(G)))
  info("[pgx.add_GMT] total number of GO terms = ",num_goterms)
  
  ## -----------------------------------------------------------
  ##  Prioritize gene sets by fast rank-correlation
  ## -----------------------------------------------------------
  ## NEED RETHINK!! IK. Probably not needed anymore with generalized
  ## features. Generally G and X are not aligned anymore.
  ## !!!!!!!!!!!!
  ## NOTE: this can be replace by PLAID??

  if (is.null(max.genesets)) max.genesets <- 20000
  if (max.genesets < 0) max.genesets <- 20000
  if (!is.null(G) && ncol(G) > max.genesets) {
    message("[pgx.add_GMT] Matching gene set matrix...")
    # we use SYMBOL as rownames
    gX <- pgx$X
    if (!all(rownames(gX) %in% pgx$genes$symbol)) {
      gX <- rename_by(gX, pgx$genes, "symbol", unique = TRUE)
    }

    ## if reduced samples
    ss <- rownames(pgx$model.parameters$exp.matrix)
    if (!is.null(ss)) {
      gX <- gX[, ss, drop = FALSE]
    }

    ## Align the GENESETxGENE matrix with genes in X_geneset
    gg <- rownames(gX)
    ii <- intersect(gg, rownames(G))
    G <- G[ii, , drop = FALSE]
    ## gX <- gX[ii, , drop = FALSE]
    xx <- setdiff(gg, rownames(G))
    matX <- Matrix::Matrix(0, nrow = length(xx), ncol = ncol(G), sparse = TRUE)
    rownames(matX) <- xx
    colnames(matX) <- colnames(G)
    G <- rbind(G, matX)
    G <- G[match(gg, rownames(G)), , drop = FALSE]
    rownames(G) <- rownames(gX) ## must be symbol

    ## Prioritize gene sets by fast rank-correlation
    message("[pgx.add_GMT] Reducing gene set matrix... ")
    ## Reduce gene sets by selecting top varying genesets. We use the
    ## very fast sparse rank-correlation for approximate single sample
    ## geneset activation.
    cX <- gX - rowMeans(gX, na.rm = TRUE) ## center!
    cX <- t(matrixStats::colRanks(cX))
    if (ncol(cX) <= 5000) {
      gsetX <- qlcMatrix::corSparse(G, cX)
    } else { ## split into chuncks. faster & needs less memory.
      index <- unique(c(seq(1, ncol(cX), by = round(ncol(cX) / 10, 0)), ncol(cX)))
      i <- 1
      LL.cor <- list()
      for (i in 1:(length(index) - 1)) {
        if (index[i] == 1) jj <- 1:index[i + 1]
        if (index[i] > 1) jj <- (index[i] + 1):index[i + 1]
        LL.cor[[i]] <- qlcMatrix::corSparse(G, cX[, jj, drop = FALSE])
      }
      gsetX <- do.call(cbind, LL.cor)
      rm(index, LL.cor)
    }
    grp <- pgx$model.parameters$group
    gsetX.bygroup <- NULL
    ## If groups/conditions are present we calculate the SD by group
    if (!is.null(grp)) {
      gsetX.bygroup <- tapply(1:ncol(gsetX), grp, function(i) rowMeans(gsetX[, i, drop = FALSE], na.rm = TRUE))
      gsetX.bygroup <- do.call(cbind, gsetX.bygroup)
      ## sdx <- apply(gsetX.bygroup, 1, stats::sd, na.rm = TRUE)
      sdx <- matrixStats::rowSds(gsetX.bygroup, na.rm = TRUE)
    } else {
      sdx <- matrixStats::rowSds(gsetX, na.rm = TRUE)
    }
    names(sdx) <- colnames(G)
    jj <- Matrix::head(order(-sdx), max.genesets)
    must.include <- "hallmark|kegg|^go|^celltype|^pathway|^custom|^metabo"
    jj <- unique(c(jj, grep(must.include, colnames(G), ignore.case = TRUE)))
    jj <- jj[order(colnames(G)[jj])] ## sort alphabetically
    G <- G[, jj, drop = FALSE]
    rm(gsetX.bygroup, gsetX)
  }

  ## -----------------------------------------------------------------------
  ## Clean up and return pgx object
  ## -----------------------------------------------------------------------

  # final check: drop genesets in G based on geneset size
  if (!is.null(G)) {
    gmt.size <- Matrix::colSums(G != 0)
    has.metabolites <- sum(grepl("^[0-9]+$|CHEBI|LIPID", rownames(G))) >= 10
    has.metabolites
    ## if (pgx$datatype %in% c("metabolomics","multi-omics")) {
    if (has.metabolites) {
      # metabolomics genesets are MUCH smaller than transcriptomics,
      # metabolomics have usually less features, so we need to reduce
      # the min size
      size.ok <- which(gmt.size >= 3 & gmt.size <= 400)
    } else {
      size.ok <- which(gmt.size >= 10 & gmt.size <= 400)
    }

    # add all custom genesets to size.ok
    idx_custom_gmt <- grep("CUSTOM", colnames(G))
    # make sure we dont miss CUSTOM genesets due to size.ok exclusion
    if (length(idx_custom_gmt) > 0) {
      names(idx_custom_gmt) <- colnames(G)[idx_custom_gmt]
      size.ok <- union(size.ok, idx_custom_gmt)
    }
    G <- G[, size.ok, drop = FALSE]
  }

  # add random genesets if G is too small
  if (is.null(G) || ncol(G) < 30 || nrow(G) < 3) {
    add.gmt <- NULL
    rr <- sample(3:400, 50)
    gg <- pgx$genes$symbol
    random.gmt <- lapply(rr, function(n) head(sample(gg), min(n, length(gg) / 2)))
    names(random.gmt) <- paste0("TEST:random_geneset.", 1:length(random.gmt))
    # Extreme low feature count control, avoids crash
    if (all(lapply(random.gmt, length) |> unlist() < 3)) {
      min.geneset.size <- 1
    } else {
      min.geneset.size <- 3
    }

    G <- .append_gmt_to_matrix(
      random.gmt, G,
      all_genes = unique(pgx$genes$symbol),
      minsize = min.geneset.size,
      maxsize = 400
    )
  }

  # normalize columns (required for some methods downstream)log2foldchange
  G <- normalize_cols(G)

  pgx$GMT <- G
  pgx$custom.geneset <- custom.geneset
  message(glue::glue("[pgx.add_GMT] Final GMT: {nrow(G)} x {ncol(G)}"))
  rm(G)

  gc()
  return(pgx)
}


## ----------------------------------------------------------------------
## -------------------------- end of file -------------------------------
## ----------------------------------------------------------------------
