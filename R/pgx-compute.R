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
#'   [pgx.createPGX()] and from there to
#'   [playbase.preprocess::pgx.preprocess()]. NULL (the default)
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
#' counts <- playbase.ingest::example_file("counts.csv")
#' contrasts <- playbase.ingest::example_file("contrasts.csv")
#' samples <- playbase.ingest::example_file("samples.csv")
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
  counts <- playbase.ingest::read.as_matrix(counts.file)

  ## compile sample table
  samples <- playbase.ingest::read.as_matrix(samples.file)
  samples <- data.frame(samples, check.names = FALSE)

  ## parse requested phenotypes
  if (!is.null(contrasts.file) && file.exists(contrasts.file)) {
    message("reading contrasts file ", contrasts.file)
    contrasts <- playbase.ingest::read.as_matrix(contrasts.file)
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

#' Create a PGX object
#' This function creates a pgx object, which is the core object in the
#' OmicsPlayground.
#' @param counts Matrix of count data with genes as rows and samples as columns.
#' @param samples Data frame containing sample information.
#' @param organism Default "Human", it indicates the species used
#' for the gene annotation table and the probe to symbol conversion.
#' @param contrasts Data frame defining sample contrasts.
#' @param X (Optional) Matrix of normalized expression data. If NULL, will be calculated from counts.
#' @param preprocess (Optional) Named list of preprocessing options. Whenever
#'   `X` is NULL, `X` is built from `counts` via
#'   [playbase.preprocess::pgx.preprocess()]. A NULL value selects explicit
#'   legacy-compatible log2 conversion; supplied options select leaf pipeline
#'   stages. Batch correction is selected only by `batch.correct.method` and
#'   runs inside that pipeline. The legacy conversion requires at least one
#'   finite positive count.
#' @param is.logx Deprecated upload hint. Confirmed log2 uploads must be
#'   back-transformed to counts before this function is called; they then pass
#'   through the complete preprocessing pipeline like any other count input.
#' @param dotimeseries Logical indicating if timeseries analysis has been activated by the user at upload
#' @param batch.correct.method BC method. Default is "no_batch_correct" (meaning no batch correction).
#'   For multi-omics, a named vector selects one method per layer prefix,
#'   e.g. `c(gx = "limma", px = "no_batch_correct")`.
#' @param batch.pars BC variable. Default "autodetect" as per QC/BC tab in upload.
#' @param covariates variables to regress out. Valid only for linear model-based tests.
#' @param dma Differential methylation analysis. If datatype=="methylomics", can be DMP (default) vs. DMR. Else NULL.
#' @param remove.xy.probes Logical. Only activated when datatype=="methylomics". Remove X- and Y-linked CpG probes.
#' @param meth_type Type of array: 450K array, EPIC array or EPIC v2 array
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
                          celltype_ref = "pbmcref",
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

  if (datatype == "scRNA-seq") {
    pgx <- pgx.createSingleCellPGX(
      counts = counts,
      samples = samples,
      contrasts = contrasts,
      organism = organism,
      celltype_ref = celltype_ref,
      batch = NULL,
      sc_compute_settings = sc_compute_settings
    )
    return(pgx)
  }

  ## -------------------------------------------------------------------
  ## Establish the source and analysis sample axes
  ## -------------------------------------------------------------------
  samples <- as.data.frame(samples, drop = FALSE)
  counts <- as.matrix(counts)
  if (is.null(contrasts)) contrasts <- samples[, 0]

  ## EPIC v2 ids carry a design suffix (cg00000029_TC21), and ~5,000 CpGs
  ## have several replicate probes. Keep the recommended replicate per CpG,
  ## renamed to the bare id, so the pgx only ever holds bare ids - the ids the
  ## annotation, clocks, cell references, EWAS catalog and masks key on.
  if (identical(datatype, "methylomics") &&
    any(grepl(EPICV2_SUFFIX, rownames(counts)))) {
    source.ids <- rownames(counts)
    counts <- .pgx_collapse_epicv2(counts)
    if (!is.null(X)) X <- .pgx_collapse_epicv2(X)
    if (!is.null(annot_table)) {
      ## The rows the collapse kept, in its (source) order.
      kept <- source.ids[!source.ids %in% attr(counts, "epicv2_dropped")]
      annot_table <- annot_table[match(kept, rownames(annot_table)), , drop = FALSE]
      rownames(annot_table) <- rownames(counts)
    }
  }
  contrasts <- contrasts.convertToLabelMatrix(contrasts, samples)
  contrasts <- fixContrastMatrix(contrasts)
  if (dotimeseries) {
    contrasts <- contrasts.addTimeInteraction(contrasts, samples)
  }
  contrasts[contrasts %in% c("", " ", "NA")] <- NA

  source.samples <- intersect(colnames(counts), rownames(samples))
  samples <- samples[source.samples, , drop = FALSE]
  samples <- utils::type.convert(samples, as.is = TRUE)
  rownames(samples) <- source.samples
  if (all(source.samples %in% rownames(contrasts))) {
    contrasts <- contrasts[source.samples, , drop = FALSE]
  }
  analysis.samples <- source.samples
  used.samples <- rownames(contrasts)[rowSums(!is.na(contrasts)) > 0]
  if (prune.samples) {
    analysis.samples <- intersect(analysis.samples, used.samples)
  }
  if (!length(analysis.samples)) {
    info("[createPGX] FATAL. no analysis samples")
    return(NULL)
  }

  ## -------------------------------------------------------------------
  ## Build the analysis matrix once
  ## -------------------------------------------------------------------
  preprocess.metadata <- NULL
  ran.preprocess <- is.null(X)
  if (ran.preprocess) {
    message("[pgx.createPGX] building X via playbase.preprocess::pgx.preprocess()")
    preprocess <- .pgx_preprocess_options(
      preprocess,
      counts,
      norm_method,
      average.duplicated,
      batch.correct.method
    )
    preprocess <- .pgx_resolve_batch_metadata(
      preprocess,
      counts,
      samples,
      contrasts,
      analysis.samples,
      batch.pars
    )

    analysis.index <- match(analysis.samples, colnames(counts))
    pp <- playbase.preprocess::pgx.preprocess(
      counts = counts[, analysis.index, drop = FALSE],
      samples = samples[analysis.samples, , drop = FALSE],
      contrasts = contrasts[analysis.samples, , drop = FALSE],
      annot = annot_table,
      options = preprocess
    )
    pp$alignment$cols <- analysis.index[pp$alignment$cols]
    settings$options <- pp$options
    preprocess.metadata <- .pgx_preprocess_metadata(pp)
    settings$preprocess <- preprocess.metadata
    X <- pp$X
    if (!is.null(annot_table)) annot_table <- pp$annot
  } else {
    if (!all(unlist(batch.correct.method) == "no_batch_correct")) {
      stop(
        "[pgx.createPGX] batch correction requires X = NULL",
        call. = FALSE
      )
    }
    prior <- 1
    preprocess.metadata <- .pgx_identity_preprocess_metadata(
      counts,
      X,
      space = "log2",
      prior = prior
    )
    ndup <- sum(duplicated(rownames(counts)))
    if (ndup > 0) {
      method <- if (average.duplicated) "average" else "unique"
      action <- if (average.duplicated) "Averaging" else "Making unique"
      message(
        "[pgx.createPGX] ",
        ndup,
        " duplicated feature(s) detected. ",
        action,
        " X while preserving source counts."
      )
      deduplicated <- playbase.preprocess::pp.deduplicate(
        X,
        method = method,
        space = "log2"
      )
      X <- deduplicated$X
      preprocess.metadata <- .pgx_group_preprocess_rows(
        preprocess.metadata,
        deduplicated$row_groups
      )
    }
    settings$preprocess <- preprocess.metadata
  }

  ## Keep analysis outputs on the final sample axis.
  keep.x <- which(colnames(X) %in% analysis.samples)
  X <- X[, keep.x, drop = FALSE]
  preprocess.metadata <- .pgx_subset_preprocess_cols(
    preprocess.metadata,
    keep.x
  )
  settings$preprocess <- preprocess.metadata

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

  playbase.preprocess::pp.alignCounts(
    counts,
    preprocess.metadata$alignment,
    X = X
  )

  if (datatype == "multi-omics") {
    has.colons <- mean(grepl("[:]", rownames(counts)), na.rm = TRUE) > 0.9
    if (!has.colons) stop("[pgx.createPGX] FATAL: features must have multi-omics prefix\n")
  }

  X <- as.matrix(X)

  ## -------------------------------------------------------------------
  ## Auto-scaling (scale down huge values, often in proteomics)
  ## -------------------------------------------------------------------
  counts_multiplier <- Inf

  ## -------------------------------------------------------------------
  ## Check the final analysis matrix
  ## -------------------------------------------------------------------
  message("[createPGX] conforming matrices...")
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
    xrows <- .pgx_first_source_rows(preprocess.metadata)
    newnames <- annotate_phospho_residue(rownames(counts))
    newnames <- make_unique(newnames)
    rownames(X) <- newnames[xrows]
    if (!is.null(annot_table)) {
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

  ## Appending the symbol rewrites cg00000029 -> cg00000029_RBL2, and every
  ## methylation lookup keys on the bare probe ID: the epigenetic clocks come
  ## back with zero coverage and no ages at all. Off for methylomics whatever
  ## the caller asked for, and corrected here so settings record what happened.
  if (identical(datatype, "methylomics") && isTRUE(convert.hugo)) {
    message("[pgx.createPGX] methylomics: not appending symbols to probe IDs")
    convert.hugo <- FALSE
  }

  ## add to setting info
  settings$filter.genes <- filter.genes
  settings$exclude.genes <- exclude.genes
  settings$only.known <- only.known
  settings$only.proteincoding <- only.proteincoding
  settings$convert.hugo <- convert.hugo
  settings$custom.geneset <- !is.null(custom.geneset)
  settings$batch.correct.method <- batch.correct.method

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
    xrows <- .pgx_first_source_rows(preprocess.metadata)
    reordered <- vapply(seq_len(nrow(counts)), function(i) {
      reorder_uniprots(rownames(counts)[i], feature.lengths[i])$feature
    }, character(1))
    rownames(X) <- reordered[xrows]
  }
  settings$preprocess <- preprocess.metadata

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
  pgx$genes <- playbase.annotation::pgx.annotate(
    organism = pgx$organism,
    probes = rownames(pgx$counts),
    datatype = pgx$datatype,
    meth_type = meth_type,
    probetype = pgx$probe_type,
    ortholog_species = pgx$ortholog_species,    
    annot_table = annot_table
  )$genes

  ## Array and genome build, so consumers never re-guess them from probe ids.
  ## 450K and EPIC annotate against hg19, EPIC v2 against hg38.
  if (identical(datatype, "methylomics")) {
    array <- if (isTRUE(meth_type %in% c("EPIC array", "EPIC v2 array"))) meth_type else "450K array"
    pgx$meth$array <- array
    pgx$meth$genome <- attr(pgx$genes, "genome") %||%
      (if (array == "EPIC v2 array") "hg38" else "hg19")
  }

  ## Reorder uniprots in pgx$genes. Valid for all datatypes.
  message("[pgx.createPGX] Reordering uniprot column in pgx$genes")
  hh <- grep("uniprot", tolower(colnames(pgx$genes)))
  if (length(hh) > 0) {
    feature.lengths <- NULL
    kk <- grep("length|size", tolower(colnames(pgx$genes)))
    if (length(kk) > 0) {
      feature.lengths <- as.character(pgx$genes[, kk[1]])
    }
    for (i in 1:nrow(pgx$genes)) {
      pgx$genes[i, hh[1]] <- reorder_uniprots(pgx$genes[i, hh[1]], feature.lengths[i])$feature
    }
  }

  if (is.null(pgx$genes)) {
    stop("[pgx.createPGX] FATAL: Could not build gene annotation")
  }

  if (!"symbol" %in% colnames(pgx$genes) && "gene_name" %in% colnames(pgx$genes)) {
    dbg("[pgx.createPGX] WARNING! no symbol column. copying deprecated gene_name column as symbol")
    pgx$genes$symbol <- pgx$genes$gene_name
  }

  if (all(is.na(pgx$genes$symbol))) {
    dbg("[pgx.createPGX] WARNING! all symbol NA. copying rownames as symbol")
    pgx$genes$symbol <- gsub(".*:|[.].*", "", rownames(pgx$genes))
  }

  ## -------------------------------------------------------------------
  ## Filter the analysis rows while retaining the source counts and alignment.
  if (filter.genes) {
    nexpr <- sum(rowSums(pgx$counts, na.rm = TRUE) == 0)
    message("[pgx.createPGX] Filtering out ", nexpr, " not-expressed genes...")
    filtered <- playbase.preprocess::pp.filterFeatures(
      pgx$counts,
      method = "zero_counts"
    )
    keep.x <- .pgx_rows_for_source_rows(
      pgx$settings$preprocess,
      filtered$keep_rows
    )
    pgx$X <- pgx$X[keep.x, , drop = FALSE]
    pgx$settings$preprocess <- .pgx_subset_preprocess_rows(
      pgx$settings$preprocess,
      keep.x
    )
  }

  ## -------------------------------------------------------------------
  ## Filter genes
  ## -------------------------------------------------------------------
  do.filter <- (only.known || only.proteincoding || !is.null(exclude.genes))
  if (do.filter) {
    source.keep <- seq_len(nrow(pgx$genes))
    if (only.known) {
      message("[pgx.createPGX] Removing genes without symbol...")
      symbol <- pgx$genes$symbol[source.keep]
      no.symbol <- is.na(symbol) | symbol %in% c("", "-")
      source.keep <- source.keep[which(!no.symbol)]
    }

    if (only.proteincoding) {
      message("[pgx.createPGX] Removing Rik/ORF/LOC genes...")
      symbol <- pgx$genes$symbol[source.keep]
      is.unknown <- grepl("^rik|^loc|^orf", tolower(symbol))
      is.unknown <- is.unknown & !is.na(symbol)
      source.keep <- source.keep[which(!is.unknown)]
    }

    if (!is.null(exclude.genes)) {
      message("[pgx.createPGX] excluding genes: ", exclude.genes)
      exstr <- strsplit(tolower(exclude.genes), split = "[ ,]")[[1]]
      exexpr <- paste(c(paste0("^", exstr), paste0(exstr, "$")), collapse = "|")
      exgene <- grepl(exexpr, tolower(pgx$genes$symbol[source.keep]))
      if (sum(exgene)) {
        source.keep <- source.keep[which(!exgene)]
      }
    }

    keep.x <- .pgx_rows_for_source_rows(
      pgx$settings$preprocess,
      source.keep
    )
    pgx$X <- pgx$X[keep.x, , drop = FALSE]
    pgx$settings$preprocess <- .pgx_subset_preprocess_rows(
      pgx$settings$preprocess,
      keep.x
    )
  }

  ## Methylomics arrays: if user-specified, remove X- & Y-linked CpG probes.
  if (pgx$datatype == "methylomics" & remove.xy.probes) {
    kk <- intersect(c("chr", "map"), colnames(pgx$genes))[1]
    if (length(kk) > 0) {
      jj <- grep("chrX|chrY|^X|^Y", pgx$genes[, kk], ignore.case = TRUE)
      if (length(jj) > 0) {
        message("[pgx.createPGX] Methylomics: removing ", length(jj), " X- & Y-linked CpG probes...")
        keep <- setdiff(seq_len(nrow(pgx$counts)), jj)
        keep.x <- .pgx_rows_for_source_rows(
          pgx$settings$preprocess,
          keep
        )
        pgx$X <- pgx$X[keep.x, , drop = FALSE]
        pgx$settings$preprocess <- .pgx_subset_preprocess_rows(
          pgx$settings$preprocess,
          keep.x
        )
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

    xrows <- .pgx_first_source_rows(pgx$settings$preprocess)
    rownames(pgx$genes) <- new.names
    pgx$genes$gene_name <- new.names ## gene_name should also be renamed??
    pgx$genes$feature <- new.names ## feature should also be renamed??
    ## Apply source annotation names to the aligned analysis rows.
    rownames(pgx$X) <- new.names[xrows]
  }

  ## -------------------------------------------------------------------
  ## Infer cell cycle/gender
  ## -------------------------------------------------------------------
  ## Both signature callers use the explicit count-scale policy helper.
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

  ## Build the model design in processed-X sample order.
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
    ## Methylomics clusters on M, not on beta. Beta is a bounded ratio and is
    ## heteroscedastic at both ends, so the top-SD feature selection inside
    ## pgx.clusterBigMatrix is dominated by probes sitting near 0.5 whatever
    ## their behaviour, while a probe moving 0.02 -> 0.10 - a fourfold change
    ## in methylation - barely registers. M is the convention for every
    ## distance-based method on arrays (Du 2010), and it is the scale the EWAS
    ## is fitted on, so the embedding and the model see the same data.
    ## Converted from the space preprocessing declared, never guessed from
    ## the range, so a pgx already processed as M is not transformed twice.
    Xclust <- NULL
    if (isTRUE(pgx$datatype == "methylomics")) {
      Xclust <- playbase.preprocess::pp.convertSpace(
        pgx$X,
        from = pgx$settings$preprocess$space %||% "beta",
        to = "mvalue"
      )
    }
    pgx <- pgx.clusterSamples(pgx, dims = c(2, 3), perplexity = NULL, X = Xclust, methods = mm)
  }

  ## pgx.initialize() lists tsne2d in obj.needed and returns NULL when it is
  ## absent, so a dataset computed with do.cluster = FALSE - which is how
  ## methylomics is computed now - passes pgx.checkObject(), reaches the
  ## Library, and then fails to open with "ERROR in object initialization".
  ## GMT already had this problem and solves it two hundred lines up by
  ## storing an explicitly empty matrix; tsne2d had no such fallback.
  ##
  ## NA coordinates rather than an empty matrix, and rownames kept, so that a
  ## consumer indexing by sample name gets NA instead of a subscript error.
  ## Every reader is a Dashboard board, and no datatype computed without
  ## clustering opens one.
  if (is.null(pgx$tsne2d)) {
    pgx$tsne2d <- matrix(
      NA_real_, nrow = ncol(pgx$X), ncol = 2,
      dimnames = list(colnames(pgx$X), c("tsne2d.1", "tsne2d.2"))
    )
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
    if (pgx$datatype == "scRNAseq") {
      mm <- c("pca", "tsne", "umap")
    }
    pgx <- pgx.clusterGenes(pgx, methods = mm, level = "gene")
  }

  if (max.genes > 0 && nrow(pgx$X) > max.genes) {
    filtered <- playbase.preprocess::pp.filterFeatures(
      pgx$X,
      method = "low_variance",
      n = as.integer(max.genes),
      rank_on = .pgx_count_scale_matrix(pgx)
    )
    dropped <- setdiff(seq_len(nrow(pgx$X)), filtered$keep_rows)
    pgx$filtered[["low.variance"]] <- paste(
      rownames(pgx$X)[dropped],
      collapse = ";"
    )
    pgx$X <- filtered$X
    pgx$settings$preprocess <- .pgx_subset_preprocess_rows(
      pgx$settings$preprocess,
      filtered$keep_rows
    )
  }

  pgx$timings <- c()
  GENETEST.METHODS <- c(
    "ttest", "ttest.welch", "ttest.rank",
    "voom.limma", "trend.limma", "notrend.limma",
    "edger.qlf", "edger.lrt", "deseq2.wald", "deseq2.lrt"
  )
  GENESETTEST.METHODS <- c(
    "fisher", "gsva", "ssgsea", "spearman",
    "camera", "fry", "fgsea", "globaltest"
  ) ## no GSEA, too slow...

  ## ------------------ gene level tests ---------------------
  if (!is.null(progress)) {
    progress$inc(0.1, detail = "testing genes")
  }

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
  if (!is.null(progress)) {
    progress$inc(0.2, detail = "testing gene sets")
  }

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
  if (!is.null(progress)) {
    progress$inc(0.3, detail = "extra modules")
  }
  message(
    "[pgx.computePGX] computing extra modules: ",
    paste0(extra.methods, collapse = "; ")
  )
  pgx <- compute_extra(
    pgx,
    extra = extra.methods,
    pgx.dir = pgx.dir,
    libx.dir = libx.dir,
    user_input_dir = user_input_dir
  )

  ## methylomics: ensure all OPG graphics & tables use beta.
  if (pgx$datatype == "methylomics") {
    pgx$X <- playbase.preprocess::pp.convertSpace(
      pgx$X,
      from = pgx$settings$preprocess$space,
      to = "beta",
      prior = pgx$settings$preprocess$prior,
      layers = pgx$settings$preprocess$layers
    )
    pgx$settings$preprocess$space <- "beta"

    ## Everything below reads beta, which is why it follows the conversion.
    require_epigenetics()

    ## Epigenetic clocks, fitted here rather than in the app: they are a pure
    ## function of the beta matrix, so once it exists there is nothing left to
    ## wait for, and this is the first line at which it does. Ten clocks take
    ## ~40s, which is a fine one-off at dataset creation and a bad thing to
    ## make every session pay. Fitted with no coverage floor and no clock
    ## selection - both are display choices the app applies to `cov`, not
    ## inputs that change any age.
    ##
    ## A failure here must not lose the dataset: the slot is left NULL and a
    ## warning is raised.
    message("[pgx.computePGX] fitting epigenetic clocks...")
    pgx$meth$clocks <- tryCatch(
      playbase.epigenetics::compute_clocks(pgx$X),
      error = function(e) {
        warning("[pgx.computePGX] epigenetic clocks failed: ", conditionMessage(e))
        NULL
      }
    )

    ## Cell composition, for the same reason and on the same terms: the
    ## proportions are a pure function of the beta matrix and a panel name, so
    ## nothing about them needs a session. All nine reference panels are fitted
    ## rather than just the blood default - the panel follows the tissue, not
    ## the user, and at ~56 KB on a 150-sample cohort storing the lot costs
    ## about a tenth of a percent of the pgx while making every panel switch in
    ## the app instant. ~93s at that size.
    message("[pgx.computePGX] estimating cell composition...")
    pgx$meth$cells <- tryCatch(
      playbase.epigenetics::compute_cell_counts(pgx$X),
      error = function(e) {
        warning("[pgx.computePGX] cell composition failed: ", conditionMessage(e))
        NULL
      }
    )

    ## Per-sample QC ledger (bimodality, DNAm age, predicted vs recorded sex,
    ## imprint drift, verdict), on the same terms as the clocks: a pure function
    ## of the beta matrix, stored once so the app never refits it. Reuses the
    ## stored Horvath age rather than fitting it a second time. The annotation
    ## is aligned to X because only its chromosome column is read, to decide
    ## whether sex can be predicted from the probes X still has.
    message("[pgx.computePGX] computing methylation sample QC...")
    pgx$meth$qc <- tryCatch(
      playbase.epigenetics::sample_qc(
        pgx$X,
        pgx$genes[match(rownames(pgx$X), rownames(pgx$genes)), , drop = FALSE],
        pgx$samples,
        ## methylclock names it Horvath, the wateRmelon fallback horvath.
        ages = pgx$meth$clocks$age$Horvath %||% pgx$meth$clocks$age$horvath
      ),
      error = function(e) {
        warning("[pgx.computePGX] methylation sample QC failed: ", conditionMessage(e))
        NULL
      }
    )
  }

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

## Forwards to playbase.annotation::pgx.addGMT().
pgx.add_GMT <- function(pgx,
                        custom.geneset = NULL,
                        max.genesets = 20000,
                        include_default_gmt = TRUE,
                        include_iea = TRUE, species_go = NULL) {
  playbase.annotation::pgx.addGMT(
    pgx = pgx,
    custom.geneset = custom.geneset,
    max.genesets = max.genesets,
    include_default_gmt = include_default_gmt,
    include_iea = include_iea,
    species_go = species_go
  )
}


## ----------------------------------------------------------------------
## -------------------------- end of file -------------------------------
## ----------------------------------------------------------------------
