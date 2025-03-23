##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

#' Convert GMT file to binary matrix
#'
#' This function converts a GMT file (Gene Matrix Transposed) to a binary matrix,
#' where rows represent genes and columns represent gene sets. The binary matrix
#' indicates the presence or absence of genes in each gene set.
#'
#' @param gmt A list representing the GMT file, where each element is a character
#'            vector representing a gene set.
#' @param max.genes The maximum number of genes to include in the binary matrix.
#'                  Defaults to -1, which includes all genes.
#' @param ntop The number of top genes to consider for each gene set. Defaults to -1,
#'             which includes all genes.
#' @param sparse A logical value indicating whether to create a sparse matrix
#'               (from the 'Matrix' package) or a dense matrix. Defaults to `TRUE`.
#' @param bg A character vector representing the background set of genes. Defaults to `NULL`,
#'           which considers all unique genes from the gene sets.
#' @param use.multicore A logical value indicating whether to use parallel processing
#'                     (via the 'parallel' package) for faster execution. Defaults to `TRUE`.
#'
#' @export
#'
#' @return A binary matrix representing the presence or absence of genes in each gene set.
#'         Rows correspond to genes, and columns correspond to gene sets.
#'
gmt2mat <- function(gmt, max.genes = -1, ntop = -1, sparse = TRUE,
                    bg = NULL, use.multicore = TRUE) {
  gmt <- gmt[order(-sapply(gmt, length))]
  gmt <- gmt[!duplicated(names(gmt))]
  if (ntop > 0) {
    gmt <- lapply(gmt, utils::head, n = ntop)
  }
  if (is.null(names(gmt))) names(gmt) <- paste("gmt.", 1:length(gmt), sep = "")
  if (is.null(bg)) {
    bg <- names(sort(table(unlist(gmt)), decreasing = TRUE))
  }
  if (max.genes < 0) max.genes <- length(bg)
  gg <- bg
  gg <- Matrix::head(bg, n = max.genes)
  gmt <- lapply(gmt, function(s) intersect(gg, s))
  kk <- unique(names(gmt))
  if (sparse) {
    D <- Matrix::Matrix(0, nrow = length(gg), ncol = length(kk), sparse = TRUE)
  } else {
    D <- matrix(0, nrow = length(gg), ncol = length(kk))
  }

  rownames(D) <- gg
  colnames(D) <- kk
  j <- 1
  if (use.multicore) {
    idx <- lapply(gmt, function(s) match(s, gg))
    idx[sapply(idx, length) == 0] <- 0
    idx <- sapply(1:length(idx), function(i) rbind(idx[[i]], i))
    idx <- matrix(unlist(idx[]), byrow = TRUE, ncol = 2)
    idx <- idx[!is.na(idx[, 1]), ]
    idx <- idx[idx[, 1] > 0, ]
    D[idx] <- 1
  } else {
    for (j in 1:ncol(D)) {
      k0 <- match(kk[j], names(gmt))
      ii0 <- which(gg %in% gmt[[k0]])
      if (length(ii0) > 0) D[ii0, j] <- +1
    }
  }
  D <- D[order(-Matrix::rowSums(D != 0, na.rm = TRUE)), ]
  D
}

#' Convert binary matrix to GMT list
#'
#' This function converts binary matrix to a GMT (Gene Matrix
#' Transposed) list, The binary matrix indicates the presence or
#' absence of genes in each gene set, where rows represent genes and
#' columns represent gene sets.
#'
#' @param mat A matrix with non-zero entries representing genes in
#'   each gene set where rows represent genes and columns represent
#'   gene sets.
#'
#' @export
#'
#' @return A list of vector representing each gene set. Each list
#'   element correspond to a gene set and is a vector of genes
#'
mat2gmt <- function(mat) {
  idx <- Matrix::which(mat != 0, arr.ind = TRUE)
  gmt <- tapply(rownames(idx), idx[, 2], list)
  names(gmt) <- colnames(mat)[as.integer(names(gmt))]
  gmt
}

#' Read data from a GMT file
#'
#' This function reads data from a GMT (Gene Matrix Transposed) file format.
#' The GMT format is commonly used to store gene sets or gene annotations.
#'
#' @param gmt.file The path to the GMT file.
#' @param dir (Optional) The directory where the GMT file is located.
#' @param add.source (Optional) Specifies whether to include the source information in the gene sets' names.
#' @param nrows (Optional) The number of rows to read from the GMT file.
#'
#' @export
#'
#' @return A list of gene sets, where each gene set is represented as a character vector of gene names.
#'
read.gmt <- function(gmt.file, dir = NULL, add.source = FALSE, nrows = -1) {
  f0 <- gmt.file
  if (strtrim(gmt.file, 1) == "/") dir <- NULL
  if (!is.null(dir)) f0 <- paste(sub("/$", "", dir), "/", gmt.file, sep = "")
  gmt <- utils::read.csv(f0, sep = "!", header = FALSE, comment.char = "#", nrows = nrows)[, 1]
  gmt <- as.character(gmt)
  gmt <- sapply(gmt, strsplit, split = "\t")
  names(gmt) <- NULL
  gmt.name <- sapply(gmt, "[", 1)
  gmt.source <- sapply(gmt, "[", 2)
  gmt.genes <- sapply(gmt, function(x) {
    if (length(x) < 3) {
      return("")
    }
    paste(x[3:length(x)], collapse = " ")
  })
  gset <- strsplit(gmt.genes, split = "[ \t]")
  gset <- lapply(gset, function(x) setdiff(x, c("", "NA", NA)))
  names(gset) <- gmt.name

  if (add.source) {
    names(gset) <- paste0(names(gset), " (", gmt.source, ")")
  }
  gset
}


#' Convert gene sets to GMT format
#'
#' This function converts a set of genes and their corresponding gene set names to the GMT (Gene Matrix Transposed) format.
#' The GMT format is commonly used to store gene sets or gene annotations.
#'
#' @param genes A character vector of gene names.
#' @param gs_name A character vector of gene set names corresponding to the genes.
#'
#' @export
#'
#' @return A list of gene sets in GMT format, where each gene set is represented as a character vector of gene names.
#'
convert.gmt <- function(genes, gs_name) {
  if (length(genes) != length(gs_name)) {
    stop("genes and gs_name must be the same length")
  }
  gmt <- tapply(genes, gs_name, list)
  gmt <- lapply(gmt, as.character)
  return(gmt)
}

#' Write gene sets to GMT file
#'
#' This function writes gene sets in GMT (Gene Matrix Transposed) format to a file.
#' The GMT format is commonly used to store gene sets or gene annotations.
#'
#' @param gmt A list of gene sets in GMT format, where each gene set is represented as a character vector of gene names.
#' @param file The file path to write the GMT file.
#' @param source A character vector specifying the source of each gene set. If not provided, the names of the gene sets are used as the source.
#'
#' @export
#'
#' @return NULL
#'
write.gmt <- function(gmt, file, source = NA) {
  gg <- lapply(gmt, paste, collapse = "\t")
  if (is.na(source)) source <- names(gmt)
  ee <- paste(names(gmt), "\t", source, "\t", gg, sep = "")
  write(ee, file = file)
}









#' Run Gene Set Enrichment Analysis (GSEA)
#'
#' This function performs Gene Set Enrichment Analysis (GSEA) on the given gene expression data.
#'
#' @param X The gene expression data matrix.
#' @param y The phenotype vector.
#' @param gmt The gene set collection in GMT format.
#' @param output.dir Directory where the output files will be saved. Default is NULL, which creates a temporary directory.
#' @param fdr False discovery rate (FDR) threshold for filtering results. Default is 0.25.
#' @param set.min Minimum size of gene sets to consider. Default is 15.
#' @param set.max Maximum size of gene sets to consider. Default is 500.
#' @param topgs Number of top gene sets to report. Default is 100.
#' @param nperm Number of permutations for statistical significance estimation. Default is 1000.
#' @param permute Permutation type, either "phenotype" or "gene_set". Default is "phenotype".
#' @param scoring The scoring scheme for ranking genes. Default is "weighted".
#' @param do.leading.edge Whether to compute the leading edge genes for each gene set. Default is FALSE.
#' @param rpt.label The label to use for the output files. Default is "gsea".
#' @param sort.y The sorting order of the phenotype classes. Default is "ref.last".
#' @param ref.type The reference class or classes in the phenotype vector. Default is c("0", "NT", "REF", "DMSO", "WT", "LOW").
#' @param make.sets Whether to create gene sets using the gene expression data. Default is TRUE.
#' @param clean.files Whether to remove old files in the output directory. Default is TRUE.
#' @param xmx Maximum memory allocation for GSEA in gigabytes. Default is 10.
#' @param force.permute Whether to force gene_set permutation when there are fewer than 7 samples. Default is FALSE.
#' @param metric The metric for ranking genes. Default is c("Signal2Noise", "Diff_of_Classes", "Pearson").
#' @param gsea.program Path to the GSEA program JAR file. Default is "/opt/GSEA/gsea-3.0.jar".
#' @param quiet Whether to suppress console output during the analysis. Default is FALSE.
#'
#' @return A data frame containing the GSEA results, including normalized enrichment score (NES), nominal p-value (NOM p-val), and false discovery rate (FDR) q-value.
#'
#' @export
run.GSEA <- function(X, y, gmt, output.dir = NULL, fdr = 0.25, set.min = 15,
                     set.max = 500, topgs = 100, nperm = 1000, permute = "phenotype",
                     scoring = "weighted", do.leading.edge = FALSE, rpt.label = "gsea",
                     sort.y = "ref.last", ref.type = c("0", "NT", "REF", "DMSO", "WT", "LOW"),
                     make.sets = TRUE, clean.files = TRUE, xmx = 10, force.permute = FALSE,
                     metric = c("Signal2Noise", "Diff_of_Classes", "Pearson"),
                     gsea.program = "/opt/GSEA/gsea-3.0.jar", quiet = FALSE) {
  if (!file.exists(gsea.program)) {
    stop("run.GSEA:: cannot find GSEA program ", gsea.program, ". Please install GSEA.\n")
  }

  ## prepare folders and filenames
  if (is.null(output.dir)) output.dir <- tempfile()
  if (output.dir %in% c(".", "./")) {
    stop("run.GSEA:: output.dir cannot be current directory\n")
  }
  output.dir <- gsub("[ -]", "_", output.dir)
  rpt.label <- gsub("[ -]", "_", rpt.label)
  if (is.null(output.dir)) {
    output.dir <- tempdir()
  } else {
    output.dir <- paste(sub("/$", "", output.dir), "/", sep = "")
  }
  output.dir <- sub("//$", "/", paste(output.dir, "/", sep = ""))
  if (!quiet) cat("output.dir=", output.dir, "\n")
  system(paste("mkdir -p", output.dir))
  if (length(dir(output.dir, pattern = "gsea")) > 0) {
    cat("warning:: removing old files in output folder\n")
    system(paste("rm -rf ", output.dir, "gsea_*", sep = ""))
  }

  ## omit NA samples or genesets without name??
  y <- as.vector(unlist(y))
  jj <- which(!is.na(y))
  X <- X[, jj]
  y <- y[jj]

  ## order samples??
  is.categorical <- (length(unique(y)) == 2)
  if (is.categorical) {
    y <- as.character(y)
    ref.type <- toupper(ref.type)
    ref.match <- paste(ref.type, collapse = "|")
    has.ref <- (length(intersect(toupper(y), ref.type)) > 0)
    ## has.ref <- grep(ref.match
    if (sort.y %in% c("ref.last") & has.ref) {
      ref <- intersect(toupper(y), ref.type)
      if (!quiet) cat("setting reference to", ref, " (last class in cls file)\n")
      jj <- order(1 * (toupper(y) %in% ref), colnames(X))
      X <- X[, jj]
      y <- y[jj]
    } else if (sort.y %in% c("no", "none", "no.sorting", FALSE)) {
      cat("no sorting\n")
    } else {
      jj <- order(y, decreasing = TRUE) ## reverse alfabetically.
      X <- X[, jj]
      y <- y[jj]
    }
  } else {
    ## assume continuous
  }

  ## SNR ranking metric (see http://software.broadinstitute.org/gsea/doc/)
  if (is.categorical) {
    y.classes <- unique(y) ## GSEA sorts like this
    j0 <- which(y == y.classes[2]) ## "reference"
    j1 <- which(y == y.classes[1])
    mu0 <- rowMeans(X[, j0, drop = FALSE], na.rm = TRUE)
    mu1 <- rowMeans(X[, j1, drop = FALSE], na.rm = TRUE)
    ## sd0 <- apply(X[, j0, drop = FALSE], 1, stats::sd, na.rm = TRUE)
    sd0 <- matrixStats::rowSds(X[, j0, drop = FALSE], na.rm = TRUE)
    ## sd1 <- apply(X[, j1, drop = FALSE], 1, stats::sd, na.rm = TRUE)
    sd1 <- matrixStats::rowSds(X[, j1, drop = FALSE], na.rm = TRUE)
    sd0 <- pmax(sd0, 0.2 * pmax(mu0, 1))
    sd1 <- pmax(sd1, 0.2 * pmax(mu1, 1))
    rnk <- (mu1 - mu0) / (sd0 + sd1)
  } else {
    rnk <- stats::cor(t(X), y, use = "pairwise")[, 1]
  }

  if (!(permute %in% c("phenotype", "gene_set"))) {
    stop("permute must be phenotype or gene_set")
  }
  if (min(table(y), na.rm = TRUE) < 7 && permute == "phenotype" && !force.permute) {
    if (!quiet) cat("WARNING: less than 7 samples. Switching to gene_set permutation\n")
    permute <- "gene_set"
  }
  if (!quiet) cat("permute.type = ", permute, "\n")

  if (is.categorical) {
    rpt.label2 <- paste(rpt.label, "_", paste(toupper(unique(y)), collapse = "_vs_"), sep = "")
  } else {
    rpt.label2 <- paste0(rpt.label, "_correlation")
  }

  if (!quiet) cat("rpt.label = ", rpt.label2, "\n")

  ## prepare gct file
  gct.file <- paste(output.dir, rpt.label2, ".gct", sep = "")
  cls.file <- paste(output.dir, rpt.label2, ".cls", sep = "")
  gmt.file <- paste(output.dir, rpt.label2, ".gmt", sep = "")
  rpt.file <- paste(output.dir, rpt.label2, ".gsea_report.txt", sep = "")
  gX <- data.frame(NAME = rownames(X), DESCRIPTION = NA, X)
  write("#1.2", file = gct.file)
  suppressWarnings(write(paste(nrow(X), ncol(X), sep = "\t"), file = gct.file, append = TRUE))
  suppressWarnings(utils::write.table(format(gX[, ], digits = 3),
    file = gct.file, append = TRUE,
    row.names = FALSE, quote = FALSE, sep = "\t"
  ))
  remove(gX)

  ## prepare class file
  if (is.categorical) {
    write(paste(ncol(X), "2 1", sep = " "), file = cls.file, append = FALSE)
    cls.types <- paste("#RED BLUE")
    cls.types <- paste("#", paste(unique(y), collapse = " "), sep = "")
    suppressWarnings(write(cls.types, file = cls.file, append = TRUE))
    suppressWarnings(write(paste(y, collapse = " "), file = cls.file, append = TRUE))
  } else {
    write("#numeric", file = cls.file, append = FALSE)
    suppressWarnings(write(paste0("#", rpt.label2), file = cls.file, append = TRUE))
    suppressWarnings(write(paste(y, collapse = " "), file = cls.file, append = TRUE))
  }

  ## prepare geneset file (cleanup names, GSEA does not like them...)
  gmt <- gmt[which(!is.na(names(gmt)))]
  gmt <- gmt[which(!(names(gmt) %in% c("", NA, "NA", " ")))]
  gmt <- gmt[order(-sapply(gmt, length))]
  gmt.origname <- names(gmt)
  altname <- toupper(gsub("[-+: *?!$,'|\\]", "_", names(gmt)))
  altname <- gsub("[)(]", "", altname)
  altname <- gsub("__|___|____", "_", altname)
  names(gmt.origname) <- names(gmt) <- altname
  if (sum(duplicated(altname) > 0)) {
    cat("warning:", sum(duplicated(altname)), "duplicated gene sets removed\n")
  }
  jj <- which(!duplicated(altname))
  gmt <- gmt[jj]
  altname <- altname[jj]
  gmt.origname <- gmt.origname[jj]
  gmt.sets <- sapply(altname, function(s) paste(s, "\t", s, "\t", paste(gmt[[s]], collapse = "\t"), sep = ""))
  write(gmt.sets, file = gmt.file)
  if (nperm < 0) {
    return(NULL)
  }
  if (!quiet) cat("analyzing", length(gmt.sets), "gene sets and", nrow(X), "genes\n")

  ## define proper metrics
  metric <- metric[1]
  if (!is.categorical) metric <- "Pearson"

  ## set command
  gsea.cmd <- paste("java -Xmx%XMXm -cp", gsea.program, "xtools.gsea.Gsea")
  gsea.par <- "-res %GCT -cls %CLS -gmx %GMX -mode Max_probe -norm meandiv -nperm %NPERM -permute %PERMUTE -scoring_scheme %SCORING -metric %METRIC -rpt_label %RPTLABEL -make_sets %MAKESETS -plot_top_x %TOPGS -rnd_seed timestamp -set_max %GSMAX -set_min %GSMIN -zip_report false -out %OUTDIR -gui false -collapse false -chip /opt/GSEA/chip_files/GENE_SYMBOL.chip"

  gsea.cmd <- sub("%XMX", xmx * 1024, gsea.cmd)
  gsea.par <- sub("%GCT", gct.file, gsea.par)
  gsea.par <- sub("%CLS", cls.file, gsea.par)
  gsea.par <- sub("%NPERM", nperm, gsea.par)
  gsea.par <- sub("%PERMUTE", permute, gsea.par)
  gsea.par <- sub("%MAKESETS", tolower(make.sets), gsea.par)
  gsea.par <- sub("%TOPGS", topgs, gsea.par)
  gsea.par <- sub("%GMX", gmt.file, gsea.par)
  gsea.par <- sub("%OUTDIR", output.dir, gsea.par)
  gsea.par <- sub("%RPTLABEL", rpt.label2, gsea.par)
  gsea.par <- sub("%GSMIN", set.min, gsea.par)
  gsea.par <- sub("%GSMAX", set.max, gsea.par)
  gsea.par <- sub("%SCORING", scoring, gsea.par)
  gsea.par <- sub("%METRIC", metric, gsea.par)
  gsea.cmd <- paste(gsea.cmd, gsea.par)
  gsea.cmd
  if (length(grep("%", gsea.cmd)) > 0) {
    stop("error: gsea.cmd")
  }

  ## run GSEA
  out.file <- paste(output.dir, "gsea_run.out", sep = "")
  if (!quiet) cat("executing GSEA...\n", gsea.cmd, "\n")
  system(paste(gsea.cmd, " > ", out.file, sep = ""))
  if (!quiet) cat("\t[OK]\n")

  ## gather results from GSEA output folder
  t1 <- paste(rpt.label2, ".Gsea[.]", sep = "")
  d1 <- dir(output.dir, pattern = t1, full.names = TRUE)
  d1
  d1 <- d1[length(d1)] ## take last report
  d1 <- sub("//", "/", d1)
  rr <- dir(d1, pattern = "^gsea_report_for_.*xls$")
  rr
  if (length(rr) == 0) {
    cat("Error: could not find report files in", d1, "\n")
    stop("exiting\n")
  }
  r1 <- utils::read.csv(paste(d1, rr[1], sep = "/"),
    sep = "\t",
    header = TRUE, row.names = 1, check.names = FALSE
  )
  r2 <- utils::read.csv(paste(d1, rr[2], sep = "/"),
    sep = "\t",
    header = TRUE, row.names = 1, check.names = FALSE
  )
  res <- rbind(r1, r2)
  res <- res[, which(colnames(res) != "")]
  dim(res)
  if (sum(!(rownames(res) %in% names(gmt.origname))) > 0) {
    cat("warning: could not map some altnames back\n")
  }

  ## clean up, extract any prefix DB-ID
  res$"GS DETAILS" <- NULL
  colnames(res) <- sub("GS<br> follow link to MSigDB", "GS", colnames(res))
  res$GS <- as.character(res$GS)
  ii <- match(rownames(res), names(gmt.origname))
  res$GS <- gmt.origname[ii]
  db <- as.character(rep("", nrow(res)))
  res <- cbind(DB = db, res)
  res$DB <- as.character(res$DB)
  jj <- grep(":", rownames(res))
  jj
  if (length(jj) > 0) {
    rjj <- rownames(res)[jj]
    res$DB[jj] <- sapply(strsplit(rjj, split = ":"), "[", 1)
  }

  ## filter on FDR, reorder
  res <- res[order(abs(res$NES), decreasing = TRUE), ]
  if (fdr < 1) {
    if ("FDR.q.val" %in% colnames(res)) jj <- which(res$FDR.q.val <= fdr)
    if ("FDR q-val" %in% colnames(res)) jj <- which(res[, "FDR q-val"] <= fdr)
    res <- res[jj, ]
  }

  ## annotate with leading edge genes
  if (nrow(res) > 0) {
    edge.genes <- rep("", nrow(res))
    i <- 1
    for (i in 1:nrow(res)) {
      nes <- res$NES[i]
      nes
      if (!is.na(nes)) {
        ri <- res[i, "RANK AT MAX"]
        sorted.rnk <- sort(rnk)
        if (nes < 0) topgenes <- Matrix::head(names(sorted.rnk), ri)
        if (nes > 0) topgenes <- Matrix::head(rev(names(sorted.rnk)), ri)
        k <- names(gmt.origname)[match(rownames(res)[i], gmt.origname)]
        gg <- intersect(topgenes, gmt[[k]])
        rnk0 <- rnk - stats::median(rnk)
        edge.genes[i] <- paste(gg[order(-abs(rnk0[gg]))], collapse = "|")
      }
    }
    res <- cbind(res, "LEADING GENES" = edge.genes)
  }

  ## write
  utils::write.table(res, file = rpt.file, quote = FALSE, sep = "\t", col.names = NA)

  ## LeadingEdge analysis
  if (do.leading.edge == TRUE) {
    gsea.LeadingEdgeAnalysis(output.dir)
  }

  ## cleanup
  if (clean.files) {
    unlink(gmt.file)
    unlink(gct.file)
    if (substring(output.dir, 1, 4) == "/tmp") {
      unlink(output.dir, recursive = TRUE, force = TRUE)
    }
  }
  return(res)
}


#' Run GSEA Preranked
#'
#' Run Gene Set Enrichment Analysis (GSEA) on preranked gene list.
#'
#' @param rnk A named numeric vector representing the preranked gene list.
#' @param gmt A named list where each element represents a gene set with the set name as the element name and a character vector of gene names as the element value.
#' @param output.dir The output directory where the GSEA results will be stored. If NULL, a temporary directory will be created.
#' @param fdr The false discovery rate (FDR) threshold for filtering enriched gene sets.
#' @param set.min The minimum size of gene sets to be considered in the analysis.
#' @param set.max The maximum size of gene sets to be considered in the analysis.
#' @param topgs The number of top gene sets to report.
#' @param nperm The number of permutations for computing the enrichment statistics.
#' @param rpt.label The label used for naming the GSEA report files.
#' @param make.sets Logical indicating whether to create gene set files.
#' @param clean.files Logical indicating whether to clean up intermediate files.
#' @param xmx The maximum memory allocation for the Java virtual machine (in gigabytes).
#' @param scoring The scoring scheme for ranking genes. Options are "weighted" or "classic".
#' @param chip The chip file specifying the gene symbol to probe mapping.
#' @param collapse Logical indicating whether to collapse probes to gene symbols.
#' @param do.leading.edge Logical indicating whether to perform leading edge analysis.
#' @param gsea.program The path to the GSEA program JAR file.
#' @param quiet Logical indicating whether to suppress console output.
#'
#' @return A data frame containing the GSEA results.
#'
#' @export
run.GSEA.preranked <- function(rnk, gmt, output.dir = NULL, fdr = 0.25,
                               set.min = 15, set.max = 500, topgs = 100, nperm = 1000,
                               rpt.label = "preranked", make.sets = TRUE,
                               clean.files = TRUE, xmx = 10, scoring = "weighted",
                               chip = "GENE_SYMBOL.chip", collapse = FALSE,
                               do.leading.edge = FALSE, gsea.program = "/opt/GSEA/gsea-3.0.jar",
                               quiet = FALSE) {
  if (!file.exists(gsea.program)) {
    stop("run.GSEA.preranked:: cannot find GSEA program ", gsea.program, "\n")
  }
  if (is.null(names(rnk))) {
    stop("run.GSEA.preranked:: rank vector must have names\n")
  }
  if (is.null(output.dir)) output.dir <- sub("file", "tmpdir", tempfile())
  if (output.dir %in% c(".", "./")) {
    stop("run.GSEA.preranked:: output.dir cannot be current directory\n")
  }
  output.dir <- sub("[- ]", "_", output.dir)
  output.dir <- paste(sub("/$", "", output.dir), "/", sep = "")
  output.dir <- sub("//$", "/", paste(output.dir, "/", sep = ""))
  if (!quiet) cat("setting output.dir=", output.dir, "\n")
  system(paste("mkdir -p", output.dir))
  if (length(dir(output.dir, pattern = "[.]GseaPreranked[.]")) > 0) {
    cat("warning:: removing old files in output folder\n")
    system(paste("rm -rf ", output.dir, "*GseaPreranked*", sep = ""))
    system(paste("rm -rf ", output.dir, "gsea_*", sep = ""))
    system(paste("rm -rf ", output.dir, "error_*", sep = ""))
  }
  if (substr(rpt.label, 1, 5) != "gsea_") rpt.label <- paste("gsea_", rpt.label, sep = "")
  gmt.file <- paste(output.dir, rpt.label, ".gmt", sep = "")
  rnk.file <- paste(output.dir, rpt.label, ".rnk", sep = "")
  limma.file <- paste(output.dir, rpt.label, "_limma.txt", sep = "")
  rpt.file <- paste(output.dir, rpt.label, ".gsea_report.txt", sep = "")
  rpt.label

  ## prepare rank file
  rnk <- sort(rnk) ## sorted!!
  utils::write.table(cbind(names(rnk), rnk),
    file = rnk.file, quote = FALSE,
    row.names = FALSE, col.names = FALSE, sep = "\t"
  )

  ## prepare geneset file (cleanup names, GSEA does not like them...)
  gmt0 <- gmt
  gmt <- gmt[which(!(names(gmt) %in% c("", NA, "NA", " ")))]
  gmt <- lapply(gmt, intersect, names(rnk)) ## remove genes out of ranked list
  gmt <- gmt[order(-sapply(gmt, length))]
  gmt.origname <- names(gmt)
  altname <- toupper(gsub("[-+: *?!$,'|\\]", "_", names(gmt)))
  altname <- gsub("__|___|____", "_", altname)
  names(gmt.origname) <- names(gmt) <- altname

  ## prepare geneset file (take out duplicated gene sets...)
  if (sum(duplicated(altname) > 0)) {
    if (!quiet) cat("warning:", sum(duplicated(altname)), "duplicated gene sets removed\n")
  }
  jj <- which(!duplicated(altname))
  gmt <- gmt[jj]
  altname <- altname[jj]
  gmt.origname <- gmt.origname[jj]
  gmt.sets <- sapply(altname, function(s) {
    paste(s, "\tNA\t", paste(gmt[[s]], collapse = "\t"), sep = "")
  })
  write(gmt.sets, file = gmt.file)
  if (nperm < 0) {
    return(NULL)
  }
  if (!quiet) cat("analyzing", length(gmt.sets), "gene sets and", length(rnk), "genes\n")

  chip <- paste0("/opt/GSEA/chip_files/", chip)

  ## run preranked
  gsea.cmd <- paste("java -Xmx%XMXm -cp", gsea.program, "xtools.gsea.GseaPreranked")
  gsea.par <- "-gmx %GMX -mode Max_probe -norm meandiv -nperm %NPERM -permutation_type gene_set -rnk %RNKFILE -scoring_scheme %SCORING -rpt_label %RPTLABEL -make_sets %MAKESETS -plot_top_x %TOPGS -rnd_seed timestamp -set_max %GSMAX -set_min %GSMIN -zip_report false -out %OUTDIR -gui false -collapse %COLLAPSE -chip %CHIP"

  gsea.cmd <- sub("%XMX", xmx * 1024, gsea.cmd)
  gsea.par <- sub("%NPERM", nperm, gsea.par)
  gsea.par <- sub("%MAKESETS", tolower(make.sets), gsea.par)
  gsea.par <- sub("%TOPGS", topgs, gsea.par)
  gsea.par <- sub("%GMX", gmt.file, gsea.par)
  gsea.par <- sub("%RNKFILE", rnk.file, gsea.par)
  gsea.par <- sub("%OUTDIR", output.dir, gsea.par)
  gsea.par <- sub("%RPTLABEL", rpt.label, gsea.par)
  gsea.par <- sub("%GSMIN", set.min, gsea.par)
  gsea.par <- sub("%GSMAX", set.max, gsea.par)
  gsea.par <- sub("%SCORING", scoring, gsea.par)
  gsea.par <- sub("%CHIP", chip, gsea.par)
  gsea.par <- sub("%COLLAPSE", collapse, gsea.par)
  gsea.cmd <- paste(gsea.cmd, gsea.par)

  ## run GSEA
  out.file <- paste(output.dir, "gsea_run.out", sep = "")
  if (!quiet) cat("executing GSEA...\n", gsea.cmd, "\n")
  system(paste(gsea.cmd, " > ", out.file, sep = ""))
  if (!quiet) cat("\t[OK]\n")

  ## gather results from GSEA output folder
  d1 <- dir(output.dir, pattern = "[.]GseaPreranked[.]", full.names = TRUE)
  d1 <- d1[length(d1)]
  d1 <- sub("//", "/", d1)
  d1
  rr <- dir(d1, pattern = "^gsea_report_for_.*xls$")
  rr
  if (length(rr) == 0) {
    cat("Error: could not find report files. exiting.\n")
    return(1)
  }
  r1 <- utils::read.csv(paste(d1, rr[1], sep = "/"),
    sep = "\t",
    header = TRUE, row.names = 1, check.names = FALSE
  )
  r2 <- utils::read.csv(paste(d1, rr[2], sep = "/"),
    sep = "\t",
    header = TRUE, row.names = 1, check.names = FALSE
  )
  res <- rbind(r1, r2)
  res <- res[, which(colnames(res) != "")]
  dim(res)
  if (sum(!(rownames(res) %in% names(gmt.origname))) > 0) {
    cat("warning: could not map some altnames back\n")
  }

  ## fill in missing value with closest neighbour
  na.sum <- sum(is.na(res$NES))
  na.sum
  if (na.sum > 0) {
    cat("WARNING:: GSEA returned some missing values!\n")
    jj <- which(is.na(res$NES))
    jj0 <- which(!is.na(res$NES))
    for (j in jj) {
      j1 <- jj0[utils::head(order(abs(res$ES[jj0] - res$ES[j])), 3)]
      res[j, "NES"] <- mean(res[j1, "NES"], na.rm = TRUE)
      res[j, "NOM p-val"] <- mean(res[j1, "NOM p-val"], na.rm = TRUE)
      res[j, "FDR q-val"] <- mean(res[j1, "FDR q-val"], na.rm = TRUE)
      res[j, "FWER p-val"] <- mean(res[j1, "FWER p-val"], na.rm = TRUE)
    }
  }

  ## clean up, extract any prefix DB-ID
  res$"GS DETAILS" <- NULL
  colnames(res) <- sub("GS<br> follow link to MSigDB", "GS", colnames(res))
  res$GS <- as.character(res$GS)
  ii <- match(rownames(res), names(gmt.origname))
  res$GS <- gmt.origname[ii]
  db <- as.character(rep("", nrow(res)))
  res <- cbind(DB = db, res)
  res$DB <- as.character(res$DB)
  jj <- grep(":", rownames(res))
  if (length(jj) > 0) {
    rjj <- rownames(res)[jj]
    res$DB[jj] <- sapply(strsplit(rjj, split = ":"), "[", 1)
  }

  ## filter on FDR, reorder
  res <- res[order(abs(res$NES), decreasing = TRUE), ]
  if (fdr < 1) {
    if ("FDR.q.val" %in% colnames(res)) jj <- which(res$FDR.q.val <= fdr)
    if ("FDR q-val" %in% colnames(res)) jj <- which(res[, "FDR q-val"] <= fdr)
    res <- res[jj, ]
  }

  ## annotate with leading edge genes
  if (nrow(res) > 0) {
    edge.genes <- rep("", nrow(res))
    i <- 1
    for (i in 1:nrow(res)) {
      nes <- res$NES[i]
      if (!is.na(nes)) {
        ri <- res$"RANK AT MAX"[i]
        if (nes < 0) topgenes <- names(rnk)[1:ri]
        if (nes > 0) topgenes <- rev(names(rnk))[1:ri]
        k <- names(gmt.origname)[match(rownames(res)[i], gmt.origname)]
        gg <- intersect(topgenes, gmt[[k]])
        edge.genes[i] <- paste(gg[order(-abs(rnk[gg]))], collapse = "|")
      }
    }
    res <- cbind(res, "LEADING GENES" = edge.genes)
  }

  ## write
  utils::write.table(res, file = rpt.file, quote = FALSE, sep = "\t", col.names = NA)

  ## LeadingEdge analysis
  if (do.leading.edge == TRUE) {
    gsea.LeadingEdgeAnalysis(output.dir, gsea.program = gsea.program)
  }

  ## cleanup big files
  if (clean.files) {
    unlink(gmt.file)
    if (substring(output.dir, 1, 4) == "/tmp") {
      unlink(output.dir, recursive = TRUE, force = TRUE)
    }
  }
  return(res)
}


#' Perform LeadingEdge analysis for GSEA results
#'
#' This function performs the LeadingEdge analysis for Gene Set Enrichment Analysis (GSEA) results.
#' It analyzes the leading edge genes of enriched gene sets to identify key genes and their
#' contribution to the enrichment signal.
#'
#' @param output.dir The output directory where GSEA results are stored.
#' @param ntop The number of leading edge genes to analyze (default: 100).
#' @param gsea.program The path to the GSEA program JAR file (default: "/opt/GSEA/gsea-3.0.jar").
#' @param xmx The maximum memory allocation for the GSEA program in gigabytes (default: 10).
#'
#' @return The function performs the leading edge analysis and does not return any value.
#'
#' @export
gsea.LeadingEdgeAnalysis <- function(output.dir, ntop = 100, gsea.program = "/opt/GSEA/gsea-3.0.jar", xmx = 10) {
  cat(">>> performing LeadingEdge analysis <<<\n")
  ## clean old results
  fe <- dir(output.dir, "LeadingEdgeTool",
    include.dirs = TRUE,
    full.names = TRUE
  )
  fe
  output.dir <- sub("/$", "", output.dir)
  output.dir
  if (length(fe) > 0) {
    cat("warning: removing old LEA results in", output.dir, "\n")
    rm.cmd <- paste("rm -fr ", output.dir, "/*LeadingEdgeTool*", sep = "")
    rm.cmd
    system(rm.cmd)
  }

  res.file <- dir(output.dir, "gsea_report.txt", full.names = TRUE)
  res.file
  res <- utils::read.csv(res.file, sep = "\t", row.names = 1)
  Matrix::head(res)[, 1:4]
  rpt.prefix <- gsub(".gsea_report.txt", "", res.file)
  rpt.prefix <- gsub(output.dir, "", rpt.prefix)
  rpt.prefix <- gsub("/", "", rpt.prefix)

  ## add leading-edge analysis
  f0 <- dir(output.dir, pattern = "[.]GseaPreranked[.][0.9]*")
  f1 <- dir(output.dir, pattern = "[.]Gsea[.][0.9]*")
  gsea.type <- "Gsea"
  if (length(f0) > 0) gsea.type <- "GseaPreranked"
  reg.expr <- paste(rpt.prefix, "[.]", gsea.type, "[.][0-9]*$", sep = "")
  rpt.dir <- dir(output.dir, pattern = reg.expr, full.names = TRUE)
  rpt.dir
  if (length(rpt.dir) > 1) {
    cat("warning:: multiple GSEA results in folder")
  }
  if (length(rpt.dir) == 0) {
    stop("warning:: could not find GSEA result path")
  }
  rpt.dir <- rpt.dir[1]
  rpt.label <- gsub("[.]GseaPreranked[.].*$", "", rpt.dir)
  rpt.label <- gsub("[.]Gsea[.].*$", "", rpt.label)
  rpt.label <- gsub(output.dir, "", rpt.label)
  rpt.label <- gsub("^/", "", rpt.label)
  rpt.label
  gmt.top <- Matrix::head(rownames(res), ntop)
  run.GSEA.LeadingEdge(
    rpt.path = rpt.dir, gmt = gmt.top,
    rpt.label = rpt.label, xmx = xmx,
    gsea.program = gsea.program,
    outdir = output.dir
  )
}


#' Run LeadingEdge analysis for GSEA results
#'
#' This function runs the LeadingEdge analysis for Gene Set Enrichment Analysis (GSEA) results.
#' It analyzes the leading edge genes of enriched gene sets to identify key genes and their
#' contribution to the enrichment signal.
#'
#' @param rpt.path The path to the GSEA result folder.
#' @param gmt A character vector specifying the gene sets to analyze.
#' @param xmx The maximum memory allocation for the GSEA program in gigabytes (default: 10).
#' @param rpt.label The label for the GSEA report (default: "gsea_leadingedge").
#' @param gsea.program The path to the GSEA program JAR file (default: "/opt/GSEA/gsea-3.0.jar").
#' @param outdir The output directory where the leading edge analysis results will be stored (default: ".").
#'
#' @return The function runs the leading edge analysis and does not return any value.
#'
#' @export
run.GSEA.LeadingEdge <- function(rpt.path, gmt, xmx = 10, rpt.label = "gsea_leadingedge",
                                 gsea.program = "/opt/GSEA/gsea-3.0.jar", outdir = ".") {
  gmt0 <- gmt
  gmt0 <- paste(gmt0, collapse = ",")
  gsea.cmd <- paste("java -Xmx%XMXm -cp ", gsea.program, " xtools.gsea.LeadingEdgeTool")
  gsea.par <- "-dir %PATH -gsets %GMT -out %OUT -rpt_label %RPT"
  gsea.cmd <- sub("%XMX", xmx * 1024, gsea.cmd)
  gsea.par <- sub("%PATH", rpt.path, gsea.par)
  gsea.par <- sub("%RPT", rpt.label, gsea.par)
  gsea.par <- sub("%OUT", outdir, gsea.par)
  gsea.par <- sub("%GMT", gmt0, gsea.par)
  cmd <- paste(gsea.cmd, gsea.par)
  cat(cmd, "\n")
  system(cmd)
}




#' Run GSEA analysis using pre-ranked or phenotype comparison mode
#'
#' This function performs Gene Set Enrichment Analysis (GSEA) using either pre-ranked mode or phenotype comparison mode.
#' In pre-ranked mode, the function performs GSEA on a pre-ranked gene list. In phenotype comparison mode,
#' the function compares different phenotypes to identify enriched gene sets.
#'
#' @param X The gene expression matrix or pre-ranked gene list.
#' @param gmt The GMT file or gene set collection.
#' @param Y The class labels or binary matrix indicating the phenotypes for phenotype comparison mode.
#' @param fdr The false discovery rate (FDR) threshold for selecting significant gene sets. Default is 1.
#' @param set.min The minimum gene set size. Default is 15.
#' @param set.max The maximum gene set size. Default is 500.
#' @param topgs The number of top gene sets to report. Default is 100.
#' @param nperm The number of permutations for estimating the statistical significance. Default is 1000.
#' @param rpt.label The label for the analysis report. Default is "analysis".
#' @param output.dir The output directory for saving the results. Default is NULL (results are not saved).
#' @param xmx The maximum memory allocation for the GSEA program in gigabytes (GB). Default is 10.
#' @param scoring The scoring scheme for ranking genes. Default is "weighted".
#' @param sort.y The sorting order of phenotypes for phenotype comparison mode. Default is "ref.last".
#' @param ref.type The reference type for phenotype comparison mode. Default is "0".
#' @param clean.files Logical value indicating whether to clean intermediate files. Default is TRUE.
#'
#' @return A list containing the results of the GSEA analysis, including q-values and NES (Normalized Enrichment Scores).
#'
#' @export
justGSEA <- function(X, gmt, Y = NULL, fdr = 1, set.min = 15, set.max = 500, topgs = 100, nperm = 1000,
                     rpt.label = "analysis", output.dir = NULL, xmx = 10, scoring = "weighted",
                     sort.y = "ref.last", ref.type = c("0", "nt", "ref", "dmso", "wt", "dsmo"),
                     clean.files = TRUE) {
  if (is.vector(X)) {
    nx <- names(X)
    X <- matrix(X, ncol = 1)
    rownames(X) <- nx
  }

  ## prepare matrices
  if (!is.null(output.dir)) {
    output.dir <- sub("/$", "", output.dir)
    output.dir <- sub("/$", "", output.dir)
  }

  kk <- c()
  if (is.null(Y)) {
    cat("preranked-GSEA for", NCOL(X), "phenotype(s)\n")
    nes <- matrix(NA, length(gmt), ncol(X))
    rownames(nes) <- names(gmt)
    colnames(nes) <- colnames(X)
    qval <- nes
    j <- 1
    for (j in 1:ncol(X)) {
      output.dir0 <- output.dir
      if (NCOL(X) > 1 && !is.null(output.dir)) {
        output.dir0 <- paste(output.dir, "/", colnames(X)[j], sep = "")
      }
      rnk <- X[, j]
      names(rnk) <- rownames(X)
      res <- matrix(NA, nrow = 0, ncol = 0)
      res <- run.GSEA.preranked(rnk,
        gmt = gmt,
        output.dir = output.dir0, rpt.label = rpt.label,
        fdr = fdr, set.min = set.min, set.max = set.max,
        topgs = topgs, nperm = nperm, scoring = scoring,
        clean.files = clean.files, xmx = xmx,
        do.leading.edge = FALSE, make.sets = FALSE
      )
      if (nrow(res) > 0) {
        ii <- match(res$GS, rownames(nes))
        nes[ii, j] <- res[, "NES"]
        qval[ii, j] <- res[, "FDR q-val"]
        kk <- unique(c(kk, res$GS))
      }
    }
  } else {
    if (is.vector(Y)) {
      ny <- names(Y)
      Y <- matrix(Y, ncol = 1)
      rownames(Y) <- ny
    }
    if (sum(duplicated(colnames(Y))) > 0) {
      cat("warning:: duplicated column names. results may be overwritten.")
    }
    cat("GSEA for", ncol(Y), "phenotype comparison(s)\n")
    nes <- matrix(NA, length(gmt), ncol(Y))
    rownames(nes) <- names(gmt)
    colnames(nes) <- colnames(Y)
    qval <- nes
    j <- 1
    for (j in 1:ncol(Y)) {
      output.dir0 <- output.dir
      rpt.label0 <- rpt.label
      if (NCOL(Y) > 1 && !is.null(output.dir)) {
        rpt.label0 <- colnames(Y)[j]
        output.dir0 <- paste(output.dir, "/", rpt.label0, sep = "")
      }
      output.dir0
      res <- matrix(NA, nrow = 0, ncol = 0)
      res <- run.GSEA(X, Y[, j],
        gmt = gmt,
        output.dir = output.dir0, rpt.label = rpt.label0,
        fdr = fdr, set.min = set.min, set.max = set.max,
        topgs = topgs, nperm = nperm, scoring = scoring,
        clean.files = clean.files, xmx = xmx,
        sort.y = sort.y, ref.type = ref.type,
        do.leading.edge = FALSE, make.sets = FALSE
      )
      if (nrow(res) > 0) {
        ii <- match(res$GS, rownames(nes))
        nes[ii, j] <- res[, "NES"]
        qval[ii, j] <- res[, "FDR q-val"]
        kk <- unique(c(kk, res$GS))
      }
    }
  }
  kk <- unique(kk)
  jj <- which(rownames(nes) %in% kk)
  nes <- nes[jj, , drop = FALSE]
  qval <- qval[jj, , drop = FALSE]
  if (!is.null(output.dir)) {
    system(paste("mkdir -p", output.dir))
    utils::write.csv(nes, file = paste(output.dir, "/justGSEA_", rpt.label, "-NES.csv", sep = ""))
    utils::write.csv(qval, file = paste(output.dir, "/justGSEA_", rpt.label, "-qvalue.csv", sep = ""))
  }
  res <- c()
  res$qvalue <- qval
  res$NES <- nes
  return(res)
}



## ========================================================================
## ======================== ACCESS FUNCTIONS ==============================
## ========================================================================


#' Retrieve GSEA Output
#'
#' This function retrieves the output files and data generated from a GSEA analysis.
#'
#' @param path The path to the GSEA output directory.
#' @param i The index of the GSEA output folder to retrieve. Default is 1.
#' @param raster.png Logical value indicating whether to convert PNG images to raster format using grid::rasterGrob. Default is FALSE.
#'
#' @return A list containing the retrieved GSEA output files and data.
#'
#' @export
getGseaOutput <- function(path = "../analysis_v1b/output_GSEA/Th17_mut_2h_VS_mut_ut", i = 1, raster.png = FALSE) {
  ## untangle Gsea subfolder
  gsea_dir <- dir(path, full.names = TRUE)[grep("\\.Gsea\\.", dir(path))]
  gsea_dir
  if (length(gsea_dir) == 0) {
    cat("No GSEA output folders\n")
    return(NULL)
  }
  LE_dir <- dir(path, full.names = TRUE)[grep("\\.LeadingEdgeTool\\.", dir(path))]
  gsea_dir <- gsea_dir[1]
  LE_dir <- LE_dir[1]

  topgs <- sub(".xls", "", dir(gsea_dir, pattern = "[A-Z].*xls", full.names = FALSE))
  topgs

  M <- utils::read.csv(file.path(LE_dir, "leading_edge_matrix_for_results.1.gct"),
    sep = "\t", header = TRUE, skip = 2, row.names = 1
  )
  M <- as.matrix(M[, -1])
  rownames(M) <- sub("_signal", "", rownames(M))

  output <- list()
  output$dir <- gsea_dir
  output$index <- file.path(gsea_dir, "index.html")
  reports_xls <- dir(gsea_dir, pattern = "^gsea_report.*xls", full.names = TRUE)
  report <- lapply(reports_xls, function(x) utils::read.csv(x, sep = "\t", check.names = FALSE))
  output$report <- do.call(rbind, report)
  rownames(output$report) <- report$NAME

  #
  ff <- dir(path, full.names = TRUE)
  report_name <- ff[grep("gsea_report.txt$", ff)]
  if (!is.null(report_name)) {
    output$reportx <- utils::read.csv(report_name, sep = "\t", check.names = FALSE)
  }

  output$enplots <- sapply(topgs, function(g) {
    dir(gsea_dir,
      pattern = paste0("^enplot_", g, "_[0-9]{1,4}.png"),
      full.names = TRUE
    )
  })
  output$snapshots_html <- dir(gsea_dir, pattern = "snapshot.html", full.names = TRUE)
  output$heatmaps <- sapply(topgs, function(g) {
    dir(gsea_dir, pattern = paste0("^", g, "_[0-9]{1,4}.png"), full.names = TRUE)
  })
  details_xls <- file.path(gsea_dir, paste0(topgs, ".xls"))
  output$details <- lapply(details_xls, function(f) utils::read.csv(f, sep = "\t", check.names = FALSE))
  output$topgs <- topgs
  names(output$enplots) <- topgs
  names(output$heatmaps) <- topgs
  names(output$details) <- topgs

  output$edb_dir <- file.path(gsea_dir, "edb")
  output$gmt_file <- dir(output$edb_dir, pattern = "*.gmt", full.names = TRUE)
  output$cls_file <- dir(output$edb_dir, pattern = "*.cls", full.names = TRUE)
  output$rnk_file <- dir(output$edb_dir, pattern = "*.rnk", full.names = TRUE)

  output$rnk <- utils::read.csv(output$rnk_file, sep = "\t", row.names = 1, header = FALSE)
  output$cls <- utils::read.table(output$cls_file, skip = 2)[1, ]

  dir(LE_dir)
  output$LE_index <- file.path(LE_dir, "index.html")
  output$LE_heatmap <- file.path(LE_dir, "leading_edge_heat_map_unclustered.png")
  output$LE_heatmap_clustered <- file.path(LE_dir, "leading_edge_heat_map_clustered.png")
  output$LE_matrix <- M

  if (raster.png) {
    cat("rastering PNG images...\n")


    rasterPNG <- function(p) grid::rasterGrob(grDevices::as.raster(png::readPNG(p)), interpolate = FALSE)
    output$heatmaps <- lapply(output$heatmaps, rasterPNG)
    output$enplots <- lapply(output$enplots, rasterPNG)
    output$LE_heatmap <- rasterPNG(output$LE_heatmap)
    output$LE_heatmap_clustered <- rasterPNG(output$LE_heatmap_clustered)
  }

  return(output)
}





## ========================================================================
## =========================== PLOTTING ===================================
## ========================================================================


#' @export
.bluered <- function(n = 64) {
  gplots::colorpanel(n, "royalblue3", "grey90", "indianred3")
}





#' GSEA Enrichment Plot
#'
#' Generates an enrichment plot for GSEA analysis.
#'
#' @param rnk A numeric vector representing the ranked list.
#' @param gset A character vector specifying the gene set of interest.
#' @param names A character vector of names for the ranked list.
#' @param main The main title of the plot.
#' @param decreasing Logical indicating whether the ranked list should be sorted in decreasing order.
#' @param cex The size of the text labels.
#' @param cex.main The size of the main title text.
#' @param len.main The maximum number of characters for the main title before line breaks are applied.
#' @param lab.line The position of the axis labels.
#' @param cex.lab The size of the axis labels.
#' @param main.line The position of the main title.
#' @param xlab The label for the x-axis.
#' @param res The resolution of the plot.
#' @param ylab The label for the y-axis.
#' @export
#' @return plot
gsea.enplot <- function(rnk, gset, names = NULL, main = NULL,
                        decreasing = TRUE, cex = 1, cex.main = 0.95, len.main = 40,
                        lab.line = c(0.8, 2), cex.lab = 0.85, main.line = 0.3,
                        xlab = "Rank in ordered dataset", res = 1200,
                        ylab = "Rank metric") {
  if (!is.null(names)) names(rnk) <- names
  rnk <- rnk[!is.na(rnk)]
  rnk <- rnk[order(rnk + 1e-8 * stats::rnorm(length(rnk)), decreasing = decreasing)]

  ## ranked list metric
  ii <- (1:length(rnk))
  if (length(ii) > res) ii <- ii[seq(1, length(ii), length(ii) / res)]
  qq <- stats::quantile(rnk[ii], probs = c(0.01, 0.99), na.rm = TRUE)
  y1 <- qq[2]
  y0 <- qq[1]
  dy <- 0.25 * (y1 - y0)
  plot(ii, rnk[ii],
    type = "h", col = "grey", ylim = c(y0 - dy, y1),
    xlab = NA, ylab = NA, xaxt = "n"
  )
  graphics::mtext(xlab, 1, line = lab.line[1], cex = cex.lab)
  graphics::mtext(ylab, 2, line = lab.line[2], cex = cex.lab)
  graphics::abline(h = 0, lty = 2, lwd = 0.5)

  ## gene set barcode
  jj <- match(gset, names(rnk))
  jj <- jj[!is.na(jj)]
  w1 <- ifelse(length(jj) < 100, 0.7, 0.4)
  w1 <- ifelse(length(jj) < 50, 1, w1)
  col1 <- "grey10"
  if (all(grepl("[:]", gset))) {
    dtype <- sub(":.*", "", names(rnk))
    ntype <- length(unique(dtype))
    col1 <- rainbow(ntype)[factor(dtype[jj])]
  }
  graphics::arrows(jj, (y0 - dy), jj, y0, col = col1, lwd = w1 * cex, length = 0)

  ## red/blue bar at bottom
  kk <- c(seq(1, length(rnk) * 0.99, floor(length(rnk) / 20)), length(rnk))
  length(kk)
  i <- 1
  pal <- grDevices::colorRampPalette(c(omics_colors("brand_blue"), omics_colors("grey"), omics_colors("red")))(32)
  for (i in 1:(length(kk) - 1)) {
    r <- mean(rnk[kk[c(i, i + 1)]], na.rm = TRUE)
    r1 <- (r / max(abs(rnk), na.rm = TRUE))
    r1 <- abs(r1)**0.5 * sign(r1)
    irnk <- floor(31 * (1 + r1) / 2)
    cc <- pal[1 + irnk]
    graphics::rect(kk[i], y0 - 1.05 * dy, kk[i + 1], y0 - 0.65 * dy, col = cc, border = NA)
  }

  ## weighted cumulative random walk
  x0 <- 1 * (names(rnk) %in% gset)
  x0 <- x0 * abs(rnk)
  n0 <- sum(!(names(rnk) %in% gset))
  n1 <- sum(names(rnk) %in% gset)
  r0 <- cumsum(x0 == 0) / sum(x0 == 0)
  r1 <- cumsum(x0) / (1e-4 + sum(x0))
  rnk.trace <- (r1 - r0)
  rnk.trace <- rnk.trace / max(abs(rnk.trace), na.rm = TRUE) * 0.9
  if (max(rnk.trace, na.rm = TRUE) >= abs(min(rnk.trace, na.rm = TRUE))) rnk.trace <- rnk.trace * abs(y1)
  if (max(rnk.trace, na.rm = TRUE) < abs(min(rnk.trace, na.rm = TRUE))) rnk.trace <- rnk.trace * abs(y0)
  if (!decreasing) rnk.trace <- -1 * rnk.trace
  graphics::lines(ii, rnk.trace[ii], col = "green", type = "l", lwd = 2.4)

  if (is.null(main)) main <- "Enrichment plot"
  tt.main <- as.character(main)
  if (nchar(tt.main) > len.main) {
    tt.main < breakstring(tt.main, len.main) ## pgx-funtions.R
  }
  graphics::title(main = tt.main, cex.main = cex.main, line = main.line)
}




#' gsea.ftplot
#'
#' A function to generate a GSEA (Gene Set Enrichment Analysis) plot.
#'
#' @param x The input data matrix.
#' @param rft The metric to use for row ordering. If character, the metric function will be applied to each row of the matrix \code{x}.
#' @param cft The metric to use for column ordering. If character, the metric function will be applied to each column of the matrix \code{x}.
#' @param rmax The maximum number of rows to include in the plot.
#' @param cmax The maximum number of columns to include in the plot.
#' @param main The title of the plot.
#' @param cex.main The size of the main title.
#' @param xnames The names to use for the rows of the matrix.
#' @param ynames The names to use for the columns of the matrix.
#' @param cft.lab The label for the column metric.
#' @param rft.lab The label for the row metric.
#' @param cexRow The size of the row names.
#' @param cexCol The size of the column names.
#' @param scale The scaling method for the heatmap. Options are "rows" (scale rows) and "cols" (scale columns).
#' @param rsort The method for ordering the rows. Options are "hclust" (hierarchical clustering), "metric" (order by metric value), and "none" (no ordering).
#' @param csort The method for ordering the columns. Options are "hclust" (hierarchical clustering), "metric" (order by metric value), and "none" (no ordering).
#' @param p The power to raise the heatmap values to.
#' @param mar The margin sizes for the plot.
#' @param layout.reset Whether to reset the layout after plotting.
#' @param sort.decreasing Whether to sort the metric values in decreasing order.
#' @param col The color palette for the heatmap.
#' @param legend The legend labels for the heatmap colors.
#' @param cex.legend The size of the legend labels.
#' @param layout.widths The widths of the plot layout.
#' @param layout.heights The heights of the plot layout.
#'
#' @return None
#' @export
gsea.ftplot <- function(x, rft = "var", cft = "var",
                        rmax = nrow(x), cmax = ncol(x), main = "", cex.main = 1.5,
                        xnames = NULL, ynames = NULL,
                        cft.lab = "col metric", rft.lab = "row metric",
                        cexRow = 0.65, cexCol = 0.65, scale = "rows",
                        rsort = "hclust", csort = "hclust", p = 1,
                        mar = c(8, 4, 4, 2), layout.reset = TRUE,
                        sort.decreasing = FALSE, col = gplots::bluered(64),
                        legend = NULL, cex.legend = 1,
                        layout.widths = c(0.15, 0.15, 1),
                        layout.heights = c(0.15, 0.15, 1)) {
  if (class(rft) == "character") {
    rft.lab <- rft
    rft <- apply(x, 1, function(x) do.call(rft, list(x, na.rm = TRUE)))
  }
  if (class(cft) == "character") {
    cft.lab <- cft
    cft <- apply(x, 2, function(x) do.call(cft, list(x, na.rm = TRUE)))
  }
  if (length(rft) != nrow(x) || length(cft) != ncol(x)) {
    stop("dimension mismatch")
  }
  if (is.null(xnames)) xnames <- rownames(x)
  if (is.null(ynames)) ynames <- colnames(x)
  if (length(mar) == 2) mar <- c(mar[1], 4, 4, mar[2])

  ii <- 1:nrow(x)
  jj <- 1:ncol(x)
  if (nrow(x) > rmax) {
    ii <- ii[ii %in% Matrix::head(order(-abs(rft)), rmax)]
  }
  if (ncol(x) > cmax) {
    jj <- jj[jj %in% Matrix::head(order(-abs(cft)), rmax)]
  }
  x <- x[ii, jj]
  rft <- rft[ii]
  cft <- cft[jj]
  xnames <- xnames[ii]
  ynames <- ynames[jj]
  order.sign <- c(1, -1)[1 + (sort.decreasing == TRUE)]
  if (rsort == "hclust") {
    rc <- fastcluster::hclust(stats::as.dist(1 - stats::cor(t(x), use = "pairwise")))
    rc.order <- rc$order
  } else if (rsort == "metric") {
    rc.order <- order(-rft * order.sign)
  } else {
    rc.order <- 1:length(rft)
  }
  if (csort == "hclust") {
    cc <- fastcluster::hclust(stats::as.dist(1 - stats::cor(x, use = "pairwise")))
    cc.order <- cc$order
  } else if (csort == "metric") {
    cc.order <- order(cft * order.sign)
  } else {
    cc.order <- 1:length(cft)
  }

  ## ------------------------------------------------------------
  ## start plotting area
  ## ------------------------------------------------------------
  plotly::layout(matrix(c(6, 6, 1, 6, 6, 2, 5, 4, 3), 3, 3, byrow = TRUE),
    widths = layout.widths, heights = layout.heights
  )
  ## upper dendrogram (or skip)
  if (csort == "hclust") {
    graphics::par(mar = c(0.5, 0, mar[3], mar[4] + 0))
    plot(stats::as.dendrogram(cc), leaflab = "none", yaxt = "n")
  } else {
    graphics::par(mar = c(0.5, 0, mar[3], 0))
    graphics::frame()
  }
  graphics::title(main, cex.main = cex.main * 1.5, line = 0)

  ## column-feature barplot
  graphics::par(mar = c(1, 1.4, 0, mar[4] + 1.4))
  k0 <- "grey"
  graphics::barplot(cft[cc.order],
    width = 5 / 6, las = 3,
    xaxs = "i", add = FALSE, offset = 0, col = k0, cex.names = 1,
    names.arg = "", cex.axis = 1.0, border = FALSE,
    yaxt = "n", ylim = range(cft)
  )
  graphics::axis(4, las = 1)
  graphics::mtext(cft.lab, 4, padj = 4, cex = 0.8)

  ## heatmap
  graphics::par(mar = c(mar[1], 1.4, 1, mar[4] + 1.4))
  xx <- x[rc.order, cc.order]
  if (scale == "row") xx <- t(scale(t(xx), center = TRUE))
  if (scale == "col") xx <- scale(xx, center = TRUE)
  zlim1 <- c(-1, 1) * max(abs(xx**p), na.rm = TRUE)
  Matrix::image(1:ncol(xx), 1:nrow(xx), t(abs(xx)**p * sign(xx)),
    col = col, zlim = zlim1, xaxt = "n", xlab = "", ylab = "", yaxt = "n"
  )
  if (!is.null(xnames)) {
    graphics::mtext(xnames[rc.order], 4,
      at = (1:length(xnames)), adj = 0,
      padj = 0.5, las = 1, cex = cexRow, line = 0.5
    )
  }
  if (!is.null(ynames)) {
    graphics::mtext(ynames[cc.order], 1,
      at = (1:length(ynames)) - 0.5, adj = 1,
      padj = 1, las = 3, cex = cexCol, line = 0.5
    )
  }

  ## sample metrix barplot (vertical)
  graphics::par(mar = c(mar[1], 0, 1, 0))
  graphics::barplot(rft[rc.order],
    horiz = TRUE, las = 1, yaxs = "i", names.arg = "",
    border = FALSE, xlab = "", xlim = rev(range(rft))
  )
  graphics::mtext(rft.lab, 1, padj = 4, cex = 0.8)

  ## left dendrogram
  if (rsort == "hclust") {
    graphics::par(mar = c(mar[1], 2, 0, 0.5))
    plot(stats::as.dendrogram(rc),
      horiz = TRUE, leaflab = "none",
      yaxt = "n", yaxs = "i"
    )
  } else {
    graphics::par(mar = c(0, 0, 0, 0))
    graphics::frame()
  }

  ## legend
  if (!is.null(legend)) {
    graphics::par(mar = c(0, mar[2], mar[3], 0) + 0.5)
    graphics::frame()
    legend("topleft", legend = legend, cex = cex.legend)
  }

  ## reset
  graphics::par(mar = c(4, 4, 2, 2))
  if (layout.reset) {
    graphics::par(mfrow = c(1, 1))
    graphics::par(mar = c(4, 4, 2, 2))
  }
}


#' gsea.quick_report
#'
#' A function to generate a quick report of metaGSEA results in a specified folder.
#'
#' @param output.dir The directory containing the metaGSEA results.
#' @param pattern A pattern to match specific files in the directory (optional).
#'
#' @return A data frame with the summary of metaGSEA results.
#' @export
gsea.quick_report <- function(output.dir, pattern = NULL) {
  ## --------------------------------------------------------
  ## quick report of all metaGSEA results in folder
  ## --------------------------------------------------------
  output.dir <- sub("/$", "", output.dir)
  dd <- dir(output.dir, pattern = pattern, full.names = TRUE)
  dd
  gg <- lapply(dd, dir, pattern = "gsea_report.txt", full.names = TRUE)
  jj <- which(sapply(gg, length) > 0)
  dd <- dd[jj]
  gg <- unlist(gg[jj])
  tt <- gsub(paste(output.dir, "/", sep = ""), "", dd)
  tt <- gsub(output.dir, "", tt)
  tt <- gsub("^metagsva_|^metagsea_|^gsea_", "", tt)
  tt
  rr <- lapply(gg, utils::read.csv, sep = "\t", header = TRUE, row.names = 1, nrow = -2000)
  Matrix::head(rr[[1]])[, 1:7]
  n.set <- sapply(rr, function(x) nrow(x))
  np.pos <- sapply(rr, function(x) sum(x$NOM.p.val < 0.05 & x$NES > 0))
  np.neg <- sapply(rr, function(x) sum(x$NOM.p.val < 0.05 & x$NES < 0))
  nq.pos <- sapply(rr, function(x) sum(x$FDR.q.val < 0.25 & x$NES > 0))
  nq.neg <- sapply(rr, function(x) sum(x$FDR.q.val < 0.25 & x$NES < 0))
  nq.pos2 <- sapply(rr, function(x) sum(x$FDR.q.val < 0.01 & x$NES > 0))
  nq.neg2 <- sapply(rr, function(x) sum(x$FDR.q.val < 0.01 & x$NES < 0))
  R <- data.frame(
    comparison = tt, n.set,
    np.p005 = np.pos, nn.p005 = np.neg,
    np.q025 = nq.pos, nn.q025 = nq.neg,
    np.q001 = nq.pos2, nn.q001 = nq.neg2
  )
  R
}

#' gsea.radarplot
#'
#' A function to generate a radar plot for visualization of values.
#'
#' @param values A numeric vector of values.
#' @param names A character vector of names for the values (optional).
#' @param mar A numeric vector of margins (default: c(2, 5, 3, 5) * 2).
#' @param offset A numeric value specifying the offset (default: 0.1).
#' @param clust A logical value indicating whether clustering should be performed (default: TRUE).
#' @param main A character string specifying the plot title (default: "").
#' @param cex A numeric value specifying the character expansion factor (default: 0.8).
#' @param col A character vector of colors for the radar plot (default: c("steelblue2", "tomato")).
#' @return a plot
#' @export
gsea.radarplot <- function(values, names = NULL, mar = c(2, 5, 3, 5) * 2,
                           offset = 0.1, clust = TRUE, main = "",
                           cex = 0.8, col = c("steelblue2", "tomato")) {
  mx <- values
  if (is.null(names(mx)) & is.null(names)) {
    stop("value need names")
  }
  if (!is.null(names)) names(values) <- names
  M <- rbind(pmax(-mx, 0), pmax(mx, 0))
  M <- pmax(M / max(M, na.rm = TRUE), offset)
  M <- rbind(M, 1)
  colnames(M) <- names(values)
  rownames(M) <- c("neg", "pos", "unit")
  jj <- 1:ncol(M)
  if (clust) jj <- fastcluster::hclust(stats::dist(t(M)))$order
  M <- M[, jj]
  ii <- 1:3
  if (sum(mx < 0) == 0) ii <- c(2, 3)
  c0 <- c(col, NA)
  graphics::stars(M[ii, , drop = FALSE],
    radius = TRUE, cex = cex,
    locations = c(0, 0), key.loc = c(0, 0), lty = 1,
    col.stars = c0[ii], mar = mar, scale = FALSE,
    main = main, cex.main = 1
  )
}


#' gsea.enplotPDF
#'
#' A function to generate a PDF file containing ensemble plots for GSEA results.
#'
#' @param gsea.dir A character string specifying the directory path containing GSEA results.
#' @param gs A character vector of gene sets to include in the ensemble plots.
#' @param output.pdf A character string specifying the output PDF file path.
#' @return NULL
#' @export
gsea.enplotPDF <- function(gsea.dir, gs, output.pdf) {
  gsea.dir1 <- dir(gsea.dir,
    pattern = "gsea_preranked.GseaPreranked.",
    full.names = TRUE
  )[1]
  ep0 <- dir(gsea.dir1, pattern = "enplot_.*png", full.names = TRUE)
  ep1 <- dir(gsea.dir1, pattern = "enplot_.*png", full.names = FALSE)
  ep2 <- gsub("enplot_|_[0-9].*$", "", ep1)
  gs1 <- intersect(gs, ep2)
  jj <- match(gs1, ep2)
  ee <- paste(ep0[jj], collapse = " ")
  cmd <- paste("montage", ee, "-tile 6x4 -geometry 100% -border 5")
  cmd <- paste(cmd, output.pdf)
  cmd
  system(cmd)
}


## ========================================================================
## ========================================================================
## ========================================================================
