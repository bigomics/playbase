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
    gmt <- lapply(gmt, head, n = ntop)
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
  dim(D)
  rownames(D) <- gg
  colnames(D) <- kk
  j <- 1
  if (use.multicore) {
    idx <- parallel::mclapply(gmt, function(s) match(s, gg))
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
  D <- D[order(-Matrix::rowSums(D != 0)), ]
  D
}


#' Write data to a GCT file
#'
#' This function writes data to a GCT (Gene Cluster Text) file format.
#' The GCT format is commonly used to store gene expression data in a tabular format.
#'
#' @param X The data matrix or data frame to be written to the GCT file.
#'          Rows represent genes and columns represent samples.
#' @param file The file path to save the GCT file.
#'
#' @export
#'
#' @return None
#'
write.gct <- function(X, file) {
  gX <- data.frame(NAME = rownames(X), DESCRIPTION = NA, X)
  write("#1.2", file = file)
  suppressWarnings(write(paste(nrow(X), ncol(X), sep = "\t"), file = file, append = TRUE))
  suppressWarnings(write.table(format(gX[, ], digits = 3),
    file = file, append = TRUE,
    row.names = FALSE, quote = FALSE, sep = "\t"
  ))
}


#' Write data to a CLS file
#'
#' This function writes data to a CLS (CLustering Solutions) file format.
#' The CLS format is commonly used to store class labels or clustering assignments.
#'
#' @param y The vector of class labels or clustering assignments.
#' @param file The file path to save the CLS file.
#' @param name (Optional) The name of the CLS file.
#'
#' @export
#'
#' @return None
#'
write.cls <- function(y, file, name = "") {
  ## prepare class file
  is.numeric <- (length(setdiff(unique(y), NA)) > 3)
  if (is.numeric) {
    write("#numeric", file = file, append = FALSE)
    suppressWarnings(write(paste0("#", name), file = file, append = TRUE))
    suppressWarnings(write(paste(y, collapse = " "), file = file, append = TRUE))
  } else {
    write(paste(length(y), "2 1", sep = " "), file = file, append = FALSE)
    cls.types <- paste("#", paste(unique(y), collapse = " "), sep = "")
    suppressWarnings(write(cls.types, file = file, append = TRUE))
    suppressWarnings(write(paste(y, collapse = " "), file = file, append = TRUE))
  }
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
  gmt <- read.csv(f0, sep = "!", header = FALSE, comment.char = "#", nrows = nrows)[, 1]
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



#' Clean names for GSEA analysis
#'
#' This function cleans the names of a vector or list for use in GSEA (Gene Set Enrichment Analysis) analysis.
#' It replaces special characters and consecutive underscores in the names with a single underscore.
#'
#' @param s A vector or list with names to be cleaned.
#'
#' @export
#'
#' @return A character vector or list with cleaned names.
#'
gsea.clean_names <- function(s) {
  gsub("__|___", "_", gsub("[-+: *?!$,'|\\.]", "_", names(s)))
}


#' Perform GSEA analysis for all contrasts
#'
#' This function performs GSEA (Gene Set Enrichment Analysis) analysis for all the contrasts specified in the design matrix.
#' It iterates over the contrasts, runs GSEA for each contrast, and saves the results in the specified output directory.
#'
#' @param X The expression matrix.
#' @param gmt The gene set collection in GMT format.
#' @param design The design matrix.
#' @param contr.matrix The contrast matrix.
#' @param output_dir The directory to save the GSEA results.
#' @param set.min The minimum size of gene sets to consider in the analysis.
#' @param set.max The maximum size of gene sets to consider in the analysis.
#' @param fdr The false discovery rate (FDR) threshold.
#' @param skip.done Logical value indicating whether to skip the contrasts that have already been analyzed and have existing results in the output directory.
#'
#' @export
#'
#' @return A list of GSEA analysis results for each contrast.
#'
gsea.fitAllContrasts <- function(X, gmt, design, contr.matrix, output_dir,
                                 set.min = 15, set.max = 500, fdr = 0.25, skip.done = TRUE) {
  exp.matrix <- (design %*% contr.matrix)
  length(gmt)
  comparisons <- colnames(contr.matrix)
  outputs <- list()
  k <- 1
  for (k in 1:length(comparisons)) {
    cat(">>> running GSEA for comparison: ", comparisons[k], " <<<\n")

    ## set params

    comp.name <- comparisons[k]

    comp.name

    ## check if already done
    sub_dir <- file.path(output_dir, comp.name)
    sub_dir
    if (dir.exists(sub_dir) && skip.done) {
      cat("comparison ", comp.name, "already done. skipping.\n")
      ## should we read the output file??
      next
    }

    sel <- which(exp.matrix[, comp.name] != 0)
    xx <- X[, sel]
    yy <- 1 * (exp.matrix[sel, comp.name] > 0)

    j1 <- which(contr.matrix[, comp.name] > 0)
    j0 <- which(contr.matrix[, comp.name] < 0)
    cls1 <- paste(rownames(contr.matrix)[j1], sep = "and")
    cls0 <- paste(rownames(contr.matrix)[j0], sep = "and")
    c(cls0, cls1)
    yy <- c(cls0, cls1)[1 + yy]
    yy
    ref <- cls0
    dim(xx)
    yy
    ref

    ## run GSEA
    outputs[[k]] <- run.GSEA(xx, yy, gmt,
      fdr = fdr,
      do.leading.edge = TRUE, clean.files = FALSE,
      output.dir = sub_dir,
      set.min = set.min, set.max = set.max, topgs = 100,
      ref.type = ref, permute = "gene_set"
    )
  } ## for-loop comparisons

  return(outputs)
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
    sd0 <- apply(X[, j0, drop = FALSE], 1, sd, na.rm = TRUE)
    sd1 <- apply(X[, j1, drop = FALSE], 1, sd, na.rm = TRUE)
    sd0 <- pmax(sd0, 0.2 * pmax(mu0, 1))
    sd1 <- pmax(sd1, 0.2 * pmax(mu1, 1))
    rnk <- (mu1 - mu0) / (sd0 + sd1)
  } else {
    rnk <- stats::cor(t(X), y, use = "pairwise")[, 1]
  }

  if (!(permute %in% c("phenotype", "gene_set"))) {
    stop("permute must be phenotype or gene_set")
  }
  if (min(table(y)) < 7 && permute == "phenotype" && !force.permute) {
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
  suppressWarnings(write.table(format(gX[, ], digits = 3),
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
  r1 <- read.csv(paste(d1, rr[1], sep = "/"),
    sep = "\t",
    header = TRUE, row.names = 1, check.names = FALSE
  )
  r2 <- read.csv(paste(d1, rr[2], sep = "/"),
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
        rnk0 <- rnk - median(rnk)
        edge.genes[i] <- paste(gg[order(-abs(rnk0[gg]))], collapse = "|")
      }
    }
    res <- cbind(res, "LEADING GENES" = edge.genes)
  }

  ## write
  write.table(res, file = rpt.file, quote = FALSE, sep = "\t", col.names = NA)

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
  write.table(cbind(names(rnk), rnk),
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
  r1 <- read.csv(paste(d1, rr[1], sep = "/"),
    sep = "\t",
    header = TRUE, row.names = 1, check.names = FALSE
  )
  r2 <- read.csv(paste(d1, rr[2], sep = "/"),
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
      j1 <- jj0[head(order(abs(res$ES[jj0] - res$ES[j])), 3)]
      res[j, "NES"] <- mean(res[j1, "NES"])
      res[j, "NOM p-val"] <- mean(res[j1, "NOM p-val"])
      res[j, "FDR q-val"] <- mean(res[j1, "FDR q-val"])
      res[j, "FWER p-val"] <- mean(res[j1, "FWER p-val"])
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
  write.table(res, file = rpt.file, quote = FALSE, sep = "\t", col.names = NA)

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
  res <- read.csv(res.file, sep = "\t", row.names = 1)
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


#' Calculate Signal-to-Noise ratio for GSEA
#'
#' This function calculates the Signal-to-Noise ratio (SNR) for Gene Set Enrichment Analysis (GSEA).
#' It measures the difference in the means of two classes (0 and 1) relative to the standard deviation.
#' The SNR is used to rank genes and assess their relevance in the context of gene set enrichment.
#'
#' @param X The gene expression matrix.
#' @param Y The class labels or binary matrix indicating the two classes (0 and 1).
#'
#' @return A matrix of Signal-to-Noise ratios for each gene in the expression matrix.
#'
#' @export
gsea.snr <- function(X, Y) {
  if (NCOL(Y) == 1) Y <- matrix(Y, ncol = 1)
  if (ncol(X) != nrow(Y)) {
    stop("dimension mismatch")
  }
  S <- matrix(NA, nrow(X), ncol(Y))
  j <- 1
  for (j in 1:ncol(Y)) {
    yj <- Y[, j]
    mx <- cbind(
      rowMeans(X[, yj == 0, drop = FALSE], na.rm = TRUE),
      rowMeans(X[, yj == 1, drop = FALSE], na.rm = TRUE)
    )
    sx <- cbind(
      apply(X[, yj == 0, drop = FALSE], 1, sd, na.rm = TRUE),
      apply(X[, yj == 1, drop = FALSE], 1, sd, na.rm = TRUE)
    )
    sx <- pmax(pmax(sx, 0.2 * abs(mx)), 0.2)
    sx <- rowMeans(sx, na.rm = TRUE) * ncol(sx) ## robust "sum"
    S[, j] <- (mx[, 1] - mx[, 2]) / sx
  }
  colnames(S) <- colnames(Y)
  rownames(S) <- rownames(X)
  return(S)
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
    write.csv(nes, file = paste(output.dir, "/justGSEA_", rpt.label, "-NES.csv", sep = ""))
    write.csv(qval, file = paste(output.dir, "/justGSEA_", rpt.label, "-qvalue.csv", sep = ""))
  }
  res <- c()
  res$qvalue <- qval
  res$NES <- nes
  return(res)
}


#' Run GCEA (Gene Concept Enrichment Analysis) using GSEA
#'
#' This function performs Gene Concept Enrichment Analysis (GCEA) using Gene Set Enrichment Analysis (GSEA).
#' It combines gene set analysis and concept analysis to identify enriched gene sets and concepts.
#'
#' @param X The gene expression matrix or pre-ranked gene list.
#' @param sets The gene set collection for gene set analysis.
#' @param concepts The gene concept collection for concept analysis.
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
#' @param clean.files Logical value indicating whether to clean intermediate files. Default is TRUE.
#'
#' @return A list containing the results of the GCEA analysis, including gene set analysis results (gsea) and concept analysis results (gcea).
#'
#' @export
justGCEA <- function(X, sets, concepts, Y = NULL, fdr = 1, set.min = 15, set.max = 500,
                     topgs = 100, nperm = 1000, rpt.label = "analysis", output.dir = NULL,
                     xmx = 10, scoring = "weighted", clean.files = TRUE) {
  output.dir1 <- output.dir2 <- NULL
  if (!is.null(output.dir)) {
    output.dir1 <- paste(output.dir, "/gsea_out", sep = "")
    output.dir2 <- paste(output.dir, "/gcea_out", sep = "")
  }
  cat(">>> GCEA gene concept enrichment analysis <<<\n")
  cat("running gene set analysis\n")
  res1 <- justGSEA(X,
    gmt = sets, Y = Y, output.dir = output.dir1,
    topgs = topgs, nperm = nperm, rpt.label = rpt.label,
    xmx = xmx, scoring = scoring, fdr = fdr,
    set.min = set.min, set.max = set.max
  )
  cat("running concept analysis\n")
  res2 <- justGSEA(res1$NES,
    gmt = concepts, output.dir = output.dir2,
    topgs = topgs, nperm = nperm, rpt.label = rpt.label,
    xmx = xmx, scoring = scoring, fdr = fdr,
    set.min = set.min, set.max = set.max
  )
  res <- c()
  res$gsea <- res1
  res$gcea <- res2
  res
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

  M <- read.csv(file.path(LE_dir, "leading_edge_matrix_for_results.1.gct"),
    sep = "\t", header = TRUE, skip = 2, row.names = 1
  )
  M <- as.matrix(M[, -1])
  rownames(M) <- sub("_signal", "", rownames(M))

  output <- list()
  output$dir <- gsea_dir
  output$index <- file.path(gsea_dir, "index.html")
  reports_xls <- dir(gsea_dir, pattern = "^gsea_report.*xls", full.names = TRUE)
  report <- lapply(reports_xls, function(x) read.csv(x, sep = "\t", check.names = FALSE))
  output$report <- do.call(rbind, report)
  rownames(output$report) <- report$NAME


  ff <- dir(path, full.names = TRUE)
  report_name <- ff[grep("gsea_report.txt$", ff)]
  if (!is.null(report_name)) {
    output$reportx <- read.csv(report_name, sep = "\t", check.names = FALSE)
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
  output$details <- lapply(details_xls, function(f) read.csv(f, sep = "\t", check.names = FALSE))
  output$topgs <- topgs
  names(output$enplots) <- topgs
  names(output$heatmaps) <- topgs
  names(output$details) <- topgs

  output$edb_dir <- file.path(gsea_dir, "edb")
  output$gmt_file <- dir(output$edb_dir, pattern = "*.gmt", full.names = TRUE)
  output$cls_file <- dir(output$edb_dir, pattern = "*.cls", full.names = TRUE)
  output$rnk_file <- dir(output$edb_dir, pattern = "*.rnk", full.names = TRUE)

  output$rnk <- read.csv(output$rnk_file, sep = "\t", row.names = 1, header = FALSE)
  output$cls <- read.table(output$cls_file, skip = 2)[1, ]

  dir(LE_dir)
  output$LE_index <- file.path(LE_dir, "index.html")
  output$LE_heatmap <- file.path(LE_dir, "leading_edge_heat_map_unclustered.png")
  output$LE_heatmap_clustered <- file.path(LE_dir, "leading_edge_heat_map_clustered.png")
  output$LE_matrix <- M

  if (raster.png) {
    cat("rastering PNG images...\n")


    rasterPNG <- function(p) grid::rasterGrob(as.raster(png::readPNG(p)), interpolate = FALSE)
    output$heatmaps <- lapply(output$heatmaps, rasterPNG)
    output$enplots <- lapply(output$enplots, rasterPNG)
    output$LE_heatmap <- rasterPNG(output$LE_heatmap)
    output$LE_heatmap_clustered <- rasterPNG(output$LE_heatmap_clustered)
  }

  return(output)
}


#' Plot GSEA Enrichment Plots
#'
#' Plot the GSEA enrichment plots from the GSEA analysis results.
#'
#' @param gsea The GSEA analysis results obtained from \code{\link{getGseaOutput}}.
#' @param gsets Optional character vector specifying the gene sets to plot.
#'              If NULL (default), all available gene sets will be plotted.
#' @param ncol The number of columns in the grid arrangement of plots.
#'             Defaults to 5.
#'
#' @export
#'
#' @return None
gseaPlotEnplots <- function(gsea, gsets = NULL, ncol = 5) {
  enplots <- gsea$enplots
  kk <- 1:length(enplots)
  if (!is.null(gsets)) kk <- match(gsets, names(enplots))
  imgs <- lapply(enplots[kk], function(p) {
    grid::rasterGrob(as.raster(png::readPNG(p)),
      interpolate = TRUE
    )
  })
  grid.arrange(grobs = imgs, ncol = ncol)
}


#' Leading Edge Heatmap
#'
#' Generates a heatmap visualization of the leading edge analysis results.
#'
#' @param gsea An object containing GSEA results.
#' @param maxrow Maximum number of rows to display in the heatmap.
#' @param maxcol Maximum number of columns to display in the heatmap.
#' @param gsets A character vector specifying the gene sets to include in the analysis.
#' @param render The rendering method for the heatmap. Options are "gx.heatmap", "d3heatmap", or "heatmaply".
#' @param info.text Logical indicating whether to display additional information text.
#' @export
#' @return NULL if the leading edge matrix is too small (less than or equal to 1 row or 1 column),
#' otherwise, it returns the generated heatmap visualization.
gseaLeadingEdgeHeatmap <- function(gsea, maxrow = 60, maxcol = 60, gsets = NULL,
                                   render = "gx.heatmap", info.text = TRUE) {
  if (info.text) {
    cat("<br><h4>Leading Edge heatmap</h4>")
    cat("<br>Click here for the <a href='", gsea$LE_heatmap_clustered, "'>full Leading Edge heatmap</a>. ")
    cat("Complete Leading Edge analysis results are <a href='", gsea$LE_index, "'>here</a>")
  }
  LEmat <- gsea$LE_matrix
  if (!is.null(gsets)) {
    kk <- intersect(gsets, rownames(gsea$LE_matrix))
    LEmat <- LEmat <- gsea[kk, ]
  }
  rnk <- as.matrix(read.table(gsea$rnk_file, row.names = 1))[, 1]
  nes <- gsea$report[, "NES"]
  names(nes) <- gsea$report$NAME
  LEmat[which(is.na(LEmat))] <- 0
  LEmat <- LEmat[head(order(-rowMeans(LEmat != 0)), maxrow), , drop = FALSE]
  LEmat <- LEmat[, head(order(-colMeans(LEmat != 0)), maxcol), drop = FALSE]
  LEmat <- LEmat * abs(nes[rownames(LEmat)])
  LEmat <- t(t(LEmat) * (rnk[colnames(LEmat)]))
  LEmat <- abs(LEmat)**0.5 * sign(LEmat)
  LEmat[which(is.na(LEmat))] <- 0
  if (nrow(LEmat) <= 1 || ncol(LEmat) <= 1) {
    return(NULL)
  }

  if (render == "d3heatmap") {
    d3heatmap(LEmat, yaxis_width = 600)
  } else if (render == "heatmaply") {
    heatmaply(LEmat)
  } else {
    gx.heatmap(LEmat - 1e-8,
      scale = "none",
      mar = c(8, 20), cexCol = 0.9, cexRow = 0.9,
      keysize = 0.4, key = FALSE
    )
  }
}

## ========================================================================
## =========================== PLOTTING ===================================
## ========================================================================


#' @export
.bluered <- function(n = 64) {
  gplots::colorpanel(n, "royalblue3", "grey90", "indianred3")
}



#' GSEA Barplot
#'
#' Generates a barplot visualization of GSEA scores.
#'
#' @param scores A numeric vector of GSEA scores.
#' @param names A character vector of names for the scores.
#' @param xlab The label for the x-axis.
#' @param xlim The limits for the x-axis.
#' @param cex.text The size of the text labels.
#' @param main The main title of the plot.
#' @param n The number of scores to display in the barplot.
#' @export
#' @return barplot
gsea.barplot <- function(scores, names = NULL, xlab = "score", xlim = NULL,
                         cex.text = 1, main = "enrichment", n = 16) {
  if (!is.null(names)) names(scores) <- names
  scores <- rev(Matrix::head(scores[order(-abs(scores))], n))
  col1 <- c("lightskyblue1", "rosybrown1")[1 + 1 * (scores > 0)]
  if (min(scores, na.rm = TRUE) == 0 && max(scores, na.rm = TRUE) == 0) {
    xlim <- c(0, 1)
  }
  barplot(abs(scores),
    horiz = TRUE, las = 1, width = 5 / 6, col = col1,
    border = NA, xlab = xlab, names.arg = rep("", length(scores)),
    xlim = xlim, cex.axis = 0.9
  )
  title(main)
  gs <- names(scores)
  mx <- max(abs(scores))
  text(0.02 * mx, 1:length(scores) - 0.45, gs,
    adj = 0, offset = 0,
    cex = cex.text
  )
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
                        decreasing = TRUE, cex = 1, cex.main = 0.9, len.main = 40,
                        lab.line = c(0.8, 2), cex.lab = 0.8, main.line = 0.3,
                        xlab = "Rank in ordered dataset", res = 1200,
                        ylab = "Rank metric") {
  if (!is.null(names)) names(rnk) <- names
  rnk <- rnk[!is.na(rnk)]
  rnk <- rnk[order(rnk + 1e-8 * rnorm(length(rnk)), decreasing = decreasing)]

  ## ranked list metric
  ii <- (1:length(rnk))
  if (length(ii) > res) ii <- ii[seq(1, length(ii), length(ii) / res)]
  qq <- quantile(rnk[ii], probs = c(0.01, 0.99), na.rm = TRUE)
  y1 <- qq[2]
  y0 <- qq[1]
  dy <- 0.25 * (y1 - y0)
  plot(ii, rnk[ii],
    type = "h", col = "grey", ylim = c(y0 - dy, y1),
    xlab = NA, ylab = NA, xaxt = "n"
  )
  mtext(xlab, 1, line = lab.line[1], cex = cex.lab)
  mtext(ylab, 2, line = lab.line[2], cex = cex.lab)
  abline(h = 0, lty = 2, lwd = 0.5)

  ## gene set barcode
  jj <- match(gset, names(rnk))
  w1 <- ifelse(length(jj) < 100, 0.6, 0.3)
  w1 <- ifelse(length(jj) < 50, 1, w1)
  arrows(jj, (y0 - dy), jj, y0, col = "grey10", lwd = w1 * cex, length = 0)

  ## red/blue bar at bottom
  kk <- c(seq(1, length(rnk) * 0.99, floor(length(rnk) / 20)), length(rnk))
  length(kk)
  i <- 1
  for (i in 1:(length(kk) - 1)) {
    r <- mean(rnk[kk[c(i, i + 1)]])
    r1 <- (r / max(abs(rnk), na.rm = TRUE))
    r1 <- abs(r1)**0.5 * sign(r1)
    irnk <- floor(31 * (1 + r1) / 2)
    cc <- gplots::bluered(32)[1 + irnk]
    rect(kk[i], y0 - 1.05 * dy, kk[i + 1], y0 - 0.65 * dy, col = cc, border = NA)
  }

  ## weighted cumulative random walk
  x0 <- 1 * (names(rnk) %in% gset)
  x0 <- x0 * abs(rnk)
  n0 <- sum(!(names(rnk) %in% gset))
  n1 <- sum(names(rnk) %in% gset)
  r0 <- cumsum(x0 == 0) / sum(x0 == 0)
  r1 <- cumsum(x0) / (1e-4 + sum(x0))
  rnk.trace <- (r1 - r0)
  rnk.trace <- rnk.trace / max(abs(rnk.trace)) * 0.9
  if (max(rnk.trace) >= abs(min(rnk.trace))) rnk.trace <- rnk.trace * abs(y1)
  if (max(rnk.trace) < abs(min(rnk.trace))) rnk.trace <- rnk.trace * abs(y0)
  if (!decreasing) rnk.trace <- -1 * rnk.trace
  lines(ii, rnk.trace[ii], col = "green", type = "l", lwd = 2.4)

  if (is.null(main)) main <- "Enrichment plot"
  tt.main <- as.character(main)
  if (nchar(tt.main) > len.main) {
    tt.main < breakstring(tt.main, len.main) ## pgx-funtions.R
  }
  title(main = tt.main, cex.main = cex.main, line = main.line)
}


#' GSEA Enrichment Plot (Up/Down Genes)
#'
#' Generates an enrichment plot for GSEA analysis with separate traces for up-regulated and down-regulated genes.
#'
#' @param rnk A numeric vector representing the ranked list.
#' @param gset.up A character vector specifying the gene set for up-regulated genes.
#' @param gset.dn A character vector specifying the gene set for down-regulated genes.
#' @param names A character vector of names for the ranked list.
#' @param main The main title of the plot.
#' @param decreasing Logical indicating whether the ranked list should be sorted in decreasing order.
#' @param cex.main The size of the main title text.
#' @param len.main The maximum number of characters for the main title before line breaks are applied.
#' @param sum.trace Logical indicating whether to plot the sum of the traces for up-regulated and down-regulated genes.
#' @param res The resolution of the plot.
#' @param xlab The label for the x-axis.
#' @param ylab The label for the y-axis.
#' @export
#' @return a plot
gsea.enplot.UPDN <- function(rnk, gset.up, gset.dn, names = NULL, main = NULL,
                             decreasing = TRUE, cex.main = 0.9, len.main = 40,
                             sum.trace = TRUE, res = 1200,
                             xlab = "Rank in ordered dataset",
                             ylab = "Ranked list metric") {
  if (!is.null(names)) names(rnk) <- names
  rnk <- rnk[!is.na(rnk)]
  rnk <- rnk[order(rnk + 1e-8 * rnorm(length(rnk)), decreasing = decreasing)]

  ## ranked list metric
  ii <- (1:length(rnk))
  if (length(ii) > res) ii <- ii[seq(1, length(ii), length(ii) / res)]
  qq <- quantile(rnk[ii], probs = c(0.01, 0.99), na.rm = TRUE)
  y1 <- qq[2]
  y0 <- qq[1]
  dy <- 0.25 * (y1 - y0)
  plot(ii, rnk[ii],
    type = "h", col = "grey", ylim = c(y0 - 1.25 * dy, y1),
    xlab = NA, ylab = ylab, xaxt = "n"
  )
  mtext(xlab, 1, line = 1.2, cex = 0.7)
  abline(h = 0, lty = 2, lwd = 0.5)

  ## gene set barcode
  jj <- match(gset.dn, names(rnk))
  arrows(jj, (y0 - dy), jj, y0 - 0.5 * dy, col = "grey30", lwd = 1, length = 0)
  jj <- match(gset.up, names(rnk))
  arrows(jj, (y0 - 0.5 * dy), jj, y0, col = "grey30", lwd = 1, length = 0)

  ## color legend
  kk <- c(seq(1, length(rnk) * 0.95, floor(length(rnk) / 10)), length(rnk))
  i <- 1
  for (i in 1:(length(kk) - 1)) {
    r <- mean(rnk[kk[c(i, i + 1)]])
    r1 <- (r / max(abs(rnk), na.rm = TRUE))
    r1 <- abs(r1)**0.66 * sign(r1)
    irnk <- round(10 + 10 * r1)
    cc <- gplots::bluered(20)[irnk]
    rect(kk[i], y0 - 1.30 * dy, kk[i + 1], y0 - 1.05 * dy, col = cc, border = NA)
  }

  ## weighted cumulative random walk (down genes)
  x0 <- 1 * (names(rnk) %in% gset.dn)
  x0 <- x0 * abs(rnk)
  n0 <- sum(!(names(rnk) %in% gset.dn))
  n1 <- sum(names(rnk) %in% gset.dn)
  r0 <- cumsum(x0 == 0) / sum(x0 == 0)
  r1 <- cumsum(x0) / (1e-4 + sum(x0))
  trace.dn <- (r1 - r0)
  trace.dn <- trace.dn / max(abs(trace.dn)) * 0.8
  if (max(trace.dn) >= abs(min(trace.dn))) trace.dn <- trace.dn * abs(y1)
  if (max(trace.dn) < abs(min(trace.dn))) trace.dn <- trace.dn * abs(y0)
  if (!decreasing) trace.dn <- -1 * trace.dn
  lines(ii, trace.dn[ii], col = "orange", type = "l", lwd = 2.4)

  ## weighted cumulative random walk
  x0 <- 1 * (names(rnk) %in% gset.up)
  x0 <- x0 * abs(rnk)
  n0 <- sum(!(names(rnk) %in% gset.up))
  n1 <- sum(names(rnk) %in% gset.up)
  r0 <- cumsum(x0 == 0) / sum(x0 == 0)
  r1 <- cumsum(x0) / (1e-4 + sum(x0))
  trace.up <- (r1 - r0)
  trace.up <- trace.up / max(abs(trace.up)) * 0.8
  if (max(trace.up) >= abs(min(trace.up))) trace.up <- trace.up * abs(y1)
  if (max(trace.up) < abs(min(trace.up))) trace.up <- trace.up * abs(y0)
  if (!decreasing) trace.up <- -1 * trace.up
  lines(ii, trace.up[ii], col = "green", type = "l", lwd = 2.4)

  if (sum.trace) lines(ii, trace.up[ii] + trace.dn[ii], col = "blue", type = "l", lwd = 2.4)

  if (is.null(main)) main <- "Enrichment plot"
  tt.main <- as.character(main)
  if (nchar(tt.main) > len.main) {
    tt.main <- paste(substr(tt.main, 1, len.main),
      substr(tt.main, len.main + 1, 2 * len.main),
      substr(tt.main, 2 * len.main + 1, 3 * len.main),
      sep = "\n"
    )
  }
  title(main = tt.main, cex.main = cex.main, line = 0.6)
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
    rc <- fastcluster::hclust(as.dist(1 - stats::cor(t(x), use = "pairwise")))
    rc.order <- rc$order
  } else if (rsort == "metric") {
    rc.order <- order(-rft * order.sign)
  } else {
    rc.order <- 1:length(rft)
  }
  if (csort == "hclust") {
    cc <- fastcluster::hclust(as.dist(1 - stats::cor(x, use = "pairwise")))
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
    par(mar = c(0.5, 0, mar[3], mar[4] + 0))
    plot(as.dendrogram(cc), leaflab = "none", yaxt = "n")
  } else {
    par(mar = c(0.5, 0, mar[3], 0))
    frame()
  }
  title(main, cex.main = cex.main * 1.5, line = 0)

  ## column-feature barplot
  par(mar = c(1, 1.4, 0, mar[4] + 1.4))
  k0 <- "grey"
  barplot(cft[cc.order],
    width = 5 / 6, las = 3,
    xaxs = "i", add = FALSE, offset = 0, col = k0, cex.names = 1,
    names.arg = "", cex.axis = 1.0, border = FALSE,
    yaxt = "n", ylim = range(cft)
  )
  axis(4, las = 1)
  mtext(cft.lab, 4, padj = 4, cex = 0.8)

  ## heatmap
  par(mar = c(mar[1], 1.4, 1, mar[4] + 1.4))
  xx <- x[rc.order, cc.order]
  if (scale == "row") xx <- t(scale(t(xx), center = TRUE))
  if (scale == "col") xx <- scale(xx, center = TRUE)
  zlim1 <- c(-1, 1) * max(abs(xx**p), na.rm = TRUE)
  Matrix::image(1:ncol(xx), 1:nrow(xx), t(abs(xx)**p * sign(xx)),
    col = col, zlim = zlim1, xaxt = "n", xlab = "", ylab = "", yaxt = "n"
  )
  if (!is.null(xnames)) {
    mtext(xnames[rc.order], 4,
      at = (1:length(xnames)), adj = 0,
      padj = 0.5, las = 1, cex = cexRow, line = 0.5
    )
  }
  if (!is.null(ynames)) {
    mtext(ynames[cc.order], 1,
      at = (1:length(ynames)) - 0.5, adj = 1,
      padj = 1, las = 3, cex = cexCol, line = 0.5
    )
  }

  ## sample metrix barplot (vertical)
  par(mar = c(mar[1], 0, 1, 0))
  barplot(rft[rc.order],
    horiz = TRUE, las = 1, yaxs = "i", names.arg = "",
    border = FALSE, xlab = "", xlim = rev(range(rft))
  )
  mtext(rft.lab, 1, padj = 4, cex = 0.8)

  ## left dendrogram
  if (rsort == "hclust") {
    par(mar = c(mar[1], 2, 0, 0.5))
    plot(as.dendrogram(rc),
      horiz = TRUE, leaflab = "none",
      yaxt = "n", yaxs = "i"
    )
  } else {
    par(mar = c(0, 0, 0, 0))
    frame()
  }

  ## legend
  if (!is.null(legend)) {
    par(mar = c(0, mar[2], mar[3], 0) + 0.5)
    frame()
    legend("topleft", legend = legend, cex = cex.legend)
  }

  ## reset
  par(mar = c(4, 4, 2, 2))
  if (layout.reset) {
    par(mfrow = c(1, 1))
    par(mar = c(4, 4, 2, 2))
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
  rr <- lapply(gg, read.csv, sep = "\t", header = TRUE, row.names = 1, nrow = -2000)
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
  M <- pmax(M / max(M), offset)
  M <- rbind(M, 1)
  colnames(M) <- names(values)
  rownames(M) <- c("neg", "pos", "unit")
  jj <- 1:ncol(M)
  if (clust) jj <- fastcluster::hclust(dist(t(M)))$order
  M <- M[, jj]
  ii <- 1:3
  if (sum(mx < 0) == 0) ii <- c(2, 3)
  c0 <- c(col, NA)
  stars(M[ii, , drop = FALSE],
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
