# run.GSEA
run_gsea <- function(X, y, gmt, output.dir = NULL, fdr = 0.25, set.min = 15,
                     set.max = 500, topgs = 100, nperm = 1000, permute = "phenotype",
                     scoring = "weighted", do.leading.edge = FALSE, rpt.label = "gsea",
                     sort.y = "ref.last", ref.type = c("0", "NT", "REF", "DMSO", "WT", "LOW"),
                     make.sets = TRUE, clean.files = TRUE, xmx = 10, force.permute = FALSE,
                     metric = c("Signal2Noise", "Diff_of_Classes", "Pearson"),
                     gsea.program = GSEA.JAR, quiet = FALSE) {
  if (0) {
    fdr <- 0.25
    set.min <- 15
    set.max <- 500
    output.dir <- "test/"
    topgs <- 100
    nperm <- 100
    permute <- "phenotype"
    xmx <- 16
    rpt.label <- "gsea"
    sort.y <- "ref.last"
    ref.type <- c("0", "NT", "REF", "DMSO", "WT", "LOW")
    make.sets <- TRUE
    clean.file <- TRUE
    scoring <- "weighted"
    gsea.program <- GSEA.JAR
    quiet <- FALSE
  }

  ## gsea.program <- "/opt/GSEA/gsea2-2.0.13.jar"
  ## gsea.program <- "/opt/GSEA/gsea2-2.2.4.jar"
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
    ## output.dir = paste(getwd(),output.dir,sep="/")
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
  ## rpt.file = paste(output.dir,"gsea_report_for_",rpt.label,".txt",sep="")
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
    ## write("#correlation", file=cls.file, append=TRUE)
    suppressWarnings(write(paste0("#", rpt.label2), file = cls.file, append = TRUE))
    suppressWarnings(write(paste(y, collapse = " "), file = cls.file, append = TRUE))
  }

  ## prepare geneset file (cleanup names, GSEA does not like them...)
  gmt <- gmt[which(!is.na(names(gmt)))]
  gmt <- gmt[which(!(names(gmt) %in% c("", NA, "NA", " ")))]
  gmt <- gmt[order(-sapply(gmt, length))]
  ## gmt <- gmt[!duplicated(names(gmt))]
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
  ## gmt.sets <- sapply( altname, function(s) paste(s,"\tNA\t",paste(gmt[[s]],collapse="\t"),sep=""))
  gmt.sets <- sapply(altname, function(s) paste(s, "\t", s, "\t", paste(gmt[[s]], collapse = "\t"), sep = ""))
  write(gmt.sets, file = gmt.file)
  if (nperm < 0) {
    return(NULL)
  }
  if (!quiet) cat("analyzing", length(gmt.sets), "gene sets and", nrow(X), "genes\n")

  ## define proper metrics
  ## metric = "Signal2Noise"
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
  ## gsea.par <- sub("%RNKFILE",rnk.file,gsea.par)
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
    ## res$GS[jj] <- sapply(strsplit(rjj,split=":"),"[",2)
  }

  ## filter on FDR, reorder
  res <- res[order(abs(res$NES), decreasing = TRUE), ]
  if (TRUE && fdr < 1) {
    if ("FDR.q.val" %in% colnames(res)) jj <- which(res$FDR.q.val <= fdr)
    if ("FDR q-val" %in% colnames(res)) jj <- which(res[, "FDR q-val"] <= fdr)
    res <- res[jj, ]
  }

  ## annotate with leading edge genes
  if (TRUE && nrow(res) > 0) {
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
    gsea_leading_edge_analysis(output.dir)
  }

  ## cleanup
  if (clean.files) {
    unlink(gmt.file)
    ## unlink(rnk.file)
    ## unlink(cls.file)
    unlink(gct.file)
    if (substring(output.dir, 1, 4) == "/tmp") {
      ## system(paste("rm -r",output.dir))
      unlink(output.dir, recursive = TRUE, force = TRUE)
    }
  }
  return(res)
}

# gsea.LeadingEdgeAnalysis
gsea_leading_edge_analysis <- function(output.dir, ntop = 100, gsea.program = GSEA.JAR, xmx = 10) {
  ## gsea.type="Gsea";ntop=20;xmx=16
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
  ## nset=20
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
  run_gsea_leading_edge(
    rpt.path = rpt.dir, gmt = gmt.top,
    rpt.label = rpt.label, xmx = xmx,
    gsea.program = gsea.program,
    outdir = output.dir
  )
}

# run.GSEA.LeadingEdge
run_gsea_leading_edge <- function(rpt.path, gmt, xmx = 10, rpt.label = "gsea_leadingedge",
                                  gsea.program = GSEA.JAR, outdir = ".") {
  gmt0 <- gmt
  ## gmt0 <- lapply(gmt0, function(s) paste("'",s,"'",sep=""))
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


# run.GSEA.preranked
run_gsea_preranked <- function(rnk, gmt, output.dir = NULL, fdr = 0.25,
                               set.min = 15, set.max = 500, topgs = 100, nperm = 1000,
                               rpt.label = "preranked", make.sets = TRUE,
                               clean.files = TRUE, xmx = 10, scoring = "weighted",
                               chip = "GENE_SYMBOL.chip", collapse = FALSE,
                               do.leading.edge = FALSE, gsea.program = GSEA.JAR,
                               quiet = FALSE) {
  if (0) {
    output.dir <- "test"
    fdr <- 1
    set.min <- 15
    set.max <- 500
    topgs <- 20
    nperm <- 1000
    rpt.label <- "test"
    make.sets <- TRUE
    xmx <- 16
    scoring <- "weighted"
    do.leading.edge <- TRUE
  }
  ## gsea.program <- "/opt/GSEA/gsea2-2.0.13.jar"
  ## gsea.program <- "/opt/GSEA/gsea2-2.2.4.jar"
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
  ## gmt <- gmt[!duplicated(names(gmt))]
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
  # gsea.par <- sub("%COLLAPSE",tolower(collapse,gsea.par))
  gsea.par <- sub("%COLLAPSE", collapse, gsea.par)
  gsea.cmd <- paste(gsea.cmd, gsea.par)

  ## run GSEA
  out.file <- paste(output.dir, "gsea_run.out", sep = "")
  if (!quiet) cat("executing GSEA...\n", gsea.cmd, "\n")
  ## system(gsea.cmd)
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
    ## res$GS[jj] <- sapply(strsplit(rjj,split=":"),"[",2)
  }

  ## filter on FDR, reorder
  res <- res[order(abs(res$NES), decreasing = TRUE), ]
  if (TRUE && fdr < 1) {
    if ("FDR.q.val" %in% colnames(res)) jj <- which(res$FDR.q.val <= fdr)
    if ("FDR q-val" %in% colnames(res)) jj <- which(res[, "FDR q-val"] <= fdr)
    res <- res[jj, ]
  }

  ## annotate with leading edge genes
  if (TRUE && nrow(res) > 0) {
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
    gsea_leading_edge_analysis(output.dir, gsea.program = gsea.program)
  }

  ## cleanup big files
  if (clean.files) {
    unlink(gmt.file)
    ## unlink(rnk.file)
    if (substring(output.dir, 1, 4) == "/tmp") {
      ## system(paste("rm -r",output.dir))
      unlink(output.dir, recursive = TRUE, force = TRUE)
    }
  }
  return(res)
}
