##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' @export
pgx.testTCGAsurvival <- function(sig, matrix_file, ntop = 100, deceased.only = TRUE,
                                 min.cases = 10, sortby.p = FALSE, plot = TRUE, verbose = 1) {
  if (!file.exists(matrix_file)) {
    stop("cannot find TCGA H5 matrix file")
  }
  if (is.null(names(sig))) {
    stop("sig must have names")
  }
  ## !!!!!!!!!! upper case for mouse genes !!!!!
  names(sig) <- toupper(names(sig))

  ## get the top DE genes
  sig <- sort(sig)
  genes <- c(Matrix::head(names(sig), ntop), tail(names(sig), ntop))
  Matrix::head(genes)

  ## Read the H5 matrix file

  ## aa



  ## aa.head

  if (verbose) message("[pgx.testTCGAsurvival] extracting expression from H5 matrix file")

  h5.samples <- rhdf5::h5read(matrix_file, "/meta/gdc_cases.submitter_id")
  h5.genes <- rhdf5::h5read(matrix_file, "/meta/genes")
  h5.project <- rhdf5::h5read(matrix_file, "/meta/gdc_cases.project.project_id")
  rhdf5::h5ls(matrix_file)

  sample_index <- 1:length(h5.samples)
  gene_index <- 1:length(h5.genes)


  gene_index <- which(h5.genes %in% genes)
  Matrix::head(gene_index)
  expression <- rhdf5::h5read(
    matrix_file, "data/expression",
    index = list(gene_index, sample_index)
  )
  colnames(expression) <- h5.samples[sample_index]
  rownames(expression) <- h5.genes[gene_index]
  expression <- log2(1 + expression)
  dim(expression)

  genes <- intersect(genes, rownames(expression))
  expression <- expression[genes, ]
  dim(expression)

  ## Read the survival data
  if (verbose) message("[pgx.testTCGAsurvival] reading TCGA survival data...")

  surv <- playdata::RTCGA_SURVIVAL
  Matrix::head(surv)
  surv$months <- round(surv$times / 365 * 12, 2)
  surv$status <- surv$patient.vital_status

  ## conform expression and surv matrices
  samples <- intersect(colnames(expression), rownames(surv))
  length(samples)
  expression <- expression[, samples, drop = FALSE]
  surv <- surv[samples, ]

  if (verbose) {
    message("[pgx.testTCGAsurvival] dim(expression) = ", dim(expression))
    message("[pgx.testTCGAsurvival] dim(surv) = ", dim(surv))
  }

  ## print KM survival plot for each study/cancertype
  all.studies <- sort(unique(surv$cancer_type))
  length(all.studies)
  table(all.studies)
  study <- all.studies[1]
  study <- "BRCA"

  message("[pgx.testTCGAsurvival] fitting survival probabilities...\n")

  surv.p <- rep(NA, length(all.studies))
  rho.list <- list()
  sel.list <- list()
  names(surv.p) <- all.studies
  for (study in Matrix::head(all.studies, 99)) {
    study

    ## calculate correlation with signature
    sel <- which(surv$cancer_type == study)
    if (deceased.only) {
      sel <- which(surv$cancer_type == study & surv$status == 1) ## only died people
    }
    if (length(sel) < min.cases) next()
    gg <- intersect(rownames(expression), names(sig))
    sel.X <- expression[gg, sel] - rowMeans(expression[gg, sel], na.rm = TRUE)
    rho <- stats::cor(sel.X, sig[gg], use = "pairwise")[, 1]
    sel.data <- surv[sel, ]

    ## fit survival curve on two groups
    poscor <- (rho > median(rho, na.rm = TRUE))
    table(poscor)


    sdf <- survival::survdiff(survival::Surv(months, status) ~ poscor, data = sel.data)
    p.val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
    p.val

    surv.p[study] <- p.val
    rho.list[[study]] <- rho
    sel.list[[study]] <- rownames(sel.data)
  }
  surv.p <- surv.p[names(rho.list)]

  if (plot) {
    message("[pgx.testTCGAsurvival] plotting KM curves...\n")

    if (sortby.p) {
      ii <- order(surv.p)
      surv.p <- surv.p[ii]
      rho.list <- rho.list[ii]
    }
    surv.q <- p.adjust(surv.p)
    names(surv.q) <- names(surv.p)

    i <- 1
    nplots <- length(surv.p)
    nc <- 7
    par(mfrow = c(5, nc), mar = c(3.0, 3, 1.6, 0) * 0.9, oma = c(3, 3, 0, 0), mgp = c(1.7, 0.7, 0))
    i <- 1
    for (i in 1:nplots) {
      study <- names(surv.p)[i]
      study

      ## calculate correlation with signature
      sel <- sel.list[[study]]
      sel.data <- surv[sel, ]
      rho <- rho.list[[study]]
      rho <- rho[rownames(sel.data)]

      ## fit survival curve on two groups
      poscor <- (rho > median(rho, na.rm = TRUE))
      table(poscor)

      fit <- survival::survfit(survival::Surv(months, status) ~ poscor, data = sel.data)


      legend.labs <- paste(c("rho<0", "rho>0"))
      legend.labs <- paste(c("neg.cor", "pos.cor"))

      xlab <- ylab <- ""
      is.first <- (i %% nc == 1)
      if (is.first) ylab <- "survival probability"

      last.row <- ((i - 1) %/% 7 == (nplots - 1) %/% 7)
      if (last.row) xlab <- "time (days)"
      last.row

      plot(fit,
        col = 2:3, lwd = 2, main = study,
        xlab = xlab, ylab = ylab, cex.main = 1.1
      )
      legend("bottomleft", legend.labs,
        pch = "-", lwd = 2, col = 2:3,
        cex = 0.9, y.intersp = 0.85
      )

      p.val <- round(surv.p[study], 3)
      q.val <- round(surv.q[study], 3)
      pq <- c(paste("p=", p.val), paste("q=", q.val))
      legend("topright", pq, bty = "n", cex = 0.9, y.intersp = 0.85)
    } ## end of for
  }
  return(surv.p)
}

cancertype <- "dlbc"
variables <- "OS_"
cancertype <- "brca_tcga_pub"

#' @export
pgx.selectTCGAstudies <- function(cancertype, variables) {
  ## Scan the available TCGA studies for cancertype and clinical
  ## variables.
  ##
  ##
  ##

  mycgds <- cgdsr::CGDS("http://www.cbioportal.org/")
  all.studies <- sort(cgdsr::getCancerStudies(mycgds)[, 1])
  studies <- grep(cancertype, all.studies, value = TRUE)
  clin <- list()
  samples <- list()
  studies
  mystudy <- studies[1]

  for (mystudy in studies) {
    mystudy
    myprofiles <- cgdsr::getGeneticProfiles(mycgds, mystudy)[, 1]
    myprofiles

    ## mrna datatypes
    mrna.type <- "mrna"
    if (any(grepl("rna_seq_mrna$", myprofiles))) mrna.type <- "rna_seq_mrna"
    if (any(grepl("v2_mrna$", myprofiles))) mrna.type <- "rna_seq_v2_mrna"
    pr.mrna <- grep(paste0(mrna.type, "$"), myprofiles, value = TRUE)
    pr.mrna
    if (length(pr.mrna) == 0) next()

    all.cases <- cgdsr::getCaseLists(mycgds, mystudy)[, 1]
    all.cases


    caselist <- grep(paste0(mrna.type, "$"), all.cases, value = TRUE)
    caselist
    clin0 <- cgdsr::getClinicalData(mycgds, caselist)
    Matrix::head(clin0)[, 1:4]
    rownames(clin0) <- gsub("[.]", "-", rownames(clin0)) ## correct names...
    Matrix::head(clin0)[, 1:4]
    clin[[mystudy]] <- clin0
    samples[[mystudy]] <- rownames(clin0)
  }

  sel <- sapply(clin, function(v) any(grepl(variables, colnames(v))))
  sel
  sel.studies <- studies[sel]
  sel.clin <- clin[sel]

  res <- list(
    studies = sel.studies,
    clinicalData = sel.clin
  )
  return(res)
}


#' @export
pgx.getTCGAdataset <- function(study, genes = NULL, matrix_file = NULL, from.h5 = TRUE,
                               datatype = "mrna") {
  ## For a specific TCGA study get the expression matrix and
  ## clinical data.
  ##

  ## check if H5 exists
  from.h5 <- (from.h5 && !is.null(matrix_file) && file.exists(matrix_file))



  mycgds <- cgdsr::CGDS("http://www.cbioportal.org/")

  all.studies <- sort(cgdsr::getCancerStudies(mycgds)[, 1])
  if (!all(study %in% all.studies)) {
    ss <- setdiff(study, all.studies)
    stop(ss, "is not in TCGA studies")
  }

  ## Gather data from all study
  X <- list()
  clin <- list()
  mystudy <- study[1]
  for (mystudy in study) {
    cat("getting TCGA data for", mystudy, "...\n")

    mystudy

    myprofiles <- cgdsr::getGeneticProfiles(mycgds, mystudy)[, 1]
    myprofiles

    ## datatypes
    datatype0 <- datatype
    if (datatype0 == "mrna") {
      datatype <- "_mrna"
      if (any(grepl("rna_seq_v2_mrna$", myprofiles))) datatype <- "rna_seq_v2_mrna"
      if (any(grepl("rna_seq_mrna$", myprofiles))) datatype <- "rna_seq_mrna"
      pr.datatype <- grep(paste0(datatype, "$"), myprofiles, value = TRUE)
      pr.datatype
      if (length(pr.datatype) == 0) {
        cat("WARNING:: could not find mRNA for", mystudy, "...\n")
        next()
      }
    } else {
      warning(paste("datatype", datatype0, "not yet supported"))
      return(NULL)
    }

    all.cases <- cgdsr::getCaseLists(mycgds, mystudy)[, 1]
    all.cases


    caselist <- grep(paste0(datatype, "$"), all.cases, value = TRUE)
    caselist
    samples <- NULL
    Matrix::head(genes)
    if (!is.null(genes) && !from.h5) {
      cat("downloading...\n")
      ## If only a few genes, getProfileData is a faster way
      ##
      expression <- t(cgdsr::getProfileData(mycgds, genes, pr.datatype, caselist))
      samples <- gsub("[.]", "-", colnames(expression))
      colnames(expression) <- samples
      dim(expression)
    } else {
      cat("extracting from locally stored H5 matrix...\n")
      ## For all genes, getProfileData cannot do and we use
      ## locally stored H5 TCGA data file from Archs4.
      ##
      xx <- cgdsr::getProfileData(mycgds, "---", pr.datatype, caselist)
      samples <- gsub("[.]", "-", colnames(xx))[3:ncol(xx)]
      Matrix::head(samples)



      rhdf5::h5closeAll()

      has.h5 <- file.exists(matrix_file)
      has.h5

      if (!has.h5) {
        stop("FATAL: could not find tcga_matrix.h5 matrix. Please download from Archs4.")
      } else {
        ## Retrieve information from locally stored H5 compressed data
        aa <- rhdf5::h5ls(matrix_file)[, 1:2]
        aa
        ii <- which(aa[, 1] == "/meta")[-1]
        lapply(ii, function(i) Matrix::head(rhdf5::h5read(matrix_file, paste0("/meta/", aa[i, 2]))))

        id1 <- rhdf5::h5read(matrix_file, "/meta/gdc_cases.samples.portions.submitter_id")
        id2 <- rhdf5::h5read(matrix_file, "/meta/gdc_cases.samples.submitter_id")
        id3 <- rhdf5::h5read(matrix_file, "/meta/gdc_cases.submitter_id")
        id2x <- substring(id2, 1, 15)

        h5.genes <- rhdf5::h5read(matrix_file, "/meta/genes")
        if (!is.null(genes)) h5.genes <- intersect(genes, h5.genes)
        samples <- intersect(samples, id2x)
        sample_index <- which(id2x %in% samples)
        gene_index <- 1:length(h5.genes)

        if (length(sample_index) == 0 || length(gene_index) == 0) {
          return(list(X = NULL, clin = NULL))
        }

        expression <- rhdf5::h5read(
          matrix_file, "data/expression",
          index = list(gene_index, sample_index)
        )
        rhdf5::H5close()
        dim(expression)
        colnames(expression) <- substring(id2[sample_index], 1, 15)
        rownames(expression) <- h5.genes
        expression <- expression[, order(-Matrix::colSums(expression))]
        expression <- expression[, samples]
      }
    }
    dim(expression)
    this.clin <- cgdsr::getClinicalData(mycgds, caselist)
    rownames(this.clin) <- gsub("[.]", "-", rownames(this.clin))
    this.clin <- this.clin[samples, , drop = FALSE]
    expression <- expression[, samples, drop = FALSE]
    X[[mystudy]] <- expression
    clin[[mystudy]] <- this.clin
  }

  res <- list(X = X, clin = clin)
  return(res)
}


#' @export
pgx.getTCGA.multiomics.TOBEFINISHED <- function(studies, genes = NULL, batch.correct = TRUE,
                                                tcga.only = TRUE) {
  ## Better use curatedTCGA bioconductor package!!!!
  ##



  mycgds <- cgdsr::CGDS("http://www.cbioportal.org/")

  ## for(mystudy in all.studies) {


  ## }

  GENE <- "CIITA"
  GENE <- "NLRC5"

  ## Gather data from all cancers
  all.X <- list()
  mystudy <- studies[1]
  for (mystudy in studies) {
    mystudy

    myprofiles <- cgdsr::getGeneticProfiles(mycgds, mystudy)[, 1]
    myprofiles

    ## prioritize datatypes
    pr.mrna <- grep("rna_seq_v2_mrna$|rna_seq_mrna$", myprofiles, value = TRUE)[1]

    pr.cna <- grep("_log2CNA$|_linear_CNA$", myprofiles, value = TRUE)[1]
    pr.gistic <- grep("_gistic$", myprofiles, value = TRUE)[1]
    pr.me <- grep("_methylation_hm450|_methylation_hm27", myprofiles, value = TRUE)[1]
    pr.mut <- grep("_mutations", myprofiles, value = TRUE)[1]

    all.cases <- cgdsr::getCaseLists(mycgds, mystudy)[, 1]
    all.cases
    if (!any(grepl("complete$", all.cases))) next

    caselist <- grep("complete$", all.cases, value = TRUE)
    cna <- counts <- cna.gistic <- me <- mut <- gx <- NULL
    counts <- cgdsr::getProfileData(mycgds, genes, pr.mrna, caselist)
    cna <- cgdsr::getProfileData(mycgds, GENE, pr.cna, caselist)

    cna.gistic <- cgdsr::getProfileData(mycgds, GENE, pr.gistic, caselist)
    me <- cgdsr::getProfileData(mycgds, GENE, pr.me, caselist)
    mut <- cgdsr::getProfileData(mycgds, GENE, pr.mut, caselist)
    mut <- 1 * !is.na(mut)
    gx <- log2(10 + as.matrix(counts))
    cna[is.na(cna)] <- NA
    if (grepl("linear", pr.cna)) cna <- log2(0.01 + 2 + cna) ## assume diploid


    if (!is.null(cna)) colnames(cna) <- paste0("CN:", colnames(cna))
    if (!is.null(cna.gistic)) colnames(cna.gistic) <- paste0("CNA:", colnames(cna.gistic))
    if (!is.null(me)) colnames(me) <- paste0("ME:", colnames(me))
    if (!is.null(mut)) colnames(mut) <- paste0("MT:", colnames(mut))



    xx <- list(gx, cna.gistic, me, mut)
    xx <- xx[sapply(xx, nrow) > 0]
    X <- do.call(cbind, xx)
    dim(X)

    if (!is.null(X) && ncol(X) >= 4) {
      X <- X[, colMeans(is.na(X)) < 0.5, drop = FALSE]
      X <- X[rowMeans(is.na(X)) < 0.5, , drop = FALSE]
      dim(X)
      all.X[[mystudy]] <- X
    }
  }
}

#' @export
pgx.getTCGAproteomics <- function() {
  GenomicDataCommons::status()

  qfiles <- GenomicDataCommons::files() %>% filter(~ cases.project.project_id == "TCGA-BRCA" &
    type == "gene_expression" &
    analysis.workflow_type == "HTSeq - Counts")
  manifest_df <- qfiles %>% GenomicDataCommons::manifest()
  nrow(manifest_df)
  Matrix::head(manifest_df)

  fnames <- GenomicDataCommons::gdcdata(manifest_df$id[1:2], progress = FALSE)

  resp <- GenomicDataCommons::cases() %>%
    filter(~ project.project_id == "TCGA-BRCA" &
      project.project_id == "TCGA-BRCA") %>%
    GenomicDataCommons::facet("samples.sample_type") %>%
    GenomicDataCommons::aggregations()
  resp$samples.sample_type

  resp <- GenomicDataCommons::cases() %>%
    filter(~ project.project_id == "CPTAC-3" &
      project.project_id == "CPTAC-3") %>%
    GenomicDataCommons::facet("samples.sample_type") %>%
    GenomicDataCommons::aggregations()
  resp$samples.sample_type
}
