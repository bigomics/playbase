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
  genes <- c(Matrix::head(names(sig), ntop), utils::tail(names(sig), ntop))

  ## Read the H5 matrix file

  if (verbose) message("[pgx.testTCGAsurvival] extracting expression from H5 matrix file")

  h5.samples <- rhdf5::h5read(matrix_file, "/meta/gdc_cases.submitter_id")
  h5.genes <- rhdf5::h5read(matrix_file, "/meta/genes")
  h5.project <- rhdf5::h5read(matrix_file, "/meta/gdc_cases.project.project_id")
  rhdf5::h5ls(matrix_file)

  sample_index <- 1:length(h5.samples)
  gene_index <- 1:length(h5.genes)


  gene_index <- which(h5.genes %in% genes)
  expression <- rhdf5::h5read(
    matrix_file, "data/expression",
    index = list(gene_index, sample_index)
  )
  colnames(expression) <- h5.samples[sample_index]
  rownames(expression) <- h5.genes[gene_index]
  expression <- log2(1 + expression)

  genes <- intersect(genes, rownames(expression))
  expression <- expression[genes, ]

  ## Read the survival data
  if (verbose) message("[pgx.testTCGAsurvival] reading TCGA survival data...")

  surv <- playdata::RTCGA_SURVIVAL
  surv$months <- round(surv$times / 365 * 12, 2)
  surv$status <- surv$patient.vital_status

  ## conform expression and surv matrices
  samples <- intersect(colnames(expression), rownames(surv))
  expression <- expression[, samples, drop = FALSE]
  surv <- surv[samples, ]

  if (verbose) {
    message("[pgx.testTCGAsurvival] dim(expression) = ", dim(expression))
    message("[pgx.testTCGAsurvival] dim(surv) = ", dim(surv))
  }

  ## print KM survival plot for each study/cancertype
  all.studies <- sort(unique(surv$cancer_type))
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
    poscor <- (rho > stats::median(rho, na.rm = TRUE))
    table(poscor)


    sdf <- survival::survdiff(survival::Surv(months, status) ~ poscor, data = sel.data)
    p.val <- 1 - stats::pchisq(sdf$chisq, length(sdf$n) - 1)
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
    surv.q <- stats::p.adjust(surv.p)
    names(surv.q) <- names(surv.p)

    i <- 1
    nplots <- length(surv.p)
    nc <- 7
    graphics::par(mfrow = c(5, nc), mar = c(3.0, 3, 1.6, 0) * 0.9, oma = c(3, 3, 0, 0), mgp = c(1.7, 0.7, 0))
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
      poscor <- (rho > stats::median(rho, na.rm = TRUE))
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
      graphics::legend("bottomleft", legend.labs,
        pch = "-", lwd = 2, col = 2:3,
        cex = 0.9, y.intersp = 0.85
      )

      p.val <- round(surv.p[study], 3)
      q.val <- round(surv.q[study], 3)
      pq <- c(paste("p=", p.val), paste("q=", q.val))
      graphics::legend("topright", pq, bty = "n", cex = 0.9, y.intersp = 0.85)
    } ## end of for
  }
  return(surv.p)
}

cancertype <- "dlbc"
variables <- "OS_"
cancertype <- "brca_tcga_pub"






