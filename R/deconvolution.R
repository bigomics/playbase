# pgx.multipleDeconvolution
pgx_multiple_deconvolution <- function(counts, refmat, methods = DECONV.METHODS) {
  methods
  timings <- c()
  results <- list()

  if (class(refmat) != "list") {
    refmat <- list("reference" = refmat)
  }

  refnames <- names(refmat)
  refnames
  i <- 1
  for (i in 1:length(refmat)) {
    message("[pgx.multipleDeconvolution] computing for ", refnames[i])
    ref <- refmat[[i]]
    dim(ref)
    res <- pgx_deconvolution(counts, ref = ref, methods = methods)

    if (!is.null(res)) {
      names(res)
      names(res$results)
      m <- names(refmat)[i]
      results[[m]] <- res$results
      timings <- rbind(timings, res$timings)
    }
  }
  timings
  timings0 <- NULL
  if (!is.null(timings)) {
    timings0 <- apply(timings, 2, function(x) tapply(x, rownames(timings), sum))
  }
  res2 <- list(results = results, timings = timings0)
  return(res2)
}

# pgx.deconvolution
pgx_deconvolution <- function(X, ref,
                              methods = c(
                                "I-NNLS", "CIBERSORT", "DCQ",
                                "DeconRNAseq", "EPIC", "NNLM",
                                "cor", "SingleR"
                              ),
                              add.unknown = FALSE, normalize.mat = TRUE) {
  if (max(X) < 50 || min(X) < 0) {
    cat("WARNING:: pgx.deconvolution: is X really counts? (not logarithmic)\n")
  }

  dbg("[pgx.deconvolution] called")


  ## clean up matrix, remove duplicate names
  mat <- as.matrix(X)
  rownames(mat) <- gsub(".*:", "", rownames(mat)) ## strip prefix
  rownames(mat) <- toupper(rownames(mat)) ## handle mouse??
  mat <- mat[order(-rowMeans(mat)), , drop = FALSE]
  mat <- as.matrix(mat[!duplicated(rownames(mat)), , drop = FALSE])

  ref <- as.matrix(ref)
  rownames(ref) <- toupper(rownames(ref))
  ref <- ref[order(-rowMeans(ref)), , drop = FALSE]
  ref <- as.matrix(ref[!duplicated(rownames(ref)), , drop = FALSE])

  dbg("[pgx.deconvolution] 1:")

  ## Add "unknown" class to reference matrix
  if (add.unknown) {
    dbg("[pgx.deconvolution] adding unknown group")
    gg <- intersect(rownames(ref), rownames(mat))
    x1 <- log(1 + ref[gg, , drop = FALSE])
    y1 <- log(1 + rowMeans(mat[gg, , drop = FALSE]))
    x1 <- cbind(offset = 1, x1)

    ## compute residual matrix by substracting all possible linear
    ## combinations of reference.

    x1 <- ref[gg, , drop = FALSE]
    y1 <- rowMeans(mat[gg, , drop = FALSE])
    cf <- NNLM::nnlm(x1, cbind(y1))$coefficients
    cf[is.na(cf)] <- 0
    resid <- pmax(y1 - x1 %*% cf, 0) ## residual vector
    resx <- rep(0, nrow(ref))
    names(resx) <- rownames(ref)
    resx[gg] <- resid
    ref <- cbind(ref, "other_cells" = resx)
  }

  ## normalize all matrices to CPM
  if (normalize.mat) {
    ref <- t(t(ref) / sqrt(1e-6 + Matrix::colSums(ref**2, na.rm = TRUE))) * 1e6
    mat <- t(t(mat) / sqrt(1e-6 + Matrix::colSums(mat**2, na.rm = TRUE))) * 1e6
  }

  ## add small noise, some methods need it...
  ref <- ref + 1e-2 * matrix(rnorm(length(ref)), nrow(ref), ncol(ref))
  mat <- mat + 1e-2 * matrix(rnorm(length(mat)), nrow(mat), ncol(mat))
  ref <- pmax(ref, 0)
  mat <- pmax(mat, 0)

  dbg("[pgx.deconvolution] 2:")

  gg <- intersect(rownames(ref), rownames(mat))
  Matrix::head(gg)
  if (length(gg) < 10) {
    message("WARNING:: pgx.deconvolution: no enough marker genes")
    return(NULL)
  }

  ## conform??
  if (0) {
    length(gg)
    ## ref <- ref[gg,]
    ## mat <- mat[gg,]
    qx <- limma::normalizeQuantiles(cbind(ref[gg, ], mat[gg, ]))
    ref <- qx[, colnames(ref)]
    mat <- mat[, colnames(mat)]
  }

  timings <- list()
  results <- list()
  methods
  CIBERSORT.code <- "/opt/CIBERSORT/CIBERSORTmat.R"
  if ("CIBERSORT" %in% methods && file.exists(CIBERSORT.code)) {
    ## CIBERSORT
    ## source(file.path(FILES,CIBERSORT.code))
    source(CIBERSORT.code)
    cat("starting deconvolution using CIBERSORT...\n")
    ciber.out <- NULL
    stime <- system.time(
      ## try( ciber.out <- CIBERSORT(ref, mat, perm=0, QN=TRUE) )
      try(ciber.out <- CIBERSORT(ref, mat, perm = 0, QN = FALSE))
    )
    if (!is.null(ciber.out)) {
      timings[["CIBERSORT"]] <- stime
      cat("deconvolution using CIBERSORT took", stime[3], "s\n")
      ciber.out <- ciber.out[, !(colnames(ciber.out) %in% c("P-value", "Correlation", "RMSE"))]
      results[["CIBERSORT"]] <- ciber.out
    } else {
      cat("WARNING:: CIBERSORT failed\n")
    }
  }

  if ("EPIC" %in% methods) {
    ## EPIC
    ## devtools::install_github("GfellerLab/EPIC", build_vignettes=TRUE)
    dbg("[pgx.deconvolution] calculating EPIC...")

    out <- NULL
    gg <- intersect(rownames(ref), rownames(mat))
    ref1 <- ref
    if ("other_cells" %in% colnames(ref1)) {
      ref1 <- ref1[, setdiff(colnames(ref1), "other_cells")]
    }
    ref.list <- list(refProfiles = ref1, sigGenes = gg)
    mat1 <- mat
    colnames(mat1) <- 1:ncol(mat) ## EPIC doesnt like duplicated column names...
    stime <- system.time(try(
      out <- EPIC(bulk = mat1, reference = ref.list)
    ))
    remove(mat1)
    if (!is.null(out)) {
      message(paste0("deconvolution using EPIC took", stime[3], "s\n"))
      timings[["EPIC"]] <- stime
      out.mat <- out[["cellFractions"]]
      colnames(out.mat) <- sub("otherCells", "other_cells", colnames(out.mat))
      rownames(out.mat) <- colnames(mat)
      results[["EPIC"]] <- out.mat
    } else {
      mssage("WARNING:: EPIC fail (no pun intended...)\n")
    }
  }

  if ("DeconRNAseq" %in% methods) {
    ## DeconRNAseq
    ## BiocManager::install("DeconRNASeq", version = "3.8")
    ## BiocManager::install("DeconRNASeq")
    if ("package:Seurat" %in% search()) detach("package:Seurat", unload = TRUE)
    dbg("[pgx.deconvolution] calculating DeconRNAseq...")

    ## uses psych::pca() from pcaMethods
    require(pcaMethods) ## uses pcaMethods::prep
    drs <- NULL
    stime <- system.time(suppressMessages(suppressWarnings(
      drs <- try(DeconRNASeq::DeconRNASeq(
        data.frame(mat, check.names = FALSE),
        data.frame(ref, check.names = FALSE)
      )$out.all)
    )))
    timings[["DeconRNAseq"]] <- stime
    class(drs)
    if (!is.null(drs) && class(drs) != "try-error") {
      cat("deconvolution using DeconRNAseq took", stime[3], "s\n")
      rownames(drs) <- colnames(mat)
      ## colnames(drs) <- colnames(ref)
      results[["DeconRNAseq"]] <- drs
    } else {
      cat("WARNING:: DeconRNAseq failed\n")
    }
  }

  ## --------- DCQ from ComICS ---------------------
  if ("DCQ" %in% methods) {
    ## DCQ seems to work in logX, so we use log-transform
    ## install.packages("ComICS")
    dbg("[pgx.deconvolution] calculating DCQ...")

    res.dcq <- NULL
    stime <- system.time(
      res.dcq <- try(
        ComICS::dcq(
          reference_data = log2(1 + as.matrix(ref)),
          mix_data = log2(1 + as.matrix(mat)), ## log data OK??
          ## marker_set = matrix(rownames(ref),ncol=1),
          marker_set = cbind(intersect(rownames(ref), rownames(mat))),
          alpha_used = 0.05, lambda_min = 0.2, number_of_repeats = 3,
          precent_of_data = 1.0
        )
      )
    )
    if (!is.null(res.dcq) && class(res.dcq) != "try-error") {
      timings[["DCQ"]] <- stime
      cat("deconvolution using DCQ took", stime[3], "s\n")
      results[["DCQ"]] <- res.dcq$average
    } else {
      cat("WARNING:: DCQ failed\n")
    }
  }

  if ("I-NNLS" %in% methods) {
    ## !!!!!!!!!!!!!! WARNING:: needs more testing/validation !!!!!!!!!!!!!
    ## ----- Constrained iterative non-negative least squares (Abbas et al. 2009) ----
    dbg("[pgx.deconvolution] calculating I-NNLS...")

    GetFractions.Abbas <- function(XX, y, w = NA) {
      ## XX is immune expression data
      ## y is cancer expression data
      ss.remove <- c()
      ss.names <- colnames(XX)
      while (TRUE) {
        if (length(ss.remove) == 0) {
          tmp.XX <- XX
        } else {
          if (is.null(ncol(tmp.XX))) {
            return(rep(0, ncol(XX)))
          }
          tmp.XX <- tmp.XX[, -ss.remove, drop = FALSE]
        }
        if (length(ss.remove) > 0) {
          ss.names <- ss.names[-ss.remove]
          if (length(ss.names) == 0) {
            return(rep(0, ncol(XX)))
          }
        }
        if (is.na(w[1])) {
          tmp <- lsfit(tmp.XX, y, intercept = FALSE)
        } else {
          tmp <- lsfit(tmp.XX, y, w, intercept = FALSE)
        }
        if (ncol(tmp.XX) == 1) {
          tmp.beta <- tmp$coefficients[1]
        } else {
          tmp.beta <- tmp$coefficients[1:(ncol(tmp.XX) + 0)]
        }
        if (min(tmp.beta > 0)) break ## break if coefs are all positive
        ss.remove <- which.min(tmp.beta) ## removes most negative coeff
      }
      tmp.F <- rep(0, ncol(XX))
      names(tmp.F) <- colnames(XX)
      tmp.F[ss.names] <- tmp.beta
      return(tmp.F)
    }
    gg <- intersect(rownames(ref), rownames(mat))
    XX <- ref[gg, , drop = FALSE]
    YY <- mat[gg, , drop = FALSE]
    res.abbas <- NULL
    stime <- system.time(
      res.abbas <- try(apply(YY, 2, function(y) GetFractions.Abbas(XX, y, w = NA)))
    )
    timings[["I-NNLS"]] <- stime
    cat("deconvolution using I-NNLS took", stime[3], "s\n")
    if (!is.null(res.abbas) && class(res.abbas) != "try-error") {
      rownames(res.abbas) <- colnames(ref)
      results[["I-NNLS"]] <- t(res.abbas)
    }
  }

  ## Own NNLM (non-negative linear modeling)...
  if ("NNLM" %in% methods) {
    ## NNLM
    ## install.packages("NNLM")
    dbg("[pgx.deconvolution] calculating NNLM...")

    x1 <- log2(1 + ref[gg, , drop = FALSE])
    x2 <- log2(1 + mat[gg, , drop = FALSE])
    x1 <- cbind(offset = 1, x1)
    stime <- system.time(
      cf <- NNLM::nnlm(x1, x2)$coefficients[-1, , drop = FALSE]
    )
    timings[["NNLM"]] <- stime
    cat("deconvolution using NNLM took", stime[3], "s\n")
    results[["NNLM"]] <- t(cf)

    ## very much the same as I-NNLS
    results[["NNLM.lin"]] <- t(NNLM::nnlm(ref[gg, , drop = FALSE], mat[gg, , drop = FALSE])$coefficients)

    r1 <- apply(ref[gg, , drop = FALSE], 2, rank, na.last = "keep")
    r2 <- apply(mat[gg, , drop = FALSE], 2, rank, na.last = "keep")
    r1 <- cbind(offset = 1, r1)
    cf <- NNLM::nnlm(r1, r2)$coefficients[-1, , drop = FALSE]
    results[["NNLM.rnk"]] <- t(cf)
  }

  ## Simple (rank) correlation
  if ("cor" %in% methods) {
    dbg("[pgx.deconvolution] calculating cor...")
    ## results[["cor"]] <- stats::cor(log2(1+mat[gg,]), log2(1+ref[gg,]))
    r1 <- apply(mat[gg, , drop = FALSE], 2, rank, na.last = "keep")
    r2 <- apply(ref[gg, , drop = FALSE], 2, rank, na.last = "keep")
    stime <- system.time(
      cf <- stats::cor(r1, r2, use = "pairwise")
    )
    timings[["cor"]] <- stime
    cat("deconvolution using COR took", stime[3], "s\n")
    results[["cor"]] <- cf
  }

  if ("SingleR" %in% methods) {
    dbg("[pgx.deconvolution] calculating SingleR...")
    stime <- system.time(
      sr1 <- SingleR(test = mat, ref = ref, labels = colnames(ref))
    )
    timings[["SingleR"]] <- stime
    cat("deconvolution using SingleR took", stime[3], "s\n")
    results[["SingleR"]] <- sr1$scores
  }
  ## clean up
  names(results)
  results <- results[which(!sapply(results, is.null))]
  results <- lapply(results, function(x) {
    x[is.na(x)] <- 0
    x
  })

  dbg("[pgx.deconvolution] calculating meta values...")

  ## meta
  if (length(results) > 1) {
    jj <- colnames(ref)
    norm.results <- lapply(results, function(x) x[, jj, drop = FALSE] / (1e-8 + rowSums(x[, jj, drop = FALSE])))
    lognorm.results <- lapply(norm.results, function(x) log(0.001 + pmax(x, 0)))
    res.meta1 <- Reduce("+", norm.results) / length(norm.results)
    res.meta2 <- exp(Reduce("+", lognorm.results) / length(lognorm.results))
    results[["meta"]] <- res.meta1 / (1e-8 + rowSums(res.meta1, na.rm = TRUE))
    results[["meta.prod"]] <- res.meta2 / (1e-8 + rowSums(res.meta2, na.rm = TRUE))
  }

  timings0 <- do.call(rbind, timings)
  res2 <- list(results = results, timings = timings0)

  dbg("[pgx.deconvolution] done!")

  return(res2)
}

# pgx.inferCellCyclePhase
infer_cell_cycle_phase <- function(counts) {
  ## List of cell cycle markers, from Tirosh et al, 2015
  ##
  ## cc.genes <- readLines(con = "../opt/seurat/regev_lab_cell_cycle_genes.txt")
  cc.genes <- strsplit("MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 UHRF1 MLF1IP HELLS RFC2 RPA2 NASP RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 RAD51 RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 CHAF1B BRIP1 E2F8 HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPF TACC3 FAM64A SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B HJURP CDCA3 HN1 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA", split = " ")[[1]]
  s_genes <- cc.genes[1:43]
  g2m_genes <- cc.genes[44:97]

  ## Create our Seurat object and complete the initalization steps
  rownames(counts) <- toupper(rownames(counts)) ## mouse...
  ## counts1 <- cbind(counts,counts,counts,counts,counts,counts)
  obj <- Seurat::CreateSeuratObject(counts)
  obj <- Seurat::NormalizeData(obj, verbose = 0)
  suppressWarnings(obj <- Seurat::CellCycleScoring(obj,
    s.features = s_genes,
    g2m.features = g2m_genes, set.ident = TRUE
  ))
  ## view cell cycle scores and phase assignments
  ## head(x = obj@meta.data)
  ## table(obj@meta.data$Phase)
  s.score <- obj@meta.data$S.Score
  g2m.score <- obj@meta.data$G2M.Score
  phase <- obj@meta.data$Phase
  if (is.null(phase) || length(phase) == 0) {
    return(NULL)
  }
  return(phase)
}

# pgx.inferGender
infer_gender <- function(X, gene_name = NULL) {
  ## List of cell cycle markers, from Tirosh et al, 2015
  ##
  ## cc.genes <- readLines(con = "../opt/seurat/regev_lab_cell_cycle_genes.txt")
  if (is.null(gene_name)) gene_name <- toupper(sub(".*:", "", rownames(X)))
  if (0) {
    X <- log2(1 + ngs$counts)
    gene_name <- rownames(X)
    x.genes <- rownames(ngs$genes)[grep("X", ngs$genes$chr)]
    y.genes <- rownames(ngs$genes)[grep("Y", ngs$genes$chr)]
    x.genes
    y.genes
  } else {
    y.genes <- intersect(c("DDX3Y", "RPS4Y1", "USP9Y", "KDM5D"), gene_name)
    y.genes
    x.genes <- intersect(c("XIST"), gene_name)
    x.genes
  }
  if (length(y.genes) == 0 && length(x.genes) == 0) {
    cat("warning:: could not determine sex. missing some X/Y marker genes\n")
    sex <- rep(NA, ncol(X))
    return(sex)
  }
  sex <- rep(NA, ncol(X))
  if (length(y.genes) > 0 && length(x.genes) > 0) {
    x.expr <- colMeans(X[match(x.genes, gene_name), , drop = FALSE])
    y.expr <- colMeans(X[match(y.genes, gene_name), , drop = FALSE])
    x.expr
    y.expr
    mean.expr <- colMeans(X)
    ## plot(x.expr, y.expr)
    sex <- rep(NA, ncol(X))
    sex <- ifelse(x.expr > mean.expr & y.expr < mean.expr, "F", sex)
    sex <- ifelse(y.expr > mean.expr & x.expr < mean.expr, "M", sex)
    sex
    return(sex)
  }
  return(sex)
}
