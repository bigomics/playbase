##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


#' @title Run Salmon 
#'
#' @param srr_id Character vector of SRR IDs for samples
#' @param library_layout Library layout - "SINGLE" or "PAIRED" 
#' @param index_dir Salmon index directory 
#' @param destdir Output directory 
#' @param fastq_dir Directory containing FASTQ files
#' @param use_trimmed_fastq Use trimmed FASTQ files if available. Default is FALSE.
#' @param other_opts Other Salmon options to include
#' @param nthread Number of threads for parallel processing
#'
#' @return Salmon quantification results written to output directory
#' 
#' @description Run Salmon for RNA-seq quantification on sample FASTQ files.
#'
#' @details This function runs the Salmon tool to quantify transcript abundance from RNA-seq reads.
#' It takes a vector of SRR IDs, the Salmon index directory, and FASTQ file locations as input.
#' 
#' Salmon is run in single-end or paired-end mode based on the library_layout parameter.
#' If use_trimmed_fastq=TRUE, trimmed FASTQ files will be used if available.
#' 
#' Quantification results are written to the specified output directory.
#' Parallel processing across multiple threads can be enabled via the nthread parameter.
#'
#' @export
pgx.run_salmon <- function(
    srr_id, library_layout = c("SINGLE", "PAIRED"), index_dir,
    destdir, fastq_dir, use_trimmed_fastq = FALSE, other_opts = NULL, nthread) {
  if (!dir.exists(paste0(destdir, "/salmon"))) {
    system(paste0("mkdir ", destdir, "/salmon"))
  }
  library_layout <- match.arg(library_layout, c("SINGLE", "PAIRED"))
  if (library_layout == "SINGLE") {
    if (use_trimmed_fastq) {
      system(paste0(
        "salmon quant -i ", index_dir, " -p ",
        nthread, " ", other_opts, " -l A -r ", fastq_dir,
        "/", srr_id, "_trimmed.fastq -o ", destdir, "/salmon/",
        srr_id, "_transcripts_quant"
      ))
    } else {
      system(paste0(
        "salmon quant -i ", index_dir, " -p ",
        nthread, " ", other_opts, " -l A -r ", fastq_dir,
        "/", srr_id, ".fastq -o ", destdir, "/salmon/",
        srr_id, "_transcripts_quant"
      ))
    }
  } else {
    if (use_trimmed_fastq) {
      system(paste0(
        "salmon quant -i ", index_dir, " -p ",
        nthread, " ", other_opts, " -l A -1 ", fastq_dir,
        "/", srr_id, "_trimmed_1.fastq ", "-2 ", fastq_dir,
        "/", srr_id, "_trimmed_2.fastq -o ", destdir,
        "/salmon/", srr_id, "_transcripts_quant"
      ))
    } else {
      system(paste0(
        "salmon quant -i ", index_dir, " -p ",
        nthread, " ", other_opts, " -l A -1 ", fastq_dir,
        "/", srr_id, "_1.fastq ", "-2 ", fastq_dir,
        "/", srr_id, "_2.fastq -o ", destdir, "/salmon/",
        srr_id, "_transcripts_quant"
      ))
    }
  }
  if (file.exists(paste0(destdir, "/salmon/", srr_id, "_transcripts_quant/quant.sf"))) {
    system(paste("cat ", destdir, "/salmon/", srr_id, "_transcripts_quant/quant.sf",
      "| sed -E 's/\\.[0-9]+//' > ", destdir, "/salmon/",
      srr_id, "_transcripts_quant", "/", srr_id, "_quant_new.sf",
      sep = ""
    ))
  } else {
    cat("quant.sf doesn't exist. Processing next sample.")
  }
}


#' @title Run tximport on Kallisto quantification
#'
#' @param srr_id Character vector of SRR IDs for samples
#' @param species Species name, used to get gene annotation data. Default is c("human", "mouse", "rat").
#' @param kallisto_dir Directory containing Kallisto output folders for each sample.
#'
#' @return List containing gene and transcript count matrices and tximport result objects.
#' 
#' @description Imports Kallisto transcript-level abundance estimates into R matrices at gene and transcript level using tximport.
#'
#' @details This function takes a vector of SRR IDs and the path to the Kallisto output directory containing abundance estimates (abundance.tsv files) for each sample.
#' It imports the transcript counts into R using tximport, summarizing into gene-level counts based on gene annotation data for the specified species.
#' 
#' The output is a list containing the gene counts matrix, transcript counts matrix, and the tximport result objects.
#'
#' @export
run_tximport_kallisto <- function(srr_id, species = c("human", "mouse", "rat"), kallisto_dir) {
  species <- match.arg(species, c("human", "mouse", "rat"))
  edb <- function(species) {
    if (species == "human") {
      GenomicFeatures::transcripts(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
        columns = c("tx_id", "gene_id", "gene_name"),
        return.type = "DataFrame"
      )
    } else if (species == "mouse") {
      GenomicFeatures::transcripts(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
        columns = c("tx_id", "gene_id", "gene_name"),
        return.type = "DataFrame"
      )
    } else if (species == "rat") {
      GenomicFeatures::transcripts(EnsDb.Rnorvegicus.v79::EnsDb.Rnorvegicus.v79,
        columns = c("tx_id", "gene_id", "gene_name"),
        return.type = "DataFrame"
      )
    } else {
      return(NULL)
    }
  }
  gene_ensembl <- function(species) {
    if (species == "human") {
      return(org.Hs.eg.db::org.Hs.eg.db)
    } else if (species == "mousee") {
      return(org.Mm.eg.db::org.Mm.eg.db)
    }
    else {
      return(NULL)
    }
  }
  assign("Tx.ensemble", edb(species))
  Tx.ensemble <- get("Tx.ensemble")
  tx2gene <- Tx.ensemble[, c(1, 2)]
  files <- file.path(paste0(kallisto_dir, "/", srr_id, "/abundance.tsv"))
  names(files) <- srr_id
  file.exists(files)

  cat("generating counts table\n")
  txi.t <- tximport::tximport(files,
    type = "kallisto", tx2gene = tx2gene,
    txOut = TRUE,
    dropInfReps = TRUE
  )
  txi.g <- tximport::summarizeToGene(txi.t, tx2gene,
    ignoreTxVersion = TRUE, ignoreAfterBar = TRUE
  )
  gene_counts <- txi.g$counts
  gene_counts[is.na(gene_counts)] <- 0
  colnames(gene_counts) <- srr_id
  transcript_counts <- txi.t$counts
  transcript_counts[is.na(transcript_counts)] <- 0
  colnames(transcript_counts) <- srr_id
  annot_genes <- AnnotationDbi::select(gene_ensembl(species),
    keys = rownames(gene_counts), columns = c(
      "SYMBOL", "SYMBOL",
      "GENENAME"
    ), keytype = "ENSEMBL"
  )
  annot_genes2 <- annot_genes[match(
    rownames(gene_counts),
    annot_genes[, 1]
  ), , drop = FALSE]
  gene_counts <- cbind(annot_genes2, gene_counts)
  counts <- list(
    gene_counts = gene_counts, transcript_counts = transcript_counts,
    tximport_gene_data = txi.g, tximport_transcript_data = txi.t
  )
  return(counts)
}


## -------------------------------------------------------------------------
## https://rdrr.io/github/anilchalisey/rseqR/src/R/trim_fastq.R
## -------------------------------------------------------------------------

#' An R-based wrapper for Trim Galore!
#'
#' @description Run the Trim Galore! tool
#'
#' @details This script runs the Trim Galore! tool and requires installation
#' of both Cutadapt and Trim Galore!  It is essential that Cutadapt is in the
#' executable path otherwise this tool will not work.
#'
#' @param fastq1 a character vector indicating the read files to be trimmed.
#' @param fastq2 (optional) a character vector indicating read files to be
#' trimmmed.  If specified, it is assumed the reads are paired, and this vector
#' MUST be in the same order as those listed in \code{fastq1}.  If \code{NULL}
#' then it is assumed the reads are single-end.
#' @param adapter1 a character string specifying the adapter sequence to be
#' trimmed. If not specified explicitly, Trim Galore will try to auto-detect
#' whether the Illumina universal, Nextera transposase or Illumina small RNA
#' adapter sequence was used. Also see \code{illumina}, \code{nextera} and
#' \code{small_rna} options. If no adapter can be detected within the first 1
#' million sequences of the first file specified Trim Galore defaults to
#' \code{illumina}.
#' @param adapter2 a character string specifying an optional adapter sequence to
#' be trimmed off read 2 of paired-end files. This option requires paired-end
#' reads.
#' @param illumina a logical specifying that the adapter sequence to be trimmed
#' is the first 13bp of the Illumina universal adapter AGATCGGAAGAGC instead of
#' the default auto-detection of adapter sequence.  Default: \code{FALSE}
#' @param nextera adapter sequence to be trimmed is the first 12bp of the
#' Nextera adapter CTGTCTCTTATA instead of the default auto-detection of adapter
#' sequence.
#' @param small_rna a logical specifying that the adapter sequence to be trimmed
#' is the first 12bp of the Illumina Small RNA 3' Adapter TGGAATTCTCGG instead
#' of the default auto-detection of adapter sequence.  Selecting to trim
#' smallRNA adapters will also lower the \code{length} value to 18bp. If the
#' smallRNA libraries are paired-end then \code{adapter2} will be set to the
#' Illumina small RNA 5' adapter automatically (GATCGTCGGACT) unless
#' \code{adapter2} had been defined explicitly.
#' @param minlength an integer value; reads that become shorter than this length
#' as a result of either quality or adapter trimming are discarded. A value of 0
#' effectively disables this behaviour.  Default: 20 bp.  For paired-end files,
#' both reads of a read-pair need to be longer than bp to be printed out to
#' validated paired-end files. If only one read became too short there is the
#' possibility of keeping such unpaired single-end reads (see
#' \code{retain_unpaired}). Default pair-cutoff: 20 bp.
#' @param minqual an integer value specifying the quality threshold below which
#' to trim low-quality ends from reads in addition to adapter removal. Default
#' Phred score: 20.
#' @param trimN a logical specifying whether to remove Ns from the end of reads.
#' @param retainUnpaired a logical.  If only one of the two paired-end reads
#' become too short, the longer read will be written to either .unpaired_1.fq or
#' .unpaired_2.fq output files. The length cutoff for unpaired single-end reads
#' is governed by the parameters \code{retain1length} and \code{retain2length}.
#' Default: ON.
#' @param retain1length an integer.  Unpaired single-end read length cutoff
#' needed for read 1 to be written to .unpaired_1.fq output file. These reads
#' may then be mapped in single-end mode. Default: 35 bp.
#' @param retain2length an integer.  Unpaired single-end read length cutoff
#' needed for read 2 to be written to .unpaired_1.fq output file. These reads
#' may then be mapped in single-end mode. Default: 35 bp
#' @param clipR1 an integer instructing Trim Galore to remove the specified
#' number of bp from the 5' end of read 1 (or single-end reads). This may be
#' useful if the qualities were very poor, or if there is some sort of unwanted
#' bias at the 5' end. Default: 0
#' @param clipR2 an integer instructing Trim Galore to remove the specified
#' number of bp from the 5' end of read 2 (paired-end reads only). This may be
#' useful if the qualities were very poor, or if there is some sort of unwanted
#' bias at the 5' end. Default: 0
#' @param clip3primeR1 an integer instructing Trim Galore to remove the
#' specified number of bp from the 3' end of read 1 (or single-end reads) AFTER
#' adapter/quality trimming has been performed. This may remove some unwanted
#' bias from the 3' end that is not directly related to adapter sequence or
#' basecall quality. Default: 0.
#' @param clip3primeR2 an integer instructing Trim Galore to remove the
#' specified number of bp from the 3' end of read 1 (or single-end reads) AFTER
#' adapter/quality trimming has been performed. This may remove some unwanted
#' bias from the 3' end that is not directly related to adapter sequence or
#' basecall quality. Default: 0.
#' @param robust_check a logical indicating whether to check that the paired
#' files specified are matching and have equal numbers of reads.  Default:
#' \code{FALSE}
#' @param trimgalore a character string specifying the path to the trimgalore executable.
#' On Unix systems, if the executable is in \code{$PATH}, then it may be left as
#' the default. If it is not in \code{$PATH}, then the absolute path should be given.
# ` If using the WSL on Windows 10, then the path must be the absolute path in WSL,
#' unless the system has been set up as described in the vignette.
#' @param dest.dir a character string specifying the output directory.  If NULL
#' a directory named "TRIMMED_FASTQC" is created in the current working directory
#' [DEFAULT = NULL].
#' @param threads an integer value indicating the number of parallel threads to
#' be used by FastQC. [DEFAULT = maximum number of available threads - 1].
#'
#'
#' @importFrom parallel makeCluster stopCluster
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#'
#' @export
trimgalore_fastq <- function(fastq1, fastq2 = NULL, adapter1 = NULL, adapter2 = NULL,
                             illumina = FALSE, nextera = FALSE, small_rna = FALSE,
                             minlength = 20, minqual = 20, trimN = TRUE,
                             retainUnpaired = TRUE, retain1length = 35,
                             retain2length = 35, clipR1 = NULL, clipR2 = NULL,
                             clip3primeR1 = NULL, clip3primeR2 = NULL,
                             robust_check = FALSE, dest.dir = NULL,
                             threads = NULL, do.fastqc = FALSE,
                             trimgalore = "trim_galore") {
  doParallel::registerDoParallel(cores = 2)


  cmd <- trimgalore

  if (is.null(fastq2)) {
    paired <- FALSE
  } else {
    if (length(fastq1) != length(fastq2)) {
      stop("The number of forward and reverse reads do not match")
    }
    if (robust_check) {
      fq1lengths <- lapply(fastq1, function(x) {
        sprintf("gzip -cd %s | wc -l", x)
      })
      fq2lengths <- lapply(fastq1, function(x) {
        sprintf("gzip -cd %s | wc -l", x)
      })
      fq1lengths <- unlist(lapply(fq1lengths, run_cmd, intern = TRUE))
      fq2lengths <- unlist(lapply(fq2lengths, run_cmd, intern = TRUE))
      if (!identical(fq1lengths, fq2lengths)) {
        stop(
          "One or more of the forward and reverse reads pairs have differing number of reads.\n",
          "Are you sure the two lists are in the correct paired order?"
        )
      }
    }
    paired <- TRUE
  }

  if (is.null(dest.dir)) dest.dir <- "TRIMMED_FASTQC"
  dir.create(dest.dir, showWarnings = FALSE)

  cmd <- paste(
    cmd,
    "-q", minqual, "--length", minlength, "-o", dest.dir
  )
  if (do.fastqc) cmd <- paste(cmd, "--fastqc")
  if (!is.null(adapter1)) cmd <- paste(cmd, "--adapter", adapter1)
  if (!is.null(adapter2)) cmd <- paste(cmd, "--adapter2", adapter2)
  if (illumina) cmd <- paste(cmd, "--illumina")
  if (nextera) cmd <- paste(cmd, "--nextera")
  if (small_rna) cmd <- paste(cmd, "--small_rna")
  if (trimN) cmd <- paste(cmd, "--trim-n")
  if (!is.null(clipR1)) cmd <- paste(cmd, "--clip_R1", clipR1)
  if (!is.null(clipR2)) cmd <- paste(cmd, "--clip_R2", clipR2)
  if (!is.null(clip3primeR1)) cmd <- paste(cmd, "--three_prime_clip_R1", clip3primeR1)
  if (!is.null(clip3primeR2)) cmd <- paste(cmd, "--three_prime_clip_R2", clip3primeR2)
  if (paired) {
    cmd <- paste(cmd, "--paired")
    if (retainUnpaired) {
      cmd <- paste(
        cmd, "--retain_unpaired",
        "-r1", retain1length,
        "-r2", retain2length
      )
    }
  }


  cat("\nRemoving adapters and performing quality trimming...\n\n")

  if (is.null(threads)) {
    threads <- parallel::detectCores() - 1
  }

  if (threads > 1) {
    cl <- parallel::makeCluster(threads)
    doParallel::registerDoParallel(cl)

    if (paired) {
      foreach::foreach(i = seq_along(fastq1)) %dopar% {
        tgcmd <- sprintf("%s %s %s", cmd, fastq1[[i]], fastq2[[i]])
        run_cmd <- function(cmd, intern = FALSE) {
          if (.Platform$OS.type != "windows") {
            system(command = cmd, intern = intern)
          } else {
            shell(cmd = shQuote(cmd), shell = "bash", intern = intern)
          }
        }
        run_cmd(tgcmd)
      }
    } else {
      foreach::foreach(i = seq_along(fastq1)) %dopar% {
        tgcmd <- sprintf("%s %s", cmd, fastq1[[i]])
        run_cmd <- function(cmd, intern = FALSE) {
          if (.Platform$OS.type != "windows") {
            system(command = cmd, intern = intern)
          } else {
            shell(cmd = shQuote(cmd), shell = "bash", intern = intern)
          }
        }
        run_cmd(tgcmd)
      }
    }
    parallel::stopCluster(cl)
  } else {
    if (paired) {
      for (i in seq_along(fastq1)) {
        tgcmd <- sprintf("%s %s %s", cmd, fastq1[[i]], fastq2[[i]])
        run_cmd(tgcmd)
      }
    } else {
      for (i in seq_along(fastq1)) {
        tgcmd <- sprintf("%s %s", cmd, fastq1[[i]])
        run_cmd(tgcmd)
      }
    }
  }

  if (paired) {
    trimmed.files <- list.files(path = dest.dir, pattern = "*val", full.names = TRUE)
    lapply(trimmed.files, function(x) {
      file.to <- sub(".fq", ".fastq", gsub("_val_[0-9]", "_trimmed", x))
      file.rename(from = x, to = file.to)
    })
    unpaired.files <- list.files(path = dest.dir, pattern = "*unpaired", full.names = TRUE)
    lapply(unpaired.files, function(x) {
      file.to <- sub(".fq", ".fastq", gsub("_unpaired_[0-9]", "_unpaired", x))
      file.rename(from = x, to = file.to)
    })
  } else {
    trimmed.files <- list.files(path = dest.dir, pattern = "*val", full.names = TRUE)
    lapply(trimmed.files, function(x) {
      file.to <- sub(".fq", ".fastq", gsub("_val_[0-9]", "_trimmed", x))
      file.rename(from = x, to = file.to)
    })
  }
}
