mem.proc <- function(digits = 0) {
  mem <- "[? MB]"
  if (Sys.info()["sysname"] %in% c("Linux")) {
    file <- paste("/proc", Sys.getpid(), "stat", sep = "/")
    what <- vector("character", 52)
    ## In your logging routine
    vsz <- as.numeric(scan(file, what = what, quiet = TRUE)[23])
    vsz <- vsz / (1024**2) ## MB
    mem <- paste0(round(vsz, digits), "MB")
  }
  mem
}

info <- function(..., type = "INFO") {
  dd <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg <- "some message"
  msg <- sapply(list(...), paste, collapse = " ")
  dd <- paste0("[", dd, "]")
  mm <- paste0("[", mem.proc(), "]")
  type <- paste0("[", type, "]")
  message(paste0(type, dd, mm, " --- ", sub("\n$", "", paste(msg, collapse = " "))))
}

dbg <- function(...) info(..., type = "DBUG")


#' Get mini example dataset
#'
#' @description
#' Returns a mini example dataset with sample information,
#' gene expression counts, and contrasts.
#'
#' @return List containing:
#' \itemize{
#'   \item samples - data.frame with sample info
#'   \item counts - matrix of gene counts
#'   \item contrast - data.frame with contrasts
#' }
#'
#' @examples
#' d <- get_mini_example_data()
#' @export
get_mini_example_data <- function() {
  counts <- playbase::COUNTS
  samples <- playbase::SAMPLES
  contrast <- playbase::CONTRASTS

  n_genes <- round(seq(1, nrow(counts), length.out = 100))

  # Subset each data frame to facilitate testing
  mini_counts <- counts[n_genes, c(1:3, 8:9, 11:12)]
  mini_samples <- samples[colnames(mini_counts), ]
  mini_contrast <- contrast[1:3, 1:3]

  mini_data <- list(counts = mini_counts, samples = mini_samples, contrast = mini_contrast)

  return(mini_data)
}
