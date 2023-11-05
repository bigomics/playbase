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
#' @param item Select which elements to return c("counts", "samples", "contrats", "all")
#' by default it returns a list with all the data elements to build a pgx object.
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
get_mini_example_data <- function(item = "all") {

  if (!item %in% c("counts", "samples", "contrats", "all")) {
    stop("Item not matched. Please use: 'counts', 'samples', 'contrats', or 'all'")
  }
  if (length(item) > 1L) {
    stop("Call only one item at a time. For multiple items call, use 'all'")
  }

  # Prepare input data
  counts <- playbase::COUNTS
  samples <- playbase::SAMPLES
  contrast <- playbase::CONTRASTS

  n_genes <- round(seq(1, nrow(counts), length.out = 500))

  # Subset each data frame to facilitate testing
  mini_counts <- counts[n_genes, c(1:3, 8:9, 11:12)]
  mini_samples <- samples[colnames(mini_counts), ]
  mini_contrast <- contrast[1:3, 1:3]
  mini_data <- list(counts = mini_counts,
                    samples = mini_samples,
                    contrast = mini_contrast)

  if (item == "all") {

    return(mini_data)

  } else {
    return(mini_data[[item]])
  }
}


#' Get row names for a data.table
#'
#' @description This function aims to extract the rownames of a data.table object as
#' if it would be a data.frame. The implementation is slow but it helps to maintain back
#' compatibility with data.frame objects.
#' 
#' @param dt A data.table object.
#'
#' @return Character vector of row names
#'
#' @examples
#' \dontrun{
#' dt <- data.table(rn = letters[1:5], x = 1:5)
#' rownames(dt)
#' }
#' @export
rownames.data.table <- function(dt) {
  
  if(is.character(dt[[1]])) {
    return(dt[[1]])
    
  } else if("rn" %in% names(dt)) {
    return(dt$rn)
    
  } else {
    warning("No column contains row names")
    return(NULL)
  }
}
