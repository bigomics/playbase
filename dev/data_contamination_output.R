##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

create_dir <- function(directory){
    # Check if directory exists
    if (!dir.exists(directory)) {
    # Create the directory
    dir.create(directory, recursive = TRUE)
    }
}

createPGX_safe <- purrr::safely(function(x) {
  pgx <- playbase::pgx.createPGX(
    counts = x$COUNTS,
    samples = x$SAMPLES,
    contrasts = x$CONTRASTS
  )
  pgx <- playbase::pgx.computePGX(
    pgx = pgx,
    gx.methods = c("ttest.welch", "trend.limma", "edger.qlf")[1],
    gset.methods = c("fisher", "gsva", "fgsea")[1]
    )
  pgx
})

filepaths <- list.dirs("data-test/", recursive = TRUE)

lapply(filepaths,function(x){
    any(grepl("\\..+", dir(x), ignore.case = TRUE))
}) -> directories_with_files

filepaths <- filepaths[unlist(directories_with_files)]

# same upload as opg

lapply(filepaths, function(x){
    print(x)
    #x = filepaths[[4]]
    samples_file <- list.files(x, pattern = "ampl", full.names = TRUE, ignore.case = TRUE)
    counts_file <- list.files(x, pattern = "unt", full.names = TRUE, ignore.case = TRUE)
    contrasts_file <- list.files(x, pattern = "ntra", full.names = TRUE, ignore.case = TRUE)

    SAMPLES <- playbase::read.as_matrix(samples_file)
    COUNTS <- playbase::read.as_matrix(counts_file)
    CONTRASTS <- playbase::read.as_matrix(contrasts_file)

    return(list(SAMPLES = SAMPLES, COUNTS = COUNTS, CONTRASTS = CONTRASTS))

}) -> data_files

# run pgx create

lapply(data_files, function(x){
    #x = data_files[[10]]
    res <- tryCatch(
    createPGX_safe(x),
    error = function(e) {
      e
    }
    )
    res
}) -> pgx_files

# create a table with results from error and the filename associated with it.


error_table <- lapply(1:length(pgx_files), function(x){
  #x = 10
  pgx_test <- pgx_files[[x]]
  if (!is.null(pgx_test$error)) {
      return(data.frame(filename = filepaths[[x]], playbase_error = as.character(pgx_test$error), stringsAsFactors = FALSE))
  } else {
      return(data.frame(filename = filepaths[[x]], playbase_error = "No error", stringsAsFactors = FALSE))
  }
})

# save table as csv

error_table <- do.call(rbind, error_table)

write.csv(error_table, "dev//error_table.csv", row.names = FALSE)

# save pgx files

create_dir("dev//pgx")

lapply(1:length(pgx_files), function(x){
  #x = 2
  pgx_test <- pgx_files[[x]]
  filename <- stringr::str_extract(filepaths[[x]], "(?<=/)[^/]*$")
  

  if (is.null(pgx_test$error)) {
    pgx_test$result$name <- "empty"
    pgx_test$result$description <- "empty"
    pgx_test$result$datatype <- "empty"
    pgx_test$result$creator <- "empty"
    pgx_test$result$date <- "empty"

    playbase::pgx.save(pgx_test$result, file = paste0("dev//pgx//",filename,".pgx"))
  }
})

