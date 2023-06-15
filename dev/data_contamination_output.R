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

output <- data.frame(filepaths = list.dirs("data-test/", recursive = TRUE))

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

unlist(lapply(output$filepaths,function(x){
    any(grepl("\\..+", dir(x), ignore.case = TRUE))
})) -> directories_with_files

output$directories_with_files <- directories_with_files
# same upload as opg

data_files <- list()
for (i in 1:nrow(output)) {
  #i = 11
  x <- output[i,]

  if (x[2] == TRUE) {
    samples_file <- list.files(x$filepaths, pattern = "ampl", full.names = TRUE, ignore.case = TRUE)
    counts_file <- list.files(x$filepaths, pattern = "unt", full.names = TRUE, ignore.case = TRUE)
    contrasts_file <- list.files(x$filepaths, pattern = "ntra", full.names = TRUE, ignore.case = TRUE)

    SAMPLES <- playbase::read.as_matrix(samples_file)
    COUNTS <- playbase::read.as_matrix(counts_file)
    CONTRASTS <- playbase::read.as_matrix(contrasts_file)

    # append the list of files to data_files
    data_files[[i]] <- list(SAMPLES = SAMPLES, COUNTS = COUNTS, CONTRASTS = CONTRASTS)
  }
}
# run pgx create

lapply(data_files, function(x){
    #x = data_files[[11]]
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

  if (output$directories_with_files[x] == TRUE) {
    if (!is.null(pgx_test$error)) {
        return(data.frame(filename = output$filepaths[[x]], playbase_error = as.character(pgx_test$error), stringsAsFactors = FALSE))
    } else {
        return(data.frame(filename = output$filepaths[[x]], playbase_error = "No error", stringsAsFactors = FALSE))
    }
  }
  
  if (output$directories_with_files[x] == FALSE) {
    return(data.frame(filename = output$filepaths[[x]], playbase_error = "Error: No file available", stringsAsFactors = FALSE))
  }
})

# save table as csv

error_table <- do.call(rbind, error_table)

write.csv(error_table, "dev//error_table.csv", row.names = FALSE)

# save pgx files

create_dir("dev//pgx")

lapply(1:length(pgx_files), function(x){
  #x = 171
  pgx_test <- pgx_files[[x]]
  #filename <- stringr::str_extract(output$filepaths[[x]], "(?<=/)[^/]*$")

  filename <- gsub("/", "_", output$filepaths[[x]])
  
  if (is.null(pgx_test$error) & output$directories_with_files[x] == TRUE) {
    pgx_test$result$name <- "empty"
    pgx_test$result$description <- "empty"
    pgx_test$result$datatype <- "empty"
    pgx_test$result$creator <- "empty"
    pgx_test$result$date <- "empty"

    playbase::pgx.save(pgx_test$result, file = paste0("dev//pgx//",filename,".pgx"))
  }
})



