##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

# function to create dir

create_dir <- function(directory){
    # Check if directory exists
    if (!dir.exists(directory)) {
    # Create the directory
    dir.create(directory, recursive = TRUE)
    }
}


# get template data
input_files <- list(
    samples = playbase::SAMPLES,
    counts = playbase::COUNTS[1:300,], # reduce size to speed up tests
    contrasts = playbase::CONTRASTS[,1:2] # reduce size to speed up tests
)

# special characters

characters <- list(
    digits = c(0,1,9),
    punctuation = c( ".", "-", "?", "!", ';', ":", ","),
    special = c("~", "@", "#", "$", "%", "_", "+", " "),
    escape = c("\\"),
    control = c("\n", "\t", "\r", "\b", "\a", "\f", "\v"),
    noncontrol = c("\u")
)

# generate tsv and csv files for all input files

create_dir("data-test/filetype")
lapply(names(inputs_files), function(x) {
    write.table(inputs_files[[x]], file = paste0("data-test//filetype//",x, ".tsv"), sep = "\t", quote = FALSE, row.names = TRUE)
    write.table(inputs_files[[x]], file = paste0("data-test//filetype//",x, ".csv"), sep = ",", quote = FALSE, row.names = TRUE)
})

# generating files with different header names
# special characters ok? How about spaces? Can they be? Which are OK, or must be avoided?



# generating samples, counts and contrast columns with continuous and discrete values as metadata

# generating sample values with special characters

# check first column input that has not name

# generate files with file names (count1.csv, Counts.csv...)

# generate counts with different gene names Must it be -1/0/+1 or can we use condition names? 

# check if `` fixes the char problem in R