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
inputs_files <- list(
    samples = playbase::SAMPLES,
    counts = playbase::COUNTS,
    contrasts = playbase::CONTRASTS
)

# generate tsv and csv files for all input files

create_dir("data-test/filetype")
lapply(names(inputs_files), function(x) {
    write.table(inputs_files[[x]], file = paste0("data-test//filetype//",x, ".tsv"), sep = "\t", quote = FALSE, row.names = TRUE)
    write.table(inputs_files[[x]], file = paste0("data-test//filetype//",x, ".csv"), sep = ",", quote = FALSE, row.names = TRUE)
})

# generating files with different header names
# special characters ok? How about spaces? Can they be? Which are OK, or must be avoided?



# generating samples columns with continuous and discrete values

# generating sample values with special characters

# check first column input that has not name



# samples.csv file:
# - Must it be called samples.csv? Is tab allowed?
# - Column headers - special characters ok? How about spaces? Can they be? Which are OK, or must be avoided?
# - Columns can they be continuous or must be discrete?
# - Sample values  - special characters ok? How about spaces? Can they be? Which are OK, or must be avoided?
# - Must the first column have no name? (please check!)

# counts.csv file:
# - Must it be called counts.csv? Is tab allowed?
# - Column headers - special characters ok? Can they be? Which are OK, or must be avoided?
# - Columns can they be fractional or must they be integer (raw) count?

# contrasts.csv file:
# - Must it be called counts.csv? Is tab allowed?
# - Column headers - special characters ok? How about spaces? Can they be? Which are OK, or must be avoided?
# - Column header special format? need _vs_ ? What is numerator / denominator?
# - Sample values  - special characters ok? How about spaces? Can they be? Which are OK, or must be avoided?
# - Must it be -1/0/+1 or can we use condition names? 