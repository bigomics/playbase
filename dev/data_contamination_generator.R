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

# code to add character in the middle of a string

add_character_in_second_position <- function(string, character) {
  new_string <- paste0(substring(string, 1, 1), character, substring(string, 2))
  return(new_string)
}


# get template data
input_files <- list(
    samples = playbase::SAMPLES,
    counts = playbase::COUNTS[1:300,], # reduce size to speed up tests
    contrasts = playbase::CONTRASTS # reduce size to speed up tests
)

# output_char files, each element on the list is a test
output_char <- list()

# special characters
create_dir("data-test/filetype")

characters <- list(
    digits = c(0,1,9),
    punctuation = c( ".", "-", "?", "!", ';', ":", ","),
    special = c("~", "@", "#", "$", "%", "_", "+", " "),
    escape = c("\\"),
    control = c("\\n", "\\t", "\\r", "\\b", "\\a", "\\f", "\\v"),
    noncontrol = c("\\u")
)

# generate tsv and csv files for all input files

lapply(names(inputs_files), function(x) {   
    write.table(inputs_files[[x]], file = paste0("data-test//filetype//",x, ".tsv"), sep = "\t", quote = FALSE, row.names = TRUE)
    write.table(inputs_files[[x]], file = paste0("data-test//filetype//",x, ".csv"), sep = ",", quote = FALSE, row.names = TRUE)
})

# generating files with different header names
# special characters ok? How about spaces? Can they be? Which are OK, or must be avoided?

# prepend header with special characters
lapply(unlist(characters), function (x){
    #x = unlist(characters)[1]
    sample <- input_files$samples
    contrast <- input_files$contrast
    count <- input_files$count

    colnames(sample) <- paste0(x, colnames(sample))

    colnames(contrast) <- paste0(x, colnames(contrast))

    return(list(sample = sample, contrast = contrast, count = count))

}) -> output_char$prepend_header

lapply(unlist(characters), function (x){
    return()
}) -> output_char$prepend_header

# append header with special characters

lapply(unlist(characters), function (x){
    #x = unlist(characters)[1]
    sample <- input_files$samples
    contrast <- input_files$contrast
    count <- input_files$count

    colnames(sample) <- paste0(colnames(sample),x)

    colnames(contrast) <- paste0(colnames(contrast),x)

    return(list(sample = sample, contrast = contrast, count = count))

}) -> output_char$append_header

lapply(unlist(characters), function (x){
    #x = unlist(characters)[1]
    sample <- input_files$samples
    contrast <- input_files$contrast
    count <- input_files$count

    colnames(sample) <- add_character_in_second_position(colnames(sample),x)

    colnames(contrast) <- add_character_in_second_position(colnames(contrast),x)

    return(list(sample = sample, contrast = contrast, count = count))

}) -> output_char$middle_header

# add special characters to sample metadata AND contrast file


# TODO
# TODO
# TODO



# generating sample names with special characters

lapply(unlist(characters), function (x){
    #x = unlist(characters)[1]
    sample <- input_files$samples
    contrast <- input_files$contrast
    count <- input_files$count

    rownames(sample) <- paste0(x, rownames(sample))

    colnames(count) <-  rownames(sample)
    
    return(list(sample = sample, contrast = contrast, count = count))

}) -> output_char$sample_names

# generating samples, counts and contrast columns with continuous and discrete values as metadata

samples <- input_files$samples

samples$continuous <- rnorm(nrow(samples))

output_metadata_continuous <- list(
    samples = samples,
    counts = input_files$counts,
    contrasts = input_files$contrasts)

# check first column input that has not name

# generate files with file names (count1.csv, Counts.csv...)

# generate counts with different gene names Must it be -1/0/+1 or can we use condition names? 

# check if `` fixes the char problem in R