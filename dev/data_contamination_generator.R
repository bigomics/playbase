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

characters <- list(
    digits = c(0,1,9),
    punctuation = c( ".", "-", "?", "!", ';', ":", ","),
    special = c("~", "@", "#", "$", "%", "_", "+", " "),
    quotation = c("'", '"'),
    escape = c("\\"),
    control = c("\\n", "\\t", "\\r", "\\b", "\\a", "\\f", "\\v"),
    noncontrol = c("\\u")
)

# generate tsv and csv files for all input files

create_dir("data-test/first_column_name")

lapply(names(input_files), function(x) {

    #x = input_files[[1]]

    matrix_colname <- input_files[[x]]

    rowname_df <- data.frame(row_name = rownames(matrix_colname))

    matrix_colname <- cbind(rowname_df, matrix_colname)

    write.table(matrix_colname, file = paste0("data-test//first_column_name//",x, ".csv"), sep = ",", quote = FALSE, row.names = FALSE)
})

# generate files with file names (count1.csv, Counts.csv...)

create_dir("data-test/filenames")

lapply(names(input_files), function(x) {

    # filename with first capital letter
    write.table(matrix_colname, file = paste0("data-test//filenames//",stringr::str_to_title(x), ".csv"), sep = ",", quote = FALSE, row.names = FALSE)
    # filename in all capital letters
    write.table(matrix_colname, file = paste0("data-test//filenames//",stringr::str_to_upper(x),".csv"), sep = ",", quote = FALSE, row.names = FALSE)
    # prepent before filaname
    write.table(matrix_colname, file = file.path("data-test//filenames//",paste0("lalalala01_",x,".csv")), sep = ",", quote = FALSE, row.names = FALSE)
    # append after filename
    write.table(matrix_colname, file = file.path("data-test//filenames//",paste0(x,"_lalalala01_",".csv")), sep = ",", quote = FALSE, row.names = FALSE)

})






# generate tsv and csv files for all input files

create_dir("data-test/filetype")

lapply(names(input_files), function(x) {
    write.table(input_files[[x]], file = paste0("data-test//filetype//",x, ".tsv"), sep = "\t", quote = FALSE, row.names = TRUE)
    write.table(input_files[[x]], file = paste0("data-test//filetype//",x, ".csv"), sep = ",", quote = FALSE, row.names = TRUE)
    write.table(input_files[[x]], file = paste0("data-test//filetype//",x, ".txt"), quote = FALSE, row.names = TRUE)
})

# generate counts with different gene names Must it be -1/0/+1 or can we use condition names? 


create_dir("data-test/contrastinputs")

lapply(names(input_files), function(x) {
    x = input_files[["contrast"]]

    if(x == "contrasts"){
        # prepara long contrast (-1/1)



        c1 <- 

        # prepare group-based short contrast (-1/1)

        c2 <- 
        # prepare contrast with variable names (phenotypes) in cells

    }else{
        write.table(input_files[[x]], file = paste0("data-test//filetype//",x, ".csv"), sep = ",", quote = FALSE, row.names = TRUE)
    }

    
})














# generating files with different header names

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

# add special characters to sample metadata AND respective contrast

lapply(unlist(characters), function (x){
    #x = unlist(characters)[1]
    sample <- input_files$samples
    contrast <- input_files$contrast
    count <- input_files$count

    sample[,1] <- paste0(x, sample[,1])

    rownames(contrast) <- paste0(x, rownames(contrast))

    contrast_split <- strsplit(colnames(contrast), split = "_")

    contrast_contaminated <-lapply(contrast_split, function(idx){
        idx = contrast_split[[1]]
        idx[3] <- paste0(x, idx[3])
        paste(idx, collapse = "_")
    })

    colnames(contrast) <- contrast_contaminated    
    
    return(list(sample = sample, contrast = contrast, count = count))

}) -> output_char$sample_metadata


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

