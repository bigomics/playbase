##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

filepaths <- ?list.dirs("data-test/", recursive = TRUE)

file.path(filepath[1])

lapply(filepaths,function(x){
    any(grepl("\\..+", dir(x), ignore.case = TRUE))
}) -> directories_with_files

filepath <- filepaths[unlist(directories_with_files)]

# same upload as opg

lapply(filepath, function(x){
    print(x)
    samples_file <- list.files(x, pattern = "ampl", full.names = TRUE, ignore.case = TRUE)
    counts_file <- list.files(x, pattern = "unt", full.names = TRUE, ignore.case = TRUE)
    contrasts_file <- list.files(x, pattern = "ntra", full.names = TRUE, ignore.case = TRUE)

    samples <- playbase::read.as_matrix(samples_file)
    counts <- playbase::read.as_matrix(counts_file)
    contrasts <- playbase::read.as_matrix(contrasts_file)
}) -> data_files

# run pgx create

lapply(data_files, function(x){
    pgx <- playbase::pgx.createPGX(
        counts = playbase::COUNTS,
        samples = playbase::SAMPLES,
        contrasts = playbase::CONTRASTS
        )

})

