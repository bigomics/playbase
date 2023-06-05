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



samples <- playbase::read.as_matrix(fn2)
counts <- playbase::read.as_matrix(fn2)
contrasts <- playbase::read.as_matrix(fn2)



# pgx create
filepath <- "data-test/first_column_name/"

