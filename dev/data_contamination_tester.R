##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

filepath <-dir(filepath, recursive = TRUE, full.names = TRUE)

# same upload as opg

samples <- playbase::read.as_matrix(fn2)
counts <- playbase::read.as_matrix(fn2)
contrasts <- playbase::read.as_matrix(fn2)



# pgx create
filepath <- "data-test/first_column_name/"

