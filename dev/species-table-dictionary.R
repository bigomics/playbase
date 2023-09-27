##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


SPECIES_TABLE <- read.csv("dev/SPECIES_DICTIONARY.csv", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding="UTF-8-BOM")

usethis::use_data(SPECIES_TABLE, overwrite = TRUE)
