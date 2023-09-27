##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##
library(biomaRt)
library(data.table)
ensembl <- useMart("ensembl")
species <- listDatasets(ensembl)
species <- data.table(species)
species[, species_name := sub("\\s*\\(.*\\)", "", description)]
species[, species_name := sub(" genes", "", species_name)]

write.csv(species, "dev/SPECIES_DICTIONARY.csv", row.names = FALSE, quote = FALSE)

SPECIES_TABLE = species

usethis::use_data(SPECIES_TABLE, overwrite = TRUE)
