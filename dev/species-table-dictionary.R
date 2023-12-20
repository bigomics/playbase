##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##
library(biomaRt)
library(data.table)

# get species from main ensembl vertebrate

ensembl <- useMart("ensembl")
species <- listDatasets(ensembl)
species <- data.table(species)
species[, species_name := sub("\\s*\\(.*\\)", "", description)]
species[, species_name := sub(" genes", "", species_name)]

# add col ds to species table
species$mart <- "ensembl"


# get other taxonomies from ensenbl genome

genomes <- listEnsemblGenomes()

ds_genomes <- lapply(genomes$biomart, function(x) {
    #x = genomes$biomart[1]
    ensembl <- useEnsemblGenomes(biomart = x)
    ds <- listDatasets(ensembl)
    ds$mart <- x
    return(ds)
})

ds_genomes <- do.call("rbind", ds_genomes)

ds_genomes$species_name <- sub("\\s*\\(.*", "", ds_genomes$description)
# remove " genes" fro the string
ds_genomes$species_name <- sub(" genes", "", ds_genomes$species_name)
# remoe Short Variants from the string
ds_genomes$species_name <- sub(" Short Variants", "", ds_genomes$species_name)
# remove Structural Variants from the string
ds_genomes$species_name <- sub(" Structural Variants", "", ds_genomes$species_name)
# if duplicate, keep only first entry
ds_genomes <- ds_genomes[!duplicated(ds_genomes$species_name), ]

# reorder rows of ds_genomes to match species
ds_genomes <- ds_genomes[,colnames(species) ]

# merge two tables

species <- rbind(species, ds_genomes)

# order table by Human, Mouse and Rat to appear first in species_name
species <- species[order(species_name %in% c("Human", "Mouse", "Rat"), -as.character(species_name))]
    
# reverse order of table, where lst row becomes first and so on
species <- species[rev(seq_len(nrow(species)))]


# add servers
 species[, host :=  data.table::fcase(mart ==  "ensembl", "https://www.ensembl.org",
    mart == "plants_mart", "https://plants.ensembl.org",
    mart == "protists_mart", "https://protists.ensembl.org",
    mart == "metazoa_mart", "https://metazoa.ensembl.org",
    mart == "fungi_mart", "https://fungi.ensembl.org")]

# save Rdata

SPECIES_TABLE = species

# add "No organism" as last row, for all columns
SPECIES_TABLE <- rbind(SPECIES_TABLE, data.frame(dataset = "No organism", description = "No organism", version = "No organism", species_name = "No organism", mart = "No organism", host = "ensembl"))

usethis::use_data(SPECIES_TABLE, overwrite = TRUE)

# save SPECIES_DICTIONARY as tsv
write.table(SPECIES_TABLE, file = "dev/SPECIES_TABLE.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
