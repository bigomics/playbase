##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##
library(biomaRt)
library(data.table)

# Get species from main ensembl vertebrate

ensembl <- useMart("ensembl")
species <- listDatasets(ensembl)
species <- data.table(species)
species[, species_name := sub("\\s*\\(.*\\)", "", description)]
species[, species_name := sub(" genes", "", species_name)]
species$mart <- "ensembl"

# Get other taxonomies from ensenbl genome
genomes <- listEnsemblGenomes()
ds_genomes <- lapply(genomes$biomart, function(x) {
    ensembl <- useEnsemblGenomes(biomart = x)
    ds <- listDatasets(ensembl)
    ds$mart <- x
    ds <- data.table(ds)
    return(ds)
})
ds_genomes <- rbindlist(ds_genomes)

# Clean up species names (remove genes, short variants, structural variants, etc)
ds_genomes[ ,species_name := sub("\\s*\\(.*", "", description)]
ds_genomes[ ,species_name := sub(" genes", "", species_name)]
ds_genomes[ ,species_name := sub(" Short Variants", "", species_name)]
ds_genomes[ ,species_name := sub(" Structural Variants", "", species_name)]

# Add No organism, use I() because species and ds_genomes have "AsIs" class 
no_org <- data.table(dataset = I("No organism"), 
                     description = I("No organism"), 
                     version = I("No organism"), 
                     species_name = I("No organism"), 
                     mart = "ensembl")

# Deduplicate species and Merge
ds_genomes <- ds_genomes[!duplicated(species_name), ]
species <- rbindlist(list(species, ds_genomes, no_org), use.names = TRUE)

# add servers
species[, host :=  data.table::fcase(mart ==  "ensembl", "https://www.ensembl.org",
    mart == "plants_mart", "https://plants.ensembl.org",
    mart == "protists_mart", "https://protists.ensembl.org",
    mart == "metazoa_mart", "https://metazoa.ensembl.org",
    mart == "fungi_mart", "https://fungi.ensembl.org",
    TRUE, "")]

# Order table by Human, Mouse and Rat to appear first in species_name
preferred_order <- c("Human", "Mouse", "Rat", "No organism")
species[, species_name := factor(species_name, levels = c(preferred_order, sort(setdiff(unique(species_name), preferred_order))))]
setorder(species, mart, species_name)
species

# Save as Rdata and tsv
SPECIES_TABLE <- species
usethis::use_data(SPECIES_TABLE, overwrite = TRUE)
write.table(SPECIES_TABLE, file = "dev/SPECIES_TABLE.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
