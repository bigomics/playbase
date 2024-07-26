##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##


# Get species from annotation hub/orgdb
library(data.table)

ah <- AnnotationHub::AnnotationHub() ## make global??
db <- AnnotationHub::query(ah, "OrgDb")

M1 <- AnnotationHub::mcols(db)
M1 <- M1[ M1$rdataclass == "OrgDb",]
M2 <- orthogene::map_species(method="gprofiler")
M1 <- as.data.frame(M1)
M2 <- as.data.frame(M2)

both <- intersect( M1$taxonomyid, M2$taxonomy_id )
length(both)
ah.only <- setdiff( M1$taxonomyid, M2$taxonomy_id )
gp.only <- setdiff( M2$taxonomy_id, M1$taxonomyid  )

M1.only <- M1[match(ah.only,M1$taxonomyid),]
M2.only <- M2[match(gp.only,M2$taxonomy_id),]
sort(M1.only[,"species"]) ## only in OrgDB
sort(M2.only[,"scientific_name"])   ## only in gProfiler

length(both)
M1x <- M1[match(both, M1$taxonomyid), c("taxonomyid","title","species")]
M2x <- M2[match(both, M2$taxonomy_id), c("scientific_name","id")]
species <- cbind( M1x, M2x )
colnames(species) <- c("taxonomy_id","orgdb","orgdb.species","orthogene.scientific_name","gprofiler")
rownames(species) <- NULL
head(species)
dim(species)

simple_name <- species$orthogene.scientific_name

# rename Homo sapiens to Human in the species_name column
simple_name <- sub("Homo sapiens","Human", simple_name)
simple_name <- sub("Mus musculus","Mouse", simple_name)
simple_name <- sub("Rattus norvegicus","Rat", simple_name)
simple_name <- sub("lupus familiaris","lupus_familiaris",simple_name)
simple_name <- sapply(strsplit(simple_name,split=" "),
                      function(s) paste(head(s,2),collapse=" "))
species$species_name <- simple_name

# Create a new row with "No organism" in all columns
species <- rbind( "No organism", species )
species[1,1] <- "0000"
dim(species)

# Order table by Human, Mouse and Rat to appear first in species_name
preferred_order <- unique(c("Human", "Mouse", "Rat", "No organism", sort(species$species_name)))
species <- species[ match(preferred_order, species$species_name), ]

species <- species[,c("taxonomy_id","species_name","orgdb","orgdb.species","orthogene.scientific_name","gprofiler")]
dim(species)
rownames(species) <- NULL

write.table(species, file = "SPECIES_TABLE2.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

SPECIES_TABLE2 <- species
usethis::use_data(SPECIES_TABLE2, overwrite = TRUE)

# # Get species from main ensembl vertebrate
# library(biomaRt)
# library(data.table)


# ensembl <- useMart("ensembl")
# species <- listDatasets(ensembl)
# species <- data.table(species)
# species[, species_name := sub("\\s*\\(.*\\)", "", description)]
# species[, species_name := sub(" genes", "", species_name)]
# species$mart <- "ensembl"

# # Get other taxonomies from ensenbl genome
# genomes <- listEnsemblGenomes()
# ds_genomes <- lapply(genomes$biomart, function(x) {
#     ensembl <- useEnsemblGenomes(biomart = x)
#     ds <- listDatasets(ensembl)
#     ds$mart <- x
#     ds <- data.table(ds)
#     return(ds)
# })
# ds_genomes <- rbindlist(ds_genomes)

# # Clean up species names (remove genes, short variants, structural variants, etc)
# ds_genomes[ ,species_name := sub("\\s*\\(.*", "", description)]
# ds_genomes[ ,species_name := sub(" genes", "", species_name)]
# ds_genomes[ ,species_name := sub(" Short Variants", "", species_name)]
# ds_genomes[ ,species_name := sub(" Structural Variants", "", species_name)]

# # Add No organism, use I() because species and ds_genomes have "AsIs" class
# no_org <- data.table(dataset = I("No organism"),
#                      description = I("No organism"),
#                      version = I("No organism"),
#                      species_name = I("No organism"),
#                      mart = "ensembl")

# # Deduplicate species and Merge
# ds_genomes <- ds_genomes[!duplicated(species_name), ]
# species <- rbindlist(list(species, ds_genomes, no_org), use.names = TRUE)

# # add servers
# species[, host :=  data.table::fcase(mart ==  "ensembl", "https://www.ensembl.org",
#     mart == "plants_mart", "https://plants.ensembl.org",
#     mart == "protists_mart", "https://protists.ensembl.org",
#     mart == "metazoa_mart", "https://metazoa.ensembl.org",
#     mart == "fungi_mart", "https://fungi.ensembl.org",
#     TRUE, "")]

# # Order table by Human, Mouse and Rat to appear first in species_name
# preferred_order <- c("Human", "Mouse", "Rat", "No organism")
# species[, species_name := factor(species_name, levels = c(preferred_order, sort(setdiff(unique(species_name), preferred_order))))]
# setorder(species, mart, species_name)
# species

# # Save as Rdata and tsv
# SPECIES_TABLE <- species
# usethis::use_data(SPECIES_TABLE, overwrite = TRUE)
# write.table(SPECIES_TABLE, file = "dev/SPECIES_TABLE.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
