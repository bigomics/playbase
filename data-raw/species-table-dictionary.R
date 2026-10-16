##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##

# Get species from annotation hub/orgdb
library(data.table)
setwd("..")
source("dev/include.R", chdir=TRUE)

.buildSpeciesTable <- function() {
  
  species <- data.frame(.getSpeciesTable.ANNOTHUB())
  colnames(species) <- sub("^species$","ah_species",colnames(species))
  
  # rename Homo sapiens to Human in the species_name column
  species$species_name <- species$ah_species
  species$species_name <- sub("Homo sapiens","Human",species$species_name)
  species$species_name <- sub("Mus musculus","Mouse",species$species_name)
  species$species_name <- sub("Rattus norvegicus","Rat",species$species_name)
  
  # remove duplicates, use species_name as rownames
  species <- species[!duplicated(species$species_name),]
  rownames(species) <- species$species_name
  
  # add plasmodium by hand (we add the orgdb manually)
  if(1) {
    pf <- as.character(species["Human",])
    names(pf) <- colnames(species)
    pf <- gsub("Homo sapiens|Human","Plasmodium falciparum",pf)
    pf <- sub("org.Hs.eg.db.sqlite","org.Pf.plasmo.db",pf)
    pf["ah_id"] <- NA
    pf["taxonomyid"] <- "36329"
    species <- rbind(species, "Plasmodium falciparum" = pf)
  }
  
  # union of species supported either by annothub or gprofiler
  M <- .getSpeciesTable.GPROFILER()
  all_id <- union(species$taxonomyid, M$taxonomy_id)
  species <- species[match(all_id,species$taxonomyid),]
  M <- M[match(all_id,M$taxonomy_id),]
  species$taxonomyid <- all_id
  
  # add orthogene species name for lookup. sometimes different
  species$gprofiler_species <- M$scientific_name
  species$gprofiler_id  <- M$id
  species$display_name  <- ifelse(!is.na(species$species_name), species$species_name, M$display_name)
  species$species_name  <- ifelse(!is.na(species$species_name), species$species_name, M$scientific_name)
  
  # Create a new row with "No organism" in all columns
  species <- rbind(species, "No organism" = NA)
  species['No organism',c('ah_species','species_name','display_name')] <- 'No organism'
  
  col_order <- unique(c('species_name','display_name',sort(colnames(species))))
  species <- species[,col_order]

  # Order table by Human, Mouse and Rat to appear first in species_name
  top.orgs <- c("Human", "Mouse", "Rat", "No organism")
  preferred_order <- unique(c(match(top.orgs, species$species_name),
    order(species$species_name)))
    species <- species[ preferred_order,]
  rownames(species) <- NULL
  

  species <- data.frame(species)
  return(species)
}

SPECIES_TABLE <- .buildSpeciesTable()
dim(SPECIES_TABLE)
head(SPECIES_TABLE)

ens <- read.csv("data-raw/ensembl_prefix.txt", sep="\t", header=FALSE)
uni <- read.csv("data-raw/uniprot-speclist.csv", sep=",", header=TRUE)
head(ens)
head(uni)
ens.scientific_name <- sub(" \\(.*","",ens$V2)
uni.scientific_name <- sub(" \\(.*","",uni$scientific_name)

## fill 5-letter uniprot species mnemonic. If not exists, create 3+2
## defaults.
jj <- match( SPECIES_TABLE$taxonomyid, uni$Taxon_ID )
SPECIES_TABLE$uniprot_species <- uni$Code[jj]
ii <- which(is.na(SPECIES_TABLE$uniprot_species))
mn5 <- lapply(strsplit(SPECIES_TABLE$ah_species,split=" "),substring,1,3)
mn5 <- toupper(substring(sapply(mn5,paste0,collapse=""),1,5))
SPECIES_TABLE$uniprot_species[ii] <- mn5[ii]
sum(duplicated(SPECIES_TABLE$uniprot_species))  ## ??

## fill with ENSEMBL mnemonic
jj <- match( SPECIES_TABLE$ah_species, ens.scientific_name )
ii <- which(is.na(jj))
jj[ii] <- match( SPECIES_TABLE$gprofiler_species[ii], ens.scientific_name )
SPECIES_TABLE$ensembl_species <- ens$V1[jj]

write.table(SPECIES_TABLE, file = "data-raw/SPECIES_TABLE.tsv", sep = "\t",
  quote = FALSE, row.names = FALSE)

usethis::use_data(SPECIES_TABLE, overwrite = TRUE)

