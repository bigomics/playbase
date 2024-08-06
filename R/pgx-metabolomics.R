# convert IDs to CHEBI using base R functions

#' @export
convert_probe_to_chebi <- function(ids, probe_type) {
    probe_type_dictionary <- playdata::METABOLITE_ANNOTATION
    valid_probe_types <- colnames(probe_type_dictionary)

    # for id that are "", set it to na
    ids[ids == ""] <- NA

    # check that probetype is valid
    probe_type <- match.arg(probe_type, valid_probe_types)

    if (probe_type == "ChEBI") {
        # keep only numbers in ids, as chebi ids are numeric
        chebi_ids <- gsub("[^0-9]", "", ids)
        # check that ChEBI ids are in the dictionary
        chebi_ids <- ifelse(chebi_ids %in% probe_type_dictionary$ChEBI, chebi_ids, NA)
    } else if (probe_type == "HMDB") {
        # use match
        matches <- match(ids, probe_type_dictionary$HMDB)
        chebi_ids <- probe_type_dictionary[matches, "ChEBI"]

        # return NA for unmatched IDs
        return(chebi_ids)
    } else if (probe_type == "KEGG") {
        # use match
        matches <- match(ids, probe_type_dictionary$KEGG)
        chebi_ids <- probe_type_dictionary[matches, "ChEBI"]
    } else {
        return(NA)
    }

    return(chebi_ids)
}

#' @export
pgx.addMetaboliteAnnotation <- function(counts, annotation_table) {
    return(NULL)
}
