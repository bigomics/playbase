# convert IDs to CHEBI using base R functions

#' @export
convert_probe_to_chebi <- function(ids, probe_type) {
    probe_type_dictionary <- playdata::METABOLITE_ANNOTATION
    valid_probe_types <- colnames(probe_type_dictionary)

    # check that probetype is valid
    probe_type <- match.arg(probe_type, valid_probe_types)

    if (probe_type == "ChEBI") {
        return(ids)
    } else if (probe_type == "HMDB") {
        # use match
        matches <- match(ids, probe_type_dictionary$HMDB)
        chebi_ids <- probe_type_dictionary$CHEBI[matches]

        # return NA for unmatched IDs
        return(ifelse(is.na(chebi_ids), NA, chebi_ids))
    } else if (probe_type == "KEGG") {
        matches <- match(ids, probe_type_dictionary$KEGG)
        chebi_ids <- probe_type_dictionary$CHEBI[matches]

        # return NA for unmatched IDs
        return(ifelse(is.na(chebi_ids), NA, chebi_ids))
    } else {
        return(NA)
    }
    return(chebi_ids)
}

#' @export
pgx.addMetaboliteAnnotation <- function(counts, annotation_table) {
    return(NULL)
}
