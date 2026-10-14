##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##


## ================================================================================
## ====================== GENE INFO FUNCTIONS======================================
## ================================================================================


#' Get annotation info for one feature of a PGX object
#'
#' Forwards to `playbase.annotation::pgx.getFeatureInfo()`.
#'
#' @export
pgx.getFeatureInfo <- function(pgx, feature) {
  playbase.annotation::pgx.getFeatureInfo(pgx = pgx, feature = feature)
}

