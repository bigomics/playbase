##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2023 BigOmics Analytics SA. All rights reserved.
##

# read PGX_CHECKS.csv

PGX_CHECKS <- read.csv("dev/PGX_CHECKS.csv", header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, fileEncoding="UTF-8-BOM")

# remove placeholder rows for future checks
PGX_CHECKS <- PGX_CHECKS[!PGX_CHECKS$message == "",]

usethis::use_data(PGX_CHECKS, overwrite = TRUE)
