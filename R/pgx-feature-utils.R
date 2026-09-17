## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
## It owns private feature-name parsing and downstream matrix preparation.
## Generic numerical operations delegate to playbase.preprocess.

# Extracts the text before the first feature-name colon.
# Unprefixed values are returned unchanged for downstream grouping.
# This is identifier parsing, not preprocessing dispatch.
.pgx_feature_prefix <- function(names) {
  sub(":.*", "", names)
}

# Imputes a downstream matrix with the canonical SVD2 family operation.
# Fully prefixed multi-omics rows are processed independently by layer.
# This is analysis preparation and does not mutate PGX preprocessing state.
.pgx_impute_svd2 <- function(X) {
  playbase.preprocess::pp.impute(
    X,
    layers = playbase.preprocess::pp.inferLayers(X),
    method = "SVD2"
  )
}
