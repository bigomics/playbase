## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.

## Compatibility shim for the preprocessing layer, which now lives in the
## `playbase.preprocess` leaf package. Each name that package exports is
## re-exported here under its original playbase name, so that no caller --
## playbase's own consumer files or omicsplayground -- has to change.
##
## `playbase.preprocess` is a `Depends`, not an `Imports`, so its exports also
## sit on the search path once playbase is attached; that is what keeps
## playbase's unqualified internal calls (logCPM(), imputeMissing(), ...)
## resolving without a single edited call site.
##
## `pgx.preprocess` is the one exception: it lost its `samples`/`contrasts`
## arguments on the way out, so it gets a real wrapper below rather than a bare
## delegation.

#' @importFrom playbase.preprocess betaToM
#' @export
playbase.preprocess::betaToM

#' @importFrom playbase.preprocess detectOutlierSamples
#' @export
playbase.preprocess::detectOutlierSamples

#' @importFrom playbase.preprocess getPrior
#' @export
playbase.preprocess::getPrior

#' @importFrom playbase.preprocess imputeMissing
#' @export
playbase.preprocess::imputeMissing

#' @importFrom playbase.preprocess imputeMissing.mox
#' @export
playbase.preprocess::imputeMissing.mox

#' @importFrom playbase.preprocess is.multiomics
#' @export
playbase.preprocess::is.multiomics

#' @importFrom playbase.preprocess log1s
#' @export
playbase.preprocess::log1s

#' @importFrom playbase.preprocess logCPM
#' @export
playbase.preprocess::logCPM

#' @importFrom playbase.preprocess mToBeta
#' @export
playbase.preprocess::mToBeta

#' @importFrom playbase.preprocess maxMedianNormalization
#' @export
playbase.preprocess::maxMedianNormalization

#' @importFrom playbase.preprocess maxSumNormalization
#' @export
playbase.preprocess::maxSumNormalization

#' @importFrom playbase.preprocess mofa.get_prefix
#' @export
playbase.preprocess::mofa.get_prefix

#' @importFrom playbase.preprocess nmfImpute
#' @export
playbase.preprocess::nmfImpute

#' @importFrom playbase.preprocess nmfImpute2
#' @export
playbase.preprocess::nmfImpute2

#' @importFrom playbase.preprocess normalizeExpression
#' @export
playbase.preprocess::normalizeExpression

#' @importFrom playbase.preprocess normalizeMethylation
#' @export
playbase.preprocess::normalizeMethylation

#' @importFrom playbase.preprocess normalizeMultiOmics
#' @export
playbase.preprocess::normalizeMultiOmics

#' @importFrom playbase.preprocess normalizeRLE
#' @export
playbase.preprocess::normalizeRLE

#' @importFrom playbase.preprocess normalizeTMM
#' @export
playbase.preprocess::normalizeTMM

#' @importFrom playbase.preprocess perseusImpute
#' @export
playbase.preprocess::perseusImpute

#' @importFrom playbase.preprocess pgx.countNormalization
#' @export
playbase.preprocess::pgx.countNormalization

#' @importFrom playbase.preprocess plotOutlierScores
#' @export
playbase.preprocess::plotOutlierScores

#' @importFrom playbase.preprocess referenceNormalization
#' @export
playbase.preprocess::referenceNormalization

#' @importFrom playbase.preprocess svdImpute2
#' @export
playbase.preprocess::svdImpute2


#' Reduce a contrast design to one grouping label per sample
#'
#' All the preprocessing pipeline does with the design is the group-wise
#' missingness filter, which needs a single label per sample. `contrasts` is
#' not validated by its callers and is legitimately NULL for an upload with no
#' comparisons defined yet, so an unusable design yields NULL -- meaning "no
#' design known", which pgx.preprocess() reads as one single group -- rather
#' than an error.
#'
#' @noRd
groupsFromContrasts <- function(samples, contrasts) {
  if (is.null(contrasts)) {
    return(NULL)
  }
  tryCatch(
    apply(contrasts.convertToLabelMatrix(contrasts, samples), 1, paste, collapse = "_"),
    error = function(e) {
      message(
        "[pgx.preprocess] cannot derive sample groups from contrasts: ",
        conditionMessage(e)
      )
      NULL
    }
  )
}

#' @title Preprocess uploaded counts into a normalized expression matrix
#'
#' @description
#' Turns a raw uploaded counts matrix into the normalized log-expression matrix
#' `X` used downstream. The pipeline itself lives in
#' \code{playbase.preprocess::pgx.preprocess}, which describes the design as one
#' grouping label per sample; playbase callers hold the design as a `samples`
#' table plus a `contrasts` matrix, so this wrapper reduces the two and
#' forwards the result.
#'
#' @param counts Raw counts matrix (features x samples). Not modified in place.
#' @param samples Sample annotation table, rows in column order of `counts`.
#' @param contrasts Contrast matrix defining the comparisons. NULL, or a design
#'   that cannot be reduced to sample groups, means no design is known and every
#'   sample is treated as one single group.
#' @param annot Optional annotation table aligned to `counts` rows; subset
#'   alongside NA filtering.
#' @param options Named list of preprocessing settings; see
#'   \code{\link[playbase.preprocess]{pgx.preprocess}} for the full set and its
#'   defaults.
#'
#' @return A list with `counts` (aligned to X rows/cols), `X` (normalized
#'   log-expression), `annot` (subset alongside X, or NULL) and `prior` (log2
#'   prior used).
#'
#' @seealso \code{\link[playbase.preprocess]{pgx.preprocess}}
#'
#' @export
pgx.preprocess <- function(counts,
                           samples = NULL,
                           contrasts = NULL,
                           annot = NULL,
                           options = list()) {
  playbase.preprocess::pgx.preprocess(
    counts,
    groups = groupsFromContrasts(samples, contrasts),
    annot = annot,
    options = options
  )
}
