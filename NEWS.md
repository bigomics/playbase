# playbase (development version)

## Breaking changes

* `is.xxl()` and `mofa.log1s()` removed from the public API. Both were exported
  but had zero known callers in `playbase`, `omicsplayground`, or `omicspgx`.
  If you depended on either, they are gone — there is no replacement.

## Preprocessing extraction (in progress)

This is the first entry of a multi-step port that moves the preprocessing layer
(`pgx-preprocess.R`, `pgx-normalize.R`, `pgx-impute.R`, `pgx-outlier.R`) into a
separate package, `playbase.preprocess`. Further entries will follow as that
work lands, including a `pgx.preprocess()` signature change and the extraction
itself. Neither has happened yet as of this entry.
