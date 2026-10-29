# playbase (development version)

* Numeric single-cell feature, mitochondrial, and hemoglobin QC thresholds now
  activate the corresponding cell filters in `pgx.createSeuratObject()`.

## Breaking changes

* `is.xxl()` and `mofa.log1s()` removed from the public API. Both were exported
  but had zero known callers in `playbase`, `omicsplayground`, or `omicspgx`.
  If you depended on either, they are gone — there is no replacement.
* Azimuth is no longer a dependency. `pgx.runAzimuth()` is replaced by
  `pgx.predictCellTypes(counts, organism, reference)`, which maps cells onto
  the same Azimuth reference atlases with Seurat directly. Like Azimuth it
  leaves human symbols as they are; other gene IDs (mouse, Ensembl) are
  converted with `getProbeAnnotation()` and `rename_by2()`. It now stops on
  failure instead of returning `NULL`. This drops 24 packages (~1.2 GB, including the hg38
  BSgenome and EnsDb.Hsapiens.v86) and the compiler they needed.
* `pgx.createPGX()` and `pgx.createSingleCellPGX()`: argument `azimuth_ref` is
  renamed `celltype_ref`.
* `glmGamPoi` and `SeuratObject` are now declared in Imports. Both were used
  directly but only installed because Azimuth pulled them in.

## Preprocessing extraction (in progress)

This is the first entry of a multi-step port that moves the preprocessing layer
(`pgx-preprocess.R`, `pgx-normalize.R`, `pgx-impute.R`, `pgx-outlier.R`) into a
separate package, `playbase.preprocess`. Further entries will follow as that
work lands, including a `pgx.preprocess()` signature change and the extraction
itself. Neither has happened yet as of this entry.
