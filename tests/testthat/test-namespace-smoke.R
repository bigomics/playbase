## Production reaches playbase through `playbase::pgx.createPGX()` from a bare
## Rscript. Playbase declares the leaf in Imports and calls its API through the
## explicit `playbase.preprocess::` namespace without attaching either package.
##
## This suite and `R CMD check` both attach. That is how four missing
## `@importFrom` directives sat behind 750 green tests on a code path that could
## not execute at all (code review round 1, R0). The two tests below are the
## only ones in the suite that do not attach: they shell out to a clean R
## process, exactly as production does, so attachment can never again be the
## thing holding pgx.createPGX() up.
##
## They skip under devtools::test(), where pkgload has installed nothing for a
## child process to `::` into, and run under `R CMD check`, which installs
## first, and against any installed playbase.

## A child R process with this session's library path and NOTHING attached.
run_detached <- function(code) {
  script <- tempfile(fileext = ".R")
  on.exit(unlink(script), add = TRUE)
  writeLines(code, script)
  out <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", shQuote(script)),
    stdout = TRUE,
    stderr = TRUE,
    env = paste0("R_LIBS=", paste(.libPaths(), collapse = .Platform$path.sep))
  ))
  status <- attr(out, "status")
  list(
    status = if (is.null(status)) 0L else as.integer(status),
    output = paste(out, collapse = "\n")
  )
}

## Skips only when there is no INSTALLED playbase for a detached process to
## load -- never because the call failed. A skip here means "not applicable",
## and says which case it is.
skip_unless_playbase_installed <- function() {
  dir <- find.package("playbase", quiet = TRUE)
  if (!length(dir)) {
    skip("playbase is not on this session's library path")
  }
  if (!file.exists(file.path(dir[1], "Meta", "package.rds"))) {
    skip(
      "playbase is dev-loaded; a detached process has no installed copy to ::"
    )
  }
}

test_that("playbase loads the leaf without attaching it", {
  skip_unless_playbase_installed()

  res <- run_detached(c(
    'loadNamespace("playbase")',
    'stopifnot(!"package:playbase" %in% search())',
    'stopifnot(!"package:playbase.preprocess" %in% search())',
    'fn <- getExportedValue("playbase.preprocess", "pgx.preprocess")',
    'stopifnot(is.function(fn))',
    'stopifnot(isNamespaceLoaded("playbase.preprocess"))',
    'stopifnot(!"package:playbase.preprocess" %in% search())',
    'cat("LOADED DETACHED\n")'
  ))

  expect_equal(res$status, 0L, info = res$output)
  expect_match(res$output, "LOADED DETACHED")
})

test_that("pgx.createPGX runs in a process that never attached playbase", {
  ## No skip_on_cran(): it is off by default outside devtools, so it would have
  ## skipped this test in every `R CMD check` -- the one run that matters.
  skip_unless_playbase_installed()

  res <- run_detached(c(
    '## No library() anywhere below -- that is the whole point of the test.',
    'set.seed(1)',
    'n <- 12L; m <- 4L',
    'counts <- matrix(rpois(n * m, 200), nrow = n,',
    '  dimnames = list(paste0("G", seq_len(n)), paste0("S", seq_len(m))))',
    'grp <- rep(c("a", "b"), each = m / 2)',
    'samples <- data.frame(group = grp, row.names = colnames(counts))',
    'contrasts <- data.frame(b_vs_a = grp, row.names = colnames(counts))',
    'pgx <- playbase::pgx.createPGX(',
    '  counts = counts, samples = samples, contrasts = contrasts,',
    '  organism = "No organism", datatype = "RNA-seq",',
    '  preprocess = list(remove_outliers = FALSE, impute = FALSE),',
    '  norm_method = "CPM", add.gmt = FALSE, convert.hugo = FALSE,',
    '  prune.samples = FALSE)',
    'stopifnot(is.list(pgx), !is.null(pgx$X), !is.null(pgx$counts))',
    'cat("CREATED", nrow(pgx$counts), "x", ncol(pgx$counts),',
    '  "X", nrow(pgx$X), "x", ncol(pgx$X), "\n")'
  ))

  expect_equal(res$status, 0L, info = res$output)
  expect_match(res$output, "CREATED 12 x 4 X 12 x 4")
})
