## Production reaches playbase through `playbase::pgx.createPGX()` off a bare
## Rscript -- upload_module_computepgx.R:1324 -> processx -> bin/pgxcreate_op.R
## -- and there is no library() anywhere in that chain. `Depends:` ATTACHES on
## library() but only LOADS on `pkg::fun()`, so every unqualified call into
## playbase.preprocess resolves through playbase's imports env there, and
## through the search path everywhere else.
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
    stdout = TRUE, stderr = TRUE,
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
    skip("playbase is dev-loaded; a detached process has no installed copy to ::")
  }
}

test_that("every playbase.preprocess import resolves without attaching", {
  skip_unless_playbase_installed()

  res <- run_detached(c(
    'ns  <- loadNamespace("playbase")',
    'imp <- parent.env(ns)',
    'nms <- getNamespaceImports(ns)[["playbase.preprocess"]]',
    'if (!length(nms)) stop("playbase declares no playbase.preprocess imports")',
    'bad <- nms[!vapply(nms, exists, logical(1), envir = imp, inherits = TRUE)]',
    'if (length(bad)) stop("unresolved without attaching: ", paste(bad, collapse = ", "))',
    'cat("RESOLVED", length(nms), "\n")'
  ))

  expect_equal(res$status, 0L, info = res$output)
  expect_match(res$output, "RESOLVED [0-9]+")
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

## Round 2 of the code review falsified the test above: it asserts that the
## imports playbase DECLARES resolve, never that the ones it NEEDS are
## declared. Strip all four @importFrom directives and it still passes, because
## the other 24 entries resolve fine. It cannot fail on the defect it was
## written for.
##
## This one can. It reads the source rather than the namespace: every leaf
## export called UNQUALIFIED anywhere in playbase's R/ must appear in playbase's
## imports. No child process, no installed copy, no attachment -- it is a
## statement about the source tree and it runs everywhere, including here.
test_that("every unqualified leaf call is a declared import", {
  skip_if_not_installed("playbase.preprocess")

  leaf <- getNamespaceExports("playbase.preprocess")
  rdir <- testthat::test_path("..", "..", "R")
  skip_if_not(dir.exists(rdir))

  src <- unlist(lapply(list.files(rdir, "\\.[rR]$", full.names = TRUE), readLines))
  ## drop comments and anything already qualified with a namespace
  src <- sub("#.*$", "", src)
  src <- gsub("[A-Za-z_.][A-Za-z0-9_.]*::+", "", src)

  called <- Filter(function(f) {
    any(grepl(paste0("(^|[^A-Za-z0-9_.])", gsub("\\.", "\\\\.", f), "\\s*\\("), src))
  }, leaf)

  ## From NAMESPACE, not from the loaded namespace: under pkgload the imports
  ## env does not report what an installed build would, and this must be a
  ## statement about the source tree.
  nsfile <- testthat::test_path("..", "..", "NAMESPACE")
  skip_if_not(file.exists(nsfile))
  ns <- readLines(nsfile)
  declared <- sub(
    "^importFrom\\(playbase\\.preprocess,([^)]+)\\)$", "\\1",
    grep("^importFrom\\(playbase\\.preprocess,", ns, value = TRUE)
  )
  ## playbase defines a wrapper of its own for some leaf names (pgx.preprocess
  ## is one): a call to those resolves inside playbase and needs no import.
  own <- ls(asNamespace("playbase"), all.names = TRUE)
  missing <- setdiff(setdiff(called, declared), own)

  expect_identical(
    missing, character(0),
    info = paste0(
      "called unqualified in playbase/R/ but not imported: ",
      paste(missing, collapse = ", "),
      ". These resolve only when playbase.preprocess is ATTACHED; production ",
      "reaches playbase through pkg::fun(), which only LOADS. See review R0."
    )
  )
})
