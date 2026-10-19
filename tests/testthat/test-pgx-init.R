## Run pgx.initialize on a minimal pgx up to the gene clustering step, which
## is mocked to hand back the pgx state reached so far.
init_until_clustering <- function(pgx) {
  testthat::local_mocked_bindings(
    .pb_libx = function(...) list(),
    pgx.clusterGenes = function(pgx, ...) {
      stop(structure(
        class = c("pgx_captured", "error", "condition"),
        list(message = "captured", call = NULL, pgx = pgx)
      ))
    }
  )
  tryCatch(
    pgx.initialize(pgx),
    pgx_captured = function(e) e$pgx
  )
}

minimal_pgx <- function(gene_name, organism, ...) {
  n <- length(gene_name)
  samples <- data.frame(group = rep(c("a", "b"), length.out = 4), row.names = paste0("s", 1:4))
  X <- matrix(stats::rnorm(n * 4), n, 4, dimnames = list(gene_name, rownames(samples)))
  list(
    samples = samples,
    genes = data.frame(gene_name = gene_name, gene_title = gene_name, row.names = gene_name),
    GMT = NULL, gx.meta = list(), tsne2d = NULL, X = X,
    counts = 2**X,
    model.parameters = list(group = setNames(samples$group, rownames(samples))),
    contrasts = matrix(c("a", "b", "a", "b"), 4, 1, dimnames = list(rownames(samples), "a_vs_b")),
    organism = organism, ...
  )
}

test_that("mouse pgx without version gets mapped ortholog or NA, never upper-cased", {
  testthat::local_mocked_bindings(
    getHumanOrtholog = function(organism, symbols, ...) {
      data.frame(symbol = symbols, ortholog = c("TP53", NA), orthologs = c("TP53", NA))
    }
  )
  pgx <- suppressMessages(init_until_clustering(minimal_pgx(c("Trp53", "Xyz1"), "Mouse")))
  expect_equal(pgx$genes$ortholog, c("TP53", NA))
  expect_false(identical(pgx$genes$ortholog, toupper(pgx$genes$gene_name)))
})

test_that("human pgx without version keeps gene_name as given in ortholog", {
  pgx <- suppressMessages(init_until_clustering(minimal_pgx(c("Tp53", "ACTB"), "Human")))
  expect_equal(pgx$genes$ortholog, c("Tp53", "ACTB"))
})

test_that("empty human_ortholog stays NA with a message when no ortholog is found", {
  testthat::local_mocked_bindings(
    getHumanOrtholog = function(organism, symbols, ...) {
      data.frame(symbol = symbols, ortholog = NA_character_, orthologs = NA_character_)
    }
  )
  pgx <- minimal_pgx(c("Trp53", "Xyz1"), "Mouse", version = "3")
  pgx$genes$symbol <- pgx$genes$gene_name
  pgx$genes$human_ortholog <- NA_character_
  expect_message(
    res <- init_until_clustering(pgx),
    "no human ortholog",
    all = FALSE
  )
  expect_true(all(is.na(res$genes$human_ortholog)))
})
