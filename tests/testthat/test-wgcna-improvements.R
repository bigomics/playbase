## Regression tests for bugs fixed in the "wgcna improvements" PR:
## - wgcna.getTopMEfeatures lost gene names on colSds(), breaking sortby="sd"
## - wgcna.plotModuleTrajectory crashed when the default/user start or end
##   module fell outside the selected `k` (most commonly MEgrey)

test_that("wgcna.getTopMEfeatures sorts by sd using gene names, not position", {
  genes <- paste0("g", 1:4)
  colors <- c("blue", "blue", "blue", "blue")
  names(colors) <- genes

  ## columns (genes) have deliberately different spreads, in an order that
  ## does NOT match the sd ranking, so a names-lost regression would be caught
  datExpr <- matrix(
    c(1, 1, 1, 1,
      2, 2, 2, 2,
      1, 5, 1, 1,
      9, 2, 8, 1),
    nrow = 4, byrow = TRUE
  )
  colnames(datExpr) <- genes

  M <- matrix(0.5, nrow = 4, ncol = 1, dimnames = list(genes, "MEblue"))
  wgcna <- list(net = list(colors = colors), datExpr = datExpr, stats = list(moduleMembership = M))

  top <- wgcna.getTopMEfeatures(wgcna, n = 4, sortby = "sd")

  sdx <- matrixStats::colSds(datExpr)
  names(sdx) <- colnames(datExpr)
  expected_order <- names(sort(-sdx))

  expect_equal(top[["MEblue"]], expected_order)
})

test_that("wgcna.plotModuleTrajectory default start skips grey even when it is the largest module", {
  set.seed(1)
  MEs <- matrix(rnorm(4 * 6), nrow = 6, ncol = 4)
  colnames(MEs) <- c("MEgrey", "MEblue", "MEyellow", "MEturquoise")
  rownames(MEs) <- paste0("s", 1:6)
  colors <- c(rep("grey", 20), rep("blue", 3), rep("yellow", 2), rep("turquoise", 1))
  wgcna <- list(net = list(MEs = MEs, colors = colors))

  ee <- wgcna.plotModuleTrajectory(wgcna, start = NULL, plot = FALSE)

  expect_false("MEgrey" %in% ee)
  expect_true(all(ee %in% c("MEblue", "MEyellow", "MEturquoise")))
})

test_that("wgcna.plotModuleTrajectory can be called without start at all", {
  set.seed(1)
  MEs <- matrix(rnorm(4 * 6), nrow = 6, ncol = 4)
  colnames(MEs) <- c("MEgrey", "MEblue", "MEyellow", "MEturquoise")
  rownames(MEs) <- paste0("s", 1:6)
  colors <- c(rep("grey", 20), rep("blue", 3), rep("yellow", 2), rep("turquoise", 1))
  wgcna <- list(net = list(MEs = MEs, colors = colors))

  ## start has no default in R unless the signature gives it one -- this must
  ## not error with "argument \"start\" is missing, with no default"
  ee <- wgcna.plotModuleTrajectory(wgcna, plot = FALSE)
  expect_false("MEgrey" %in% ee)
})

test_that("wgcna.plotModuleTrajectory with end does not distort the non-rho layout", {
  set.seed(1)
  MEs <- matrix(rnorm(4 * 6), nrow = 6, ncol = 4)
  colnames(MEs) <- c("MEgrey", "MEblue", "MEyellow", "MEturquoise")
  rownames(MEs) <- paste0("s", 1:6)
  colors <- c(rep("grey", 20), rep("blue", 3), rep("yellow", 2), rep("turquoise", 1))
  wgcna <- list(net = list(MEs = MEs, colors = colors))

  ## the layout positions (passed to plot()) must come from the UNMODIFIED
  ## distance matrix, not the one forced-to-999 for the TSP "end" constraint
  ME <- t(MEs)[c("MEblue", "MEyellow", "MEturquoise"), ]
  expected_pos <- svd(dist(ME))$u[, 1:2]

  captured <- NULL
  testthat::local_mocked_bindings(
    plot.default = function(x, ...) {
      captured <<- x
      graphics::plot.new()
      graphics::plot.window(range(x[, 1]), range(x[, 2]))
    },
    .package = "graphics"
  )
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)

  wgcna.plotModuleTrajectory(wgcna, start = "MEblue", end = "MEyellow",
                              posmat = "dist", plot = TRUE)

  expect_equal(unname(captured), unname(expected_pos))
})

test_that("wgcna.plotModuleTrajectory errors clearly when start/end is not in the selected modules", {
  set.seed(1)
  MEs <- matrix(rnorm(4 * 6), nrow = 6, ncol = 4)
  colnames(MEs) <- c("MEgrey", "MEblue", "MEyellow", "MEturquoise")
  rownames(MEs) <- paste0("s", 1:6)
  colors <- c(rep("grey", 20), rep("blue", 3), rep("yellow", 2), rep("turquoise", 1))
  wgcna <- list(net = list(MEs = MEs, colors = colors))

  expect_error(
    wgcna.plotModuleTrajectory(wgcna, start = "MEgrey", plot = FALSE),
    "not found among selected modules"
  )
  expect_error(
    wgcna.plotModuleTrajectory(wgcna, start = "MEblue", end = "MEgrey", plot = FALSE),
    "not found among selected modules"
  )
})

test_that("wgcna.create_module_titles picks up the extras passed in, not a hardcoded slot", {
  testthat::local_mocked_bindings(
    ai.ask = function(question, model, ...) paste("title for:", question)
  )

  extras <- list(MEblue = list(summary = "blue module summary"))
  desc <- wgcna.create_module_titles(list(), model = "x", extras = extras, from = "extras")

  expect_equal(names(desc), "MEblue")
  expect_true(grepl("blue module summary", desc[["MEblue"]]))
})
