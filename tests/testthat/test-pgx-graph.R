test_that("getGOgraph returns an igraph with many vertices", {
  skip_if_not_installed("GO.db")
  skip_if(inherits(
    tryCatch(AnnotationDbi::toTable(GO.db::GOTERM), error = function(e) e),
    "error"
  ))
  g <- playbase::getGOgraph()
  expect_s3_class(g, "igraph")
  expect_gt(igraph::vcount(g), 1000)
})
