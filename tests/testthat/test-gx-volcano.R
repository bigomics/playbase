#' Test for gx.volcanoPlot.XY
#'
#'

#' Test for gx.volcanoPlot.LIMMA
test_that("gx.volcanoPlot.LIMMA returns plot", {
  set.seed(100)
  n <- 100
  tab <- data.frame(
    logFC = rnorm(n),
    P.Value = runif(n),
    adj.P.Val = p.adjust(runif(n)),
    AveExpr = runif(n)
  )
  z <- playbase::gx.volcanoPlot.LIMMA(tab)
  z

  expect_equal(class(z), c("scatterD3", "htmlwidget"))
})
