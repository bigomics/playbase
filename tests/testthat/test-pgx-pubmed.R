#' Test for pmid.getGeneContext
#'
#'

#' Test for pmid.getPubMedContext
#'
#'

#' Test for pmid.getGeneContext
#'

#' Test for pmid.buildMatrix
#'
#'


#' Test for pubmedlink
test_that("pubmedlink generates valid hyperlink", {
  # Generate sample PMID
  pmid <- "12345678"

  # Call function
  link <- playbase::pubmedlink(pmid)

  # Check class
  expect_type(link, "character")

  # Check HTML structure
  expect_match(
    link,
    "<a href='https://www.ncbi.nlm.nih.gov/pubmed/12345678' target='_blank'>PMID:12345678</a>"
  )
})

#' Gene-name lookup used by pmid.getPubMedContext
test_that("org.Hs.eg.db maps a gene symbol to its gene name", {
  skip_if_not_installed("org.Hs.eg.db")
  skip_if(inherits(
    tryCatch(get("TP53", org.Hs.eg.db::org.Hs.egALIAS2EG), error = function(e) e),
    "error"
  ))
  gname <- get(get("TP53", org.Hs.eg.db::org.Hs.egALIAS2EG), org.Hs.eg.db::org.Hs.egGENENAME)
  expect_match(gname, "tumor protein p53", ignore.case = TRUE)
})
