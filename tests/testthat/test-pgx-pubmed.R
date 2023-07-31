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
  expect_match(link,  
  "<a href='https://www.ncbi.nlm.nih.gov/pubmed/12345678' target='_blank'>PMID:12345678</a>")
  
})
