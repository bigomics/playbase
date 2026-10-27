test_that(".clean_probe_names strips numeric prefix before UniProt accession", {
  probes <- c("3_Q13257", "1_A0A075B6H7", "Q9Y2S6;3_A0A024R1R8", "7157_at", "NM_000546", "P04637-2")
  expected <- c("Q13257", "A0A075B6H7", "Q9Y2S6", "7157", "NM_000546", "P04637")
  expect_equal(unname(playbase:::.clean_probe_names(probes)), expected)
})
