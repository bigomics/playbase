test_that(".clean_probe_names strips numeric prefix before UniProt accession", {
  probes <- c("3_Q13257", "1_A0A075B6H7", "Q9Y2S6;3_A0A024R1R8", "7157_at", "NM_000546", "P04637-2")
  expected <- c("Q13257", "A0A075B6H7", "Q9Y2S6", "7157", "NM_000546", "P04637")
  expect_equal(unname(playbase:::.clean_probe_names(probes)), expected)
})

test_that("normalizeOrganism maps all fleet aliases to canonical scientific names", {
  aliases <- c("Human", "human", "Homo sapiens", "Mouse", "mouse", "Mus musculus",
    "Rat", "rat", "Rattus norvegicus")
  expected <- c(rep("Homo sapiens", 3), rep("Mus musculus", 3), rep("Rattus norvegicus", 3))
  expect_equal(unname(normalizeOrganism(aliases)), expected)

  # regression for #522: a vector of organism aliases must stay a scalar
  # decision once reduced with all()/any(), as pgx.add_GMT, pgx.initialize
  # and compute_drugActivityEnrichment rely on.
  expect_false(all(normalizeOrganism(c("Homo sapiens", "dog")) %in% c("Homo sapiens", "Mus musculus", "Rattus norvegicus")))
  expect_true(all(normalizeOrganism(c("Homo sapiens", "Mus musculus")) %in% c("Homo sapiens", "Mus musculus", "Rattus norvegicus")))
})
