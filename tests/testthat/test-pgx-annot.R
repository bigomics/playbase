## An ALIAS key that is also another gene's official symbol must resolve to
## that gene, not to whichever gene sorts first: "TPO" is thyroid peroxidase
## (and an alias of THPO), "EBF3" is EBF3 (and an alias of MAPRE3).
test_that("an official symbol wins over another gene's alias", {
  annot <- data.frame(
    ALIAS = c("TPO", "TPO", "EBF3", "EBF3", "GAPDH", "OLDNAME", "AMBIG", "AMBIG"),
    SYMBOL = c("THPO", "TPO", "MAPRE3", "EBF3", "GAPDH", "NEWNAME", "GENEA", "GENEB"),
    GENENAME = c("thrombopoietin", "thyroid peroxidase", "MAPRE3 title", "EBF3 title",
                 "gapdh", "renamed gene", "a", "b")
  )
  out <- resolve_alias_rows(annot, "ALIAS")
  expect_identical(out$SYMBOL[out$ALIAS == "TPO"], "TPO")
  expect_identical(out$SYMBOL[out$ALIAS == "EBF3"], "EBF3")
  expect_identical(out$SYMBOL[out$ALIAS == "GAPDH"], "GAPDH")
  ## An alias of exactly one gene still resolves; an ambiguous one does not.
  expect_identical(out$SYMBOL[out$ALIAS == "OLDNAME"], "NEWNAME")
  expect_false("AMBIG" %in% out$ALIAS)
  ## Other keytypes pass through untouched.
  expect_identical(resolve_alias_rows(annot, "SYMBOL"), annot)
})

test_that("gene annotation by ALIAS keeps the official symbol's gene", {
  skip_if_not_installed("org.Hs.eg.db")
  orgdb <- org.Hs.eg.db::org.Hs.eg.db
  keys <- c("TPO", "EBF3", "TNRC18", "AHRR", "GAPDH")
  a <- suppressMessages(AnnotationDbi_select_2pass(
    orgdb, keys = keys, columns = c("SYMBOL", "GENENAME"), keytype = "ALIAS"
  ))
  expect_identical(a$SYMBOL, keys)
  expect_identical(a$GENENAME[1:2], c("thyroid peroxidase", "EBF transcription factor 3"))
})

## The loose match (case and punctuation folded) is a fallback only: an exact
## key must win, or "SF3B1" becomes SF3B2's alias "SF3b1" and "NKX6-1" the
## alias "NKX6.1".
test_that("probe names keep an exact key over a case- or punctuation-folded one", {
  skip_if_not_installed("org.Hs.eg.db")
  orgdb <- org.Hs.eg.db::org.Hs.eg.db
  keys <- c("SF3B1", "SLIT2", "TIFA", "NKX6-1", "VMA21")
  expect_identical(unname(match_probe_names(keys, orgdb, "ALIAS")), keys)
  ## Folding still rescues a key that has no exact match.
  expect_identical(unname(match_probe_names("gapdh", orgdb, "SYMBOL")), "GAPDH")
})
