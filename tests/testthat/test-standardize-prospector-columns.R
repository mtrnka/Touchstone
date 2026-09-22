make_named_row <- function(column_names) {
  result <- as.data.frame(
    matrix(seq_along(column_names), nrow = 1),
    check.names = FALSE
  )
  names(result) <- column_names
  result
}

test_that("selected paired Search Compare columns receive canonical names", {
  raw <- make_named_row(c(
    "Acc #...1", "Acc #...2",
    "Protein Name...3", "Protein Name...4",
    "Species...5", "Species...6",
    "% Bond Cleavage...7", "% Bond Cleavage...8",
    "MSMS Ions...9", "MSMS Ions...10",
    "Score Diff", "Unrelated Column"
  ))

  standardized <- touchstone:::standardizeProspectorColumns(raw)

  expect_named(standardized, c(
    "Acc.1", "Acc.2",
    "Protein.1", "Protein.2",
    "Species.1", "Species.2",
    "Perc.Bond.Cleavage.1", "Perc.Bond.Cleavage.2",
    "MSMS.Ions.1", "MSMS.Ions.2",
    "Score.Diff", "Unrelated.Column"
  ))
})

test_that("optional Search Compare columns may be omitted", {
  raw <- make_named_row(c("Acc #...1", "Acc #...2", "Score Diff"))

  standardized <- touchstone:::standardizeProspectorColumns(raw)

  expect_named(standardized, c("Acc.1", "Acc.2", "Score.Diff"))
})

test_that("Prospector peptide suffixes are recognized", {
  raw <- make_named_row(c(
    "Acc #_1", "Acc #_2",
    "Protein Name_1", "Protein Name_2"
  ))

  standardized <- touchstone:::standardizeProspectorColumns(raw)

  expect_named(
    standardized,
    c("Acc.1", "Acc.2", "Protein.1", "Protein.2")
  )
})

test_that("already canonical columns remain unchanged", {
  canonical <- make_named_row(c(
    "Acc.1", "Acc.2", "Protein.1", "Protein.2", "Score.Diff"
  ))

  expect_identical(
    names(touchstone:::standardizeProspectorColumns(canonical)),
    names(canonical)
  )
})

test_that("ambiguous required paired columns are rejected", {
  one_accession <- make_named_row(c("Acc #", "Score Diff"))
  three_accessions <- make_named_row(c(
    "Acc #...1", "Acc #...2", "Acc #...3"
  ))

  expect_error(
    touchstone:::standardizeProspectorColumns(one_accession),
    "Expected exactly two accession columns, but found 1"
  )
  expect_error(
    touchstone:::standardizeProspectorColumns(three_accessions),
    "Expected exactly two accession columns, but found 3"
  )
})

test_that("matched-intensity variants produce the same canonical feature", {
  reported_intensity <- data.frame(Match.Int = c(25, 80))
  calculated_intensity <- data.frame(
    Num.Pks = c(20, 10),
    Num.Unmat = c(15, 2)
  )

  expect_equal(calculatePercentMatched(reported_intensity)$percMatched, c(25, 80))
  expect_equal(calculatePercentMatched(calculated_intensity)$percMatched, c(0.25, 0.8))
})
