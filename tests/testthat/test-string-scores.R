test_that("STRING scores are queried once per target protein pair", {
  query.count <- 0L
  fake.db <- list(
    map = function(dat, column, removeUnmappedRows) {
      expect_identical(column, "acc")
      expect_false(removeUnmappedRows)
      tibble::tibble(
        acc = dat$acc,
        STRING_id = unname(c(
          P1 = "s1", P2 = "s2", P3 = "s3", P4 = NA_character_
        )[dat$acc])
      )
    },
    get_interactions = function(ids) {
      query.count <<- query.count + 1L
      if (setequal(ids, c("s1", "s2"))) {
        return(tibble::tibble(
          from = "s2", to = "s1", combined_score = 800
        ))
      }
      tibble::tibble(
        from = character(), to = character(), combined_score = numeric()
      )
    }
  )
  dat <- tibble::tibble(
    row = letters[1:6],
    Acc.1 = c("P1", "P2", "P1", "P1", "P4", "DECOY1"),
    Acc.2 = c("P2", "P1", "P3", "P1", "P2", "P2"),
    xlinkedProtPair = c(
      "P1::P2", "P1::P2", "P1::P3", "P1::P1", "P2::P4",
      "decoy@one::P2"
    ),
    xlinkClass = c(
      "interProtein", "interProtein", "interProtein", "intraProtein",
      "interProtein", "interProtein"
    ),
    Decoy = factor(
      c(rep("Target", 5), "Decoy"),
      levels = c("DoubleDecoy", "Decoy", "Target")
    ),
    string.score = -1
  )

  result <- .getStringScoresWithDB(dat, fake.db)

  expect_identical(result$row, dat$row)
  expect_equal(result$string.score[1:2], c(800, 800))
  expect_true(all(is.na(result$string.score[3:6])))
  expect_identical(query.count, 2L)
})

test_that("getStringScores requires existing Touchstone pair identifiers", {
  expect_error(
    getStringScores(tibble::tibble(Acc.1 = "P1", Acc.2 = "P2"), 511145),
    "Run calculatePairs"
  )
})
