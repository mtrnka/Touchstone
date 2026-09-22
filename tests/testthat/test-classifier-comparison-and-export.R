test_that("formatXLTable selects and preserves an available classifier", {
  scored <- tibble::tibble(
    xlinkedResPair = factor(c("R1", "R2", "R3")),
    Score.Diff = c(2, 10, 5),
    experimental.score = c(20, 1, 10),
    Decoy = "Target"
  )

  fallback <- formatXLTable(scored)
  expect_identical(fallback$Score.Diff, c(10, 5, 2))

  experimental <- formatXLTable(
    scored,
    classifier = "experimental.score"
  )
  expect_identical(experimental$experimental.score, c(20, 10, 1))
  expect_true("experimental.score" %in% names(experimental))
  expect_error(
    formatXLTable(scored, classifier = "missing.score"),
    "no classifier column"
  )
})

test_that("formatXLTable inherits the classifier from prepared results", {
  scored <- tibble::tibble(
    xlinkedResPair = factor(c("R1", "R2")),
    Score.Diff = c(1, 2),
    experimental.score = c(10, 5),
    Decoy = "Target"
  )
  prepared <- structure(
    list(
      data = scored,
      settings = list(classifier = "experimental.score")
    ),
    class = "touchstone_results"
  )

  output <- formatXLTable(prepared)
  expect_identical(output$experimental.score, c(10, 5))
})

test_that("makeXiNetFile supports fallback and supplied classifiers", {
  scored <- tibble::tibble(
    Score.Diff = c(2, 10),
    experimental.score = c(20, 1),
    Acc.1 = c("P1", "P2"),
    Acc.2 = c("P2", "P3"),
    XLink.AA.1 = c(10, 20),
    XLink.AA.2 = c(30, 40)
  )

  fallback <- makeXiNetFile(scored, flavor = "xiView")
  expect_identical(fallback$Score, scored$Score.Diff)
  expect_named(
    fallback,
    c("Score", "Protein1", "Protein2", "AbsPos1", "AbsPos2")
  )

  experimental <- makeXiNetFile(
    scored,
    flavor = "xiNet",
    classifier = "experimental.score"
  )
  expect_identical(experimental$Score, scored$experimental.score)
  expect_error(
    makeXiNetFile(dplyr::select(scored, -Score.Diff, -experimental.score)),
    "must contain SVM.score, Score.Diff"
  )
})

test_that("makeXiNetFile inherits a prepared classifier", {
  scored <- tibble::tibble(
    experimental.score = c(3, 4),
    Acc.1 = c("P1", "P2"),
    Acc.2 = c("P2", "P3"),
    XLink.AA.1 = c(10, 20),
    XLink.AA.2 = c(30, 40)
  )
  prepared <- structure(
    list(
      data = scored,
      settings = list(classifier = "experimental.score")
    ),
    class = "touchstone_results"
  )

  output <- makeXiNetFile(prepared)
  expect_identical(output$Score, scored$experimental.score)
})

test_that("classifier comparison reuses separated FDR error tables", {
  n <- 60
  scored <- tibble::tibble(
    xlinkedResPair = factor(paste0("R", seq_len(n))),
    Score.Diff = seq_len(n),
    experimental.score = c(seq_len(n / 2), rev(seq_len(n / 2))),
    xlinkClass = rep(c("interProtein", "intraProtein"), each = n / 2),
    Decoy = factor(
      rep(c(rep("Target", 8), "Decoy", "DoubleDecoy"), 6),
      levels = c("DoubleDecoy", "Decoy", "Target")
    )
  )

  comparison <- compareClassifiers(
    scored,
    classifiers = c("Score.Diff", "experimental.score"),
    summarizationLevel = "urp",
    scalingFactor = 1
  )

  expect_s3_class(comparison, "tbl_df")
  expect_setequal(
    comparison$classifier,
    c("Score.Diff", "experimental.score")
  )
  expect_setequal(
    comparison$xlinkClass,
    c("interProtein", "intraProtein")
  )
  expect_true(all(c("threshold", "fdr", "hits") %in% names(comparison)))
  expect_identical(attr(comparison, "scalingFactor"), 1)
  expect_s3_class(
    plotClassifierComparison(comparison, maxFDR = 1),
    "ggplot"
  )
})
