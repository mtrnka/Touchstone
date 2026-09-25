make_context_test_data <- function() {
  tibble::tibble(
    Acc.1 = c("A", "A", "decoy", "A", "B", "A"),
    Acc.2 = c("B", "B", "B", "A", "B", "C"),
    Protein.1 = c("protein A", "protein A", "protein A", "protein A", "protein B", "protein A"),
    Protein.2 = c("protein B", "protein B", "protein B", "protein A", "protein B", "protein C"),
    Species.1 = "TEST",
    Species.2 = "TEST",
    XLink.AA.1 = c(10, 20, 30, 5, 7, 40),
    XLink.AA.2 = c(15, 25, 35, 50, 60, 45),
    DB.Peptide.1 = paste0("PEPA", 1:6),
    DB.Peptide.2 = paste0("PEPB", 1:6),
    Score.Diff = c(20, 18, 8, 15, 14, 17),
    SVM.score = c(2, 1.5, 0.2, 1.2, 1.1, 1.4),
    xlinkClass = c("interProtein", "interProtein", "interProtein", "intraProtein", "intraProtein", "interProtein"),
    Decoy = factor(c("Target", "Target", "Decoy", "Target", "Target", "Target"), levels = c("DoubleDecoy", "Decoy", "Target"))
  ) %>%
    calculatePairs(scalingFactor = 1)
}

test_that("PPI context annotation fuses protein identity but retains decoy origin", {
  result <- annotatePPIContext(
    make_context_test_data(),
    coreThreshold = 1.3,
    candidateThreshold = -Inf,
    scalingFactor = 1
  )
  expect_s3_class(result, "touchstone_ppi_context")
  ab <- result$PPIs %>%
    dplyr::filter(.data$fusedProteinPair == "TEST::protein A::TEST::protein B")
  expect_true(any(ab$Decoy == "Target"))
  expect_true(any(ab$Decoy == "Decoy"))
  expect_true(all(ab$distinctURPContext))
  expect_true(all(ab$fullyDistinctURPs >= 1))
  expect_true(all(ab$bothIntraSupported))
  expect_true(ab$coreSupported[ab$Decoy == "Target"])
  expect_false(ab$coreSupported[ab$Decoy == "Decoy"])
})

test_that("network context excludes the candidate edge itself", {
  result <- annotatePPIContext(
    make_context_test_data(),
    coreThreshold = 1.3,
    scalingFactor = 1
  )
  ab <- result$PPIs %>%
    dplyr::filter(.data$fusedProteinPair == "TEST::protein A::TEST::protein B") %>%
    dplyr::slice(1)
  expect_equal(ab$coreDegreeA, 1)
  expect_equal(ab$coreDegreeB, 0)
  expect_equal(ab$commonCoreNeighbors, 0)
})

test_that("candidate and contextual-support thresholds are independent", {
  result <- annotatePPIContext(
    make_context_test_data(),
    coreThreshold = 1.3,
    candidateThreshold = 0,
    supportThreshold = 1.6,
    scalingFactor = 1
  )
  ab <- result$PPIs %>%
    dplyr::filter(.data$fusedProteinPair == "TEST::protein A::TEST::protein B") %>%
    dplyr::slice(1)
  expect_equal(ab$distinctURPs, 1)
  expect_false(ab$distinctURPContext)
  expect_equal(result$settings$candidateThreshold, 0)
  expect_equal(result$settings$supportThreshold, 1.6)
  expect_true(all(result$PPIs$contextGroup %in% c(
    "distinct", "core-connected", "intra-supported", "context-poor"
  )))
})

test_that("scaled PPI context is retained but marked experimental", {
  result <- NULL
  expect_warning(
    result <- annotatePPIContext(
      make_context_test_data(), coreThreshold = 1, scalingFactor = 5
    ),
    "scaled decoys is experimental"
  )
  expect_equal(result$settings$scalingFactor, 5)
  expect_equal(result$settings$scaledContextCalibration, "experimental")
})

test_that("density-level decoy scaling matches Touchstone count scaling", {
  one.x <- touchstone:::.scaleDecoyEvidence(
    targetDecoy = c(20, 2), doubleDecoy = c(5, 3), scalingFactor = 1
  )
  expect_equal(one.x$ftTT, c(10, 0))
  expect_equal(one.x$ffTT, c(5, 1))

  five.x <- touchstone:::.scaleDecoyEvidence(
    targetDecoy = 100, doubleDecoy = 25, scalingFactor = 5
  )
  expect_equal(five.x$ftTT, 18)
  expect_equal(five.x$ffTT, 1)
})

make_context_calibration_data <- function() {
  groups <- c("distinct", "core-connected", "intra-supported", "context-poor")
  purrr::map_dfr(groups, function(group) {
    target.score <- seq(0.5, 4, length.out = 30)
    decoy.score <- seq(-1, 1, length.out = 8)
    double.score <- seq(-1, 0, length.out = 2)
    score <- c(target.score, decoy.score, double.score)
    tibble::tibble(
      xlinkedProtPair = paste0(group, "-", seq_len(40)),
      SVM.score = score,
      Decoy = factor(
        c(rep("Target", 30), rep("Decoy", 8), rep("DoubleDecoy", 2)),
        levels = c("DoubleDecoy", "Decoy", "Target")
      ),
      xlinkClass = "interProtein",
      contextGroup = group,
      coreSupported = score >= 3
    )
  })
}

test_that("weighted context calibration is monotonic and keeps core PPIs", {
  ppis <- make_context_calibration_data()
  context <- structure(
    list(
      PPIs = ppis,
      URPs = tibble::tibble(),
      settings = list(
        classifier = "SVM.score",
        coreThreshold = 3,
        candidateThreshold = -Inf,
        supportThreshold = 0,
        scalingFactor = 1
      )
    ),
    class = "touchstone_ppi_context"
  )
  first <- classifyPPIContext(
    context,
    targetER = 0.05,
    bootstrapReplicates = 5,
    seed = 7
  )
  second <- classifyPPIContext(
    context,
    targetER = 0.05,
    bootstrapReplicates = 5,
    seed = 7
  )

  expect_s3_class(first, "touchstone_ppi_results")
  expect_true(all(c(
    "contextPEP", "contextQValue", "selectionFrequency",
    "classificationTier", "classified"
  ) %in% names(first$PPIs)))
  expect_true(all(first$PPIs$classified[first$PPIs$coreSupported]))
  expect_true(all(vapply(
    first$model$curves,
    function(curve) all(diff(curve$pep) <= sqrt(.Machine$double.eps)),
    logical(1)
  )))
  expect_identical(
    first$bootstrap$perPPI$selectionFrequency,
    second$bootstrap$perPPI$selectionFrequency
  )
  expect_equal(first$settings$calibration, "weighted-isotonic-context-PEP")

  set.seed(81)
  expected.random <- stats::runif(1)
  set.seed(81)
  invisible(classifyPPIContext(
    context,
    targetER = 0.05,
    bootstrapReplicates = 2,
    seed = 7
  ))
  expect_equal(stats::runif(1), expected.random)
})

test_that("context q-values accept or reject complete PEP plateaus", {
  pep <- c(rep(0.001, 10), rep(0.03, 100))
  decoy.class <- rep("Target", length(pep))
  q.value <- contextPEPQValues(pep, decoy.class)

  expect_equal(contextPEPThreshold(pep, decoy.class, 0.02), 0.001)
  expect_length(unique(q.value[pep == 0.03]), 1)
  expect_gt(unique(q.value[pep == 0.03]), 0.02)
})

test_that("weighted decreasing PAVA pools local score reversals", {
  fitted <- weightedDecreasingPAVA(
    c(0.4, 0.3, 0.35, 0.1),
    rep(1, 4)
  )
  expect_equal(fitted, c(0.4, 0.325, 0.325, 0.1))
  expect_true(all(diff(fitted) <= 0))
})
