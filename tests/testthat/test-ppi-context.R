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
