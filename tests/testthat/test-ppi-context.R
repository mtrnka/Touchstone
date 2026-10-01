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
  expect_true(all(c("ppiID", "proteinInferenceStatus") %in%
                    names(result$PPIs)))
  expect_true(all(result$PPIs$proteinInferenceStatus == "not-assessed"))
  expect_true(all(as.character(ab$contextGroup) ==
                    "Distinct URP + intra-supported"))
  expect_false("DB.Peptide.1" %in% names(result$PPIs))
  expect_true(all(c("ppiID", "candidateURP", "contextualSupportURP") %in%
                    names(result$URPs)))
  expect_true(all(c("ppiID", "candidateURP", "contextualSupportURP") %in%
                    names(result$CSMs)))
  ppi.csm.counts <- result$CSMs %>%
    dplyr::count(.data$ppiID, name = "expectedCSMs")
  ppi.urp.counts <- result$URPs %>%
    dplyr::count(.data$ppiID, name = "expectedURPs")
  expect_equal(
    dplyr::left_join(result$PPIs, ppi.csm.counts, by = "ppiID")$numCSMs,
    dplyr::left_join(result$PPIs, ppi.csm.counts, by = "ppiID")$expectedCSMs
  )
  expect_equal(
    dplyr::left_join(result$PPIs, ppi.urp.counts, by = "ppiID")$numURPs,
    dplyr::left_join(result$PPIs, ppi.urp.counts, by = "ppiID")$expectedURPs
  )
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
    "Distinct URP + intra-supported",
    "Core-connected + intra-supported",
    "Intra-supported",
    "No additional context"
  )))
})

test_that("PPI context groups preserve intra-support interactions", {
  groups <- assignPPIContextGroup(
    distinctURPContext = c(TRUE, TRUE, FALSE, FALSE, FALSE),
    bothCoreConnected = c(TRUE, FALSE, TRUE, FALSE, TRUE),
    bothIntraSupported = c(TRUE, TRUE, TRUE, TRUE, FALSE)
  )
  expect_identical(
    as.character(groups),
    c(
      "Distinct URP + intra-supported",
      "Distinct URP + intra-supported",
      "Core-connected + intra-supported",
      "Intra-supported",
      "No additional context"
    )
  )
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

test_that("PPI context accepts an alternate classifier", {
  data <- make_context_test_data() %>%
    dplyr::mutate(Experimental.score = .data$SVM.score) %>%
    dplyr::select(-"SVM.score")
  result <- annotatePPIContext(
    data,
    coreThreshold = 1.3,
    classifier = "Experimental.score",
    scalingFactor = 1
  )
  expect_equal(result$settings$classifier, "Experimental.score")
  expect_true("Experimental.score" %in% names(result$PPIs))
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

test_that("context target-decoy curves apply the decoy scaling factor", {
  curve <- contextTargetDecoyCurve(
    contextPEP = c(rep(0.01, 100), rep(0.02, 125)),
    decoyClass = c(
      rep("Target", 100), rep("Decoy", 100), rep("DoubleDecoy", 25)
    ),
    scalingFactor = 5
  )
  expect_equal(tail(curve$contextFDR, 1), 0.19)
  expect_equal(curve$contextQValue[[1]], 0)
})

make_context_calibration_data <- function() {
  groups <- c(
    "Distinct URP + intra-supported",
    "Core-connected + intra-supported",
    "Intra-supported",
    "No additional context"
  )
  purrr::map_dfr(groups, function(group) {
    target.score <- seq(0.5, 4, length.out = 30)
    decoy.score <- seq(-1, 1, length.out = 8)
    double.score <- seq(-1, 0, length.out = 2)
    score <- c(target.score, decoy.score, double.score)
    tibble::tibble(
      ppiID = paste0(group, "-PPI", seq_len(40)),
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

test_that("weighted context calibration is monotonic and reproducible", {
  ppis <- make_context_calibration_data()
  context <- structure(
    list(
      PPIs = ppis,
      URPs = tibble::tibble(
        ppiID = ppis$ppiID,
        xlinkedResPair = paste0("URP-", seq_len(nrow(ppis)))
      ),
      CSMs = tibble::tibble(
        ppiID = ppis$ppiID,
        Spectrum = paste0("scan-", seq_len(nrow(ppis)))
      ),
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
    "contextPEP", "contextFDR", "contextQValue", "selectionFrequency",
    "coreSupported", "contextSelected", "classificationStatus"
  ) %in% names(first$PPIs)))
  expect_false(any(c("contextQualified", "classified") %in%
                     names(first$PPIs)))
  expect_identical(
    first$PPIs$contextSelected,
    is.finite(first$PPIs$contextQValue) &
      first$PPIs$contextQValue <= first$fdr$requested
  )
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
  expect_equal(first$settings$classification,
               "direct-target-decoy-q-value")
  expect_false(first$settings$protectedCore)
  expect_equal(first$fdr$requested, 0.05)
  expect_true(is.numeric(first$fdr$estimated))
  expect_equal(calculateFDR(first), first$fdr$estimated)
  expect_equal(countDecoys(first), first$fdr$counts)

  classified <- getPPIs(first, view = "classified")
  clean <- getPPIs(first, view = "clean")
  expect_true(all(classified$contextSelected))
  expect_true(all(clean$Decoy == "Target"))
  evidence <- getPPIEvidence(first, ppiID = first$PPIs$ppiID[[1]])
  expect_equal(nrow(evidence$PPI), 1)
  expect_equal(nrow(evidence$URPs), 1)
  expect_equal(nrow(evidence$CSMs), 1)
  pair.evidence <- getPPIEvidence(
    first,
    proteinPair = as.character(first$PPIs$xlinkedProtPair[[1]])
  )
  expect_equal(pair.evidence$PPI$ppiID, evidence$PPI$ppiID)
  expect_output(print(first), "contextual PPI results")

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

test_that("target-decoy q-values accept or reject complete PEP plateaus", {
  pep <- c(rep(0.001, 10), rep(0.03, 103))
  decoy.class <- c(rep("Target", 110), rep("Decoy", 3))
  q.value <- contextPEPQValues(pep, decoy.class)

  expect_equal(contextPEPThreshold(pep, decoy.class, 0.02), 0.001)
  expect_length(unique(q.value[pep == 0.03]), 1)
  expect_gt(unique(q.value[pep == 0.03]), 0.02)
})

test_that("direct context classification can re-evaluate an ordinary core", {
  data <- tibble::tibble(
    SVM.score = c(rep(2, 100), 1, rep(0, 3)),
    contextPEP = c(rep(0.001, 100), rep(0.5, 4)),
    Decoy = factor(
      c(rep("Target", 101), rep("Decoy", 3)),
      levels = c("DoubleDecoy", "Decoy", "Target")
    ),
    xlinkClass = "interProtein",
    coreSupported = c(rep(FALSE, 100), TRUE, rep(FALSE, 3))
  )
  result <- classifyPPIContextTable(
    data,
    bootstrap = NULL,
    targetER = 0.02,
    stabilityThreshold = 0.9,
    scoreColumn = "SVM.score",
    scalingFactor = 1
  )
  expect_false(result$PPIs$contextSelected[[101]])
  expect_equal(
    as.character(result$PPIs$classificationStatus[[101]]),
    "Ordinary core not context-selected"
  )
  expect_equal(result$threshold, 0.001)
})

test_that("weighted decreasing PAVA pools local score reversals", {
  fitted <- weightedDecreasingPAVA(
    c(0.4, 0.3, 0.35, 0.1),
    rep(1, 4)
  )
  expect_equal(fitted, c(0.4, 0.325, 0.325, 0.1))
  expect_true(all(diff(fitted) <= 0))
})
