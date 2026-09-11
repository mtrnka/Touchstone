make_results_test_training <- function() {
  csms <- tibble::tibble(
    xlinkedResPair = c("R1", "R1", "R2", "R3"),
    xlinkedPepPair = c("PEP1", "PEP1", "PEP2", "PEP3"),
    xlinkedProtPair = c("P1::P2", "P1::P2", "P1::P2", "P3::P3"),
    xlinkedModulPair = c("M1::M2", "M1::M2", "M1::M2", "M3::M3"),
    Acc.1 = c("P1", "P1", "r1_P1", "P3"),
    Acc.2 = c("P2", "P2", "P2", "P3"),
    Protein.1 = c("P1", "P1", "P1", "P3"),
    Protein.2 = c("P2", "P2", "P2", "P3"),
    XLink.AA.1 = c(10, 10, 20, 30),
    XLink.AA.2 = c(20, 20, 30, 40),
    DB.Peptide.1 = c("AAAA", "AAAA", "BBBB", "CCCC"),
    DB.Peptide.2 = c("DDDD", "DDDD", "EEEE", "CCCC"),
    Score.Diff = c(20, 18, 16, 14),
    SVM.score = c(4, 3, 3, 2),
    xlinkClass = c("interProtein", "interProtein", "interProtein", "intraProtein"),
    Decoy = c("Target", "Target", "Decoy", "Target")
  )
  radial.csms <- dplyr::mutate(csms, SVM.score = .data$SVM.score + 10)
  linear <- list(
    CSMs = csms,
    URPs = bestResPair(csms),
    thresh = list(interThresh = 2, intraThresh = 1),
    kernel = "linear",
    cost = 0.01,
    gamma = NA_real_
  )
  radial <- list(
    CSMs = radial.csms,
    URPs = bestResPair(radial.csms),
    thresh = list(interThresh = 12, intraThresh = 11),
    kernel = "radial",
    cost = 0.01,
    gamma = 0.001
  )

  structure(
    list(
      recommended = linear,
      recommendedRadial = radial,
      candidates = tibble::tibble(
        index = 1:2,
        recommended = c(TRUE, FALSE),
        recommendedRadial = c(FALSE, TRUE)
      ),
      models = list(linear, radial),
      settings = list(targetER = 0.01, scalingFactor = 5)
    ),
    class = "touchstone_training"
  )
}

test_that("prepared URP results reuse training output and classify it", {
  training <- make_results_test_training()
  observed <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    countDecoys = function(datTab, threshold, scalingFactor, ...) {
      observed$threshold <- threshold
      observed$scalingFactor <- scalingFactor
      tibble::tibble(xlinkClass = "interProtein", Target = 1, FDR = 0)
    },
    calculateFDR = function(datTab, threshold, scalingFactor, ...) {
      observed$fdr.data <- datTab
      observed$fdr.threshold <- threshold
      observed$fdr.scalingFactor <- scalingFactor
      0.012
    },
    .package = "touchstone"
  )

  result <- prepareCrosslinkResults(training, summarizationLevel = "urp")

  expect_s3_class(result, "touchstone_results")
  expect_identical(result$data, training$recommended$URPs)
  expect_identical(result$sourceCSMs, training$recommended$CSMs)
  expect_identical(result$stage, "prepared")
  expect_identical(result$thresholds, training$recommended$thresh)
  expect_identical(observed$threshold, training$recommended$thresh)
  expect_identical(observed$scalingFactor, 5)
  expect_identical(observed$fdr.data, result$data)
  expect_identical(observed$fdr.threshold, result$thresholds)
  expect_identical(observed$fdr.scalingFactor, 5)
  expect_identical(result$fdr, 0.012)
  expect_identical(result$summarizationLevel, "urp")
  expect_identical(result$model$kernel, "linear")
  expect_false(any(c("classified", "reported") %in% names(result)))
})

test_that("prepared results accept and record manual thresholds", {
  training <- make_results_test_training()
  manual <- list(interThresh = 2.5, intraThresh = 1.5)

  result <- prepareCrosslinkResults(
    training,
    summarizationLevel = "urp",
    thresholds = manual,
    scalingFactor = 1
  )

  expect_identical(result$thresholds, manual)
  expect_identical(result$settings$thresholdSource, "manual")
})

test_that("prepared results use existing summarization and threshold functions", {
  training <- make_results_test_training()
  expected <- bestProtPair(training$recommended$CSMs)
  observed <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    findSeparateThresholdsModelled = function(datTab, targetER,
                                               scalingFactor, plot,
                                               classifier) {
      observed$threshold.data <- datTab
      observed$targetER <- targetER
      observed$scalingFactor <- scalingFactor
      observed$plot <- plot
      observed$classifier <- classifier
      list(interThresh = 2.5, intraThresh = 1.5)
    },
    countDecoys = function(datTab, threshold, scalingFactor, ...) {
      observed$classified.data <- datTab
      observed$classified.threshold <- threshold
      observed$count.scalingFactor <- scalingFactor
      tibble::tibble(xlinkClass = "interProtein", Target = 1, FDR = 0)
    },
    .package = "touchstone"
  )

  result <- prepareCrosslinkResults(
    training,
    summarizationLevel = "protein-pair"
  )

  expect_identical(result$data, expected)
  expect_identical(observed$threshold.data, expected)
  expect_identical(observed$classified.data, expected)
  expect_identical(observed$classified.threshold, result$thresholds)
  expect_identical(observed$targetER, 0.01)
  expect_identical(observed$scalingFactor, 5)
  expect_identical(observed$count.scalingFactor, 5)
  expect_identical(observed$classifier, "SVM.score")
  expect_false(observed$plot)
})

test_that("prepared results can use the recommended radial fit", {
  training <- make_results_test_training()

  testthat::local_mocked_bindings(
    countDecoys = function(...) {
      tibble::tibble(xlinkClass = "interProtein", Target = 1, FDR = 0)
    },
    .package = "touchstone"
  )

  result <- prepareCrosslinkResults(
    training,
    summarizationLevel = "urp",
    model = "radial"
  )

  expect_identical(result$model$kernel, "radial")
  expect_identical(result$model$index, 2L)
  expect_identical(result$data, training$recommendedRadial$URPs)
})

test_that("prepared results accept scored CSM data directly", {
  training <- make_results_test_training()
  csms <- training$recommended$CSMs
  observed <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    findSeparateThresholdsModelled = function(datTab, ...) {
      list(interThresh = 2, intraThresh = 1)
    },
    countDecoys = function(datTab, threshold, scalingFactor, ...) {
      observed$called <- TRUE
      tibble::tibble(xlinkClass = "interProtein", Target = 1, FDR = 0)
    },
    .package = "touchstone"
  )

  result <- prepareCrosslinkResults(
    csms,
    summarizationLevel = "spectra",
    targetER = 0.01,
    scalingFactor = 5
  )

  expect_identical(result$summarizationLevel, "csm")
  expect_identical(result$data, csms)
  expect_true(observed$called)
})

test_that("unscored CSM data produce an actionable error", {
  unscored <- tibble::tibble(
    xlinkClass = "interProtein",
    Decoy = "Target"
  )

  expect_error(
    prepareCrosslinkResults(unscored),
    "no classifier column named 'SVM.score'"
  )
})

test_that("prepared results consistently use an alternate classifier", {
  training <- make_results_test_training()
  training$recommended$CSMs$experimental.score <- c(10, 1, 3, 2)
  observed <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    findSeparateThresholdsModelled = function(datTab, classifier, ...) {
      observed$threshold.classifier <- classifier
      list(interThresh = 2.5, intraThresh = 1.5)
    },
    countDecoys = function(datTab, threshold, classifier, ...) {
      observed$count.classifier <- classifier
      tibble::tibble(xlinkClass = "interProtein", Target = 1, FDR = 0)
    },
    calculateFDR = function(datTab, threshold, classifier, ...) {
      observed$fdr.classifier <- classifier
      0.01
    },
    .package = "touchstone"
  )

  result <- prepareCrosslinkResults(
    training,
    summarizationLevel = "urp",
    classifier = "experimental.score"
  )

  expect_identical(result$data$experimental.score, c(10, 3, 2))
  expect_identical(result$settings$classifier, "experimental.score")
  expect_identical(observed$threshold.classifier, "experimental.score")
  expect_identical(observed$count.classifier, "experimental.score")
  expect_identical(observed$fdr.classifier, "experimental.score")
})

test_that("classification recalculates support counts after thresholding", {
  training <- make_results_test_training()
  thresholds <- list(interThresh = 2.5, intraThresh = 2.5)
  prepared <- prepareCrosslinkResults(
    training,
    summarizationLevel = "urp",
    thresholds = thresholds,
    scalingFactor = 1
  )

  result <- classifyCrosslinkResults(prepared)

  expect_identical(result$stage, "classified")
  expect_identical(result$settings$thresholdSource, "manual")
  expect_equal(nrow(result$data), 2)
  expect_setequal(result$data$numCSM, c(2, 1))
  expect_false(any(c("wtCSM", "wtURP") %in% names(result$data)))
  expect_identical(result$sourceCSMs, prepared$sourceCSMs)
  expect_identical(result$polishingAudit$rule, "threshold")
  expect_identical(result$polishingAudit$removed, 1L)
})

test_that("classification applies named Prospector-style polishing rules", {
  training <- make_results_test_training()
  csms <- training$recommended$CSMs %>%
    dplyr::mutate(
      Len.Pep.1 = c(10, 10, 10, 4),
      Len.Pep.2 = c(10, 10, 10, 4),
      Sc.1 = c(10, 4, 10, 10),
      Sc.2 = 10,
      numProdIons.1 = c(5, 5, 2, 5),
      numProdIons.2 = 5,
      ladderLen.1 = c(3, 3, 3, 3),
      ladderLen.2 = c(3, 3, 1, 3)
    )
  prepared <- prepareCrosslinkResults(
    csms,
    summarizationLevel = "urp",
    thresholds = -100,
    targetER = 0.01,
    scalingFactor = 1
  )

  result <- classifyCrosslinkResults(
    prepared,
    polishing = list(
      minPepLen = 5,
      minPepScore = 5,
      minScoreDiff = 15,
      minIons = 3,
      minLadderCoverage = 0.25
    )
  )

  expect_identical(result$stage, "polished")
  expect_equal(nrow(result$data), 1)
  expect_identical(
    result$polishingAudit$rule,
    c(
      "threshold", "minPepLen", "minPepScore", "minScoreDiff", "minIons",
      "minLadderCoverage"
    )
  )
  expect_identical(
    result$settings$polishing$minLadderCoverage,
    0.25
  )
})

test_that("classification skips and reports unavailable polishing evidence", {
  prepared <- prepareCrosslinkResults(
    make_results_test_training(),
    thresholds = -100,
    scalingFactor = 1
  )

  observed.warnings <- character()
  result <- withCallingHandlers(
    classifyCrosslinkResults(
      prepared,
      polishing = list(
        minPepScore = 5,
        minIons = 3,
        minLadderCoverage = 0.25
      )
    ),
    warning = function(warning) {
      observed.warnings <<- c(observed.warnings, conditionMessage(warning))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(observed.warnings, 3)
  expect_true(all(grepl("was skipped", observed.warnings)))
  expect_identical(result$stage, "classified")
  expect_identical(
    result$polishingAudit$rule,
    c("threshold", "minPepScore", "minIons", "minLadderCoverage")
  )
  expect_identical(
    result$polishingAudit$applied,
    c(TRUE, FALSE, FALSE, FALSE)
  )
  expect_match(result$polishingAudit$reason[[2]], "Sc.1, Sc.2")
  expect_match(result$polishingAudit$reason[[3]], "numProdIons.1")
  expect_match(result$polishingAudit$reason[[4]], "ladderLen.1")
  expect_error(
    classifyCrosslinkResults(prepared, polishing = list(minProductIons = 3)),
    "Unknown polishing option"
  )
})
