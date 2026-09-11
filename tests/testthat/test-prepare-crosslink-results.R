make_results_test_training <- function() {
  csms <- tibble::tibble(
    xlinkedResPair = c("R1", "R1", "R2", "R3"),
    xlinkedPepPair = c("PEP1", "PEP1", "PEP2", "PEP3"),
    xlinkedProtPair = c("P1::P2", "P1::P2", "P1::P2", "P3::P3"),
    xlinkedModulPair = c("M1::M2", "M1::M2", "M1::M2", "M3::M3"),
    SVM.score = c(1, 4, 3, 2),
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
