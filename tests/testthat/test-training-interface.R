test_that("linear kernels are the default tuning grid", {
  default.grid <- touchstone:::makeSVMParameterGrid(
    cost_values = c(1, 5),
    gamma_values = c(0.01, 0.1)
  )
  experimental.grid <- touchstone:::makeSVMParameterGrid(
    cost_values = c(1, 5),
    gamma_values = c(0.01, 0.1),
    kernels = c("linear", "radial")
  )

  expect_identical(default.grid$kernel, c("linear", "linear"))
  expect_true(all(is.na(default.grid$gamma)))
  expect_equal(sum(experimental.grid$kernel == "linear"), 2)
  expect_equal(sum(experimental.grid$kernel == "radial"), 4)
})

make_complexity_data <- function(n) {
  data.frame(
    Acc.1 = paste0("P", seq_len(n)),
    Acc.2 = "P1",
    Decoy = "Target"
  )
}

test_that("automatic complexity profiles use reported protein boundaries", {
  small <- touchstone:::resolveDatasetComplexity(make_complexity_data(20))
  medium.low <- touchstone:::resolveDatasetComplexity(make_complexity_data(21))
  medium.high <- touchstone:::resolveDatasetComplexity(make_complexity_data(200))
  large <- touchstone:::resolveDatasetComplexity(make_complexity_data(201))

  expect_identical(small$selected, "small")
  expect_identical(medium.low$selected, "medium")
  expect_identical(medium.high$selected, "medium")
  expect_identical(large$selected, "large")
  expect_identical(large$proteinCount, 201L)
  expect_identical(large$breaks, c(smallMax = 20L, mediumMax = 200L))
})

test_that("complexity profile and boundaries can be overridden", {
  forced <- touchstone:::resolveDatasetComplexity(
    make_complexity_data(500),
    complexity = "small"
  )
  custom <- touchstone:::resolveDatasetComplexity(
    make_complexity_data(51),
    complexityBreaks = c(50, 500)
  )

  expect_identical(forced$selected, "small")
  expect_identical(forced$requested, "small")
  expect_identical(custom$selected, "medium")
})

test_that("complexity feature profiles add higher-order features gradually", {
  columns <- c(
    "Score.Diff", "percMatched", "z", "wtCSM", "wtURP", "xlinkClass",
    "Perc.Bond.Cleavage.1", "Perc.Bond.Cleavage.2"
  )
  datTab <- as.data.frame(stats::setNames(rep(list(numeric()), length(columns)),
                                      columns))

  small <- touchstone:::complexityFeatureProfile("small", datTab)
  medium <- touchstone:::complexityFeatureProfile("medium", datTab)
  large <- touchstone:::complexityFeatureProfile("large", datTab)

  expect_false(any(c("xlinkClass", "wtURP") %in% small))
  expect_true("xlinkClass" %in% medium)
  expect_false("wtURP" %in% medium)
  expect_true(all(c("xlinkClass", "wtURP") %in% large))
  expect_true(all(c("Perc.Bond.Cleavage.1", "Perc.Bond.Cleavage.2") %in% small))
})

make_training_complexity_data <- function(n) {
  tibble::tibble(
    Acc.1 = paste0("P", seq_len(n)),
    Acc.2 = "P1",
    Decoy = "Target",
    Score.Diff = seq_len(n),
    percMatched = 0.5,
    z = 3,
    wtCSM = 1,
    wtURP = 1,
    xlinkClass = "interProtein"
  )
}

mock_tuning_result <- function(datTab, params) {
  list(list(
    CSMs = datTab,
    URPs = datTab,
    thresh = list(intraThresh = 0, interThresh = 0),
    errorTable = tibble::tibble(fdr.inter = 0.01, inter = 1),
    interInt = 1,
    interHits = 1,
    corScore = 80,
    cost = 1,
    gamma = NA_real_,
    kernel = "linear",
    params = params
  ))
}

test_that("training records automatic complexity and selected features", {
  input <- make_training_complexity_data(10)
  observed <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    tuneSVM = function(datTab, params, ...) {
      observed$params <- params
      mock_tuning_result(datTab, params)
    },
    .package = "touchstone"
  )

  training <- trainCrosslinkScore(input)

  expect_identical(training$settings$complexity$selected, "small")
  expect_identical(training$settings$complexity$proteinCount, 10L)
  expect_identical(training$settings$featureSource, "complexity-profile")
  expect_identical(training$settings$features, observed$params)
  expect_false(any(c("xlinkClass", "wtURP") %in% observed$params))
  expect_false(training$prefilter$applied)
  expect_null(training$prefilter$selectedScoreDiff)
})

test_that("large profile runs and records Score.Diff prefilter selection", {
  input <- make_training_complexity_data(201)
  observed <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    chooseScoreDiffPrefilter = function(datTab, ...) {
      observed$prefilter.called <- TRUE
      list(
        datTab = datTab[seq_len(100), , drop = FALSE],
        preFilter.summary = tibble::tibble(sd.thresh = c(0, 10)),
        bestPreFilter = 10
      )
    },
    tuneSVM = function(datTab, params, ...) {
      observed$tuning.rows <- nrow(datTab)
      observed$params <- params
      mock_tuning_result(datTab, params)
    },
    .package = "touchstone"
  )

  training <- trainCrosslinkScore(input)

  expect_true(observed$prefilter.called)
  expect_identical(training$settings$complexity$selected, "large")
  expect_true(all(c("xlinkClass", "wtURP") %in% observed$params))
  expect_true(training$prefilter$applied)
  expect_identical(training$prefilter$selectedScoreDiff, 10)
  expect_identical(training$prefilter$rowsBefore, 201L)
  expect_identical(training$prefilter$rowsAfter, 100L)
  expect_identical(observed$tuning.rows, 100L)
})

test_that("explicit params override complexity feature selection", {
  input <- make_training_complexity_data(250)
  observed <- new.env(parent = emptyenv())
  requested.params <- c("Score.Diff", "percMatched", "massError")

  testthat::local_mocked_bindings(
    chooseScoreDiffPrefilter = function(datTab, ...) {
      list(
        datTab = datTab,
        preFilter.summary = tibble::tibble(sd.thresh = 0),
        bestPreFilter = 0
      )
    },
    tuneSVM = function(datTab, params, ...) {
      observed$params <- params
      mock_tuning_result(datTab, params)
    },
    .package = "touchstone"
  )

  training <- trainCrosslinkScore(input, params = requested.params)

  expect_identical(observed$params, requested.params)
  expect_identical(training$settings$featureSource, "user")
})

test_that("candidate selection prefers correlation within near-best recovery", {
  candidates <- tibble::tibble(
    index = 1:4,
    kernel = c("linear", "linear", "radial", "radial"),
    cost = c(1, 5, 1, 5),
    gamma = c(NA_real_, NA_real_, 0.1, 0.01),
    interInt = c(40, 45, 55, 50),
    interHits = c(40, 45, 55, 50),
    scoreCorrelation = c(0.8, 0.8, 0.7, 0.75),
    eligible = c(TRUE, TRUE, TRUE, TRUE)
  )

  selection <- touchstone:::selectSVMCandidates(candidates)

  expect_identical(selection$recommended, 2L)
  expect_identical(selection$recommendedRadial, 4L)
  expect_true(selection$candidates$recommended[2])
  expect_true(selection$candidates$recommendedRadial[4])
  expect_false(selection$candidates$recommended[3])
  expect_identical(
    selection$candidates$nearBestRecovery,
    c(FALSE, TRUE, TRUE, TRUE)
  )
  expect_identical(selection$candidates$index, 1:4)
})

test_that("radial recommendation remains separate when linear is ineligible", {
  candidates <- tibble::tibble(
    index = 1:3,
    kernel = c("linear", "radial", "radial"),
    cost = c(1, 1, 5),
    gamma = c(NA_real_, 0.1, 0.01),
    interInt = c(60, 40, 45),
    interHits = c(60, 40, 45),
    scoreCorrelation = c(0.9, 0.7, 0.75),
    eligible = c(FALSE, TRUE, TRUE)
  )

  selection <- touchstone:::selectSVMCandidates(candidates)

  expect_true(is.na(selection$recommended))
  expect_identical(selection$recommendedRadial, 3L)
})

test_that("recovery fraction prevents large recovery sacrifices", {
  candidates <- tibble::tibble(
    index = 1:3,
    kernel = rep("radial", 3),
    cost = c(0.001, 0.01, 0.1),
    gamma = c(0.001, 0.001, 0.001),
    interInt = c(70, 60, 50),
    interHits = c(100, 91, 89),
    scoreCorrelation = c(0.7, 0.8, 0.99),
    eligible = rep(TRUE, 3)
  )

  selection <- touchstone:::selectSVMCandidates(
    candidates,
    recoveryFraction = 0.9
  )

  expect_true(is.na(selection$recommended))
  expect_identical(selection$recommendedRadial, 2L)
  expect_identical(
    selection$candidates$nearBestRecovery,
    c(TRUE, TRUE, FALSE)
  )
})

test_that("tuning passes the explicit decoy scaling factor to its error table", {
  mock.csms <- tibble::tibble(
    Decoy = rep("Target", 6),
    xlinkClass = c(rep("interProtein", 4), rep("intraProtein", 2)),
    Score.Diff = 6:1,
    SVM.score = 6:1
  )
  observed <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    buildSVM = function(...) mock.csms,
    bestResPair = function(datTab, ...) datTab,
    findSeparateThresholdsModelled = function(...) {
      list(interThresh = -Inf, intraThresh = -Inf)
    },
    classifyDataset = function(datTab, ...) datTab,
    removeDecoys = function(datTab) datTab,
    generateErrorTable.sep = function(datTab, scalingFactor, ...) {
      observed$scalingFactor <- scalingFactor
      tibble::tibble(fdr.inter = 0.02, inter = 4)
    },
    .package = "touchstone"
  )

  tuneSVM.helper(
    mock.csms,
    params = "Score.Diff",
    scalingFactor = 5
  )

  expect_identical(observed$scalingFactor, 5)
})

test_that("tuning evaluation respects an alternate score name", {
  mock.csms <- tibble::tibble(
    Decoy = rep("Target", 6),
    xlinkClass = c(rep("interProtein", 4), rep("intraProtein", 2)),
    Score.Diff = 6:1,
    experimental.score = 6:1
  )
  observed <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    buildSVM = function(...) mock.csms,
    bestResPair = function(datTab, classifier, ...) {
      observed$summarization.classifier <- classifier
      datTab
    },
    findSeparateThresholdsModelled = function(datTab, classifier, ...) {
      observed$threshold.classifier <- classifier
      list(interThresh = -Inf, intraThresh = -Inf)
    },
    classifyDataset = function(datTab, classifier, ...) {
      observed$classification.classifier <- classifier
      datTab
    },
    removeDecoys = function(datTab) datTab,
    generateErrorTable.sep = function(datTab, classifier, ...) {
      observed$error.classifier <- classifier
      tibble::tibble(fdr.inter = 0.02, inter = 4)
    },
    .package = "touchstone"
  )

  result <- tuneSVM.helper(
    mock.csms,
    params = "Score.Diff",
    scoreName = "experimental.score"
  )

  expect_identical(observed$summarization.classifier, "experimental.score")
  expect_identical(observed$threshold.classifier, "experimental.score")
  expect_identical(observed$classification.classifier, "experimental.score")
  expect_identical(observed$error.classifier, "experimental.score")
  expect_equal(result$corScore, 100)
})

test_that("FDR-versus-hit plots retain weak candidates in separate facets", {
  make_model <- function(kernel, gamma, cost, hits) {
    list(
      kernel = kernel,
      gamma = gamma,
      cost = cost,
      errorTable = data.frame(
        fdr.inter = c(0, 0.01, 0.03, 0.05),
        inter = hits,
        fdr.intra = c(0, 0.01, 0.03, 0.05),
        intra = hits + 10
      )
    )
  }

  models <- list(
    make_model("linear", NA_real_, 1, c(20, 18, 15, 12)),
    make_model("radial", 0.1, 1, c(2, 4, 1, 3))
  )
  training <- structure(
    list(models = models, settings = list(targetER = 0.01)),
    class = "touchstone_training"
  )

  plot <- plotSVMTuning(training)

  expect_s3_class(plot, "ggplot")
  expect_equal(nrow(plot$data), 8)
  expect_setequal(
    as.character(unique(plot$data$model)),
    c("linear", "radial (gamma = 0.1)")
  )
  radial.frontier <- plot$data$hits[as.character(plot$data$model) != "linear"]
  expect_true(all(diff(radial.frontier) >= 0))
  expect_equal(radial.frontier, c(2, 4, 4, 4))
  expect_length(plot$layers, 3)

  plot.without.raw <- plotSVMTuning(training, showRaw = FALSE)
  expect_length(plot.without.raw$layers, 2)
})
