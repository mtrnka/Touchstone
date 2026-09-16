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
    Acc.2 = paste0("P", seq_len(n)),
    Decoy = "Target",
    Score.Diff = 11,
    numCSM = 2,
    xlinkClass = "intraProtein"
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
  expect_identical(large$rawProteinCount, 201L)
  expect_identical(large$breaks, c(smallMax = 20L, mediumMax = 200L))
})

test_that("automatic complexity does not require repeated URP observations", {
  input <- make_complexity_data(201)
  input$numCSM <- NULL

  result <- touchstone:::resolveDatasetComplexity(input)

  expect_identical(result$proteinCount, 201L)
  expect_identical(result$selected, "large")
  expect_false("numCSMGreaterThan" %in% names(result$evidenceCriteria))
})

test_that("automatic complexity counts both proteins in inter-protein CSMs", {
  input <- data.frame(
    Acc.1 = paste0("P", seq_len(101)),
    Acc.2 = paste0("P", 102:202),
    Decoy = "Target",
    Score.Diff = 11,
    xlinkClass = "interProtein"
  )

  result <- touchstone:::resolveDatasetComplexity(input)

  expect_identical(result$proteinCount, 202L)
  expect_identical(result$intraProteinCount, 0L)
  expect_identical(result$highScoringCSMCount, 101L)
  expect_identical(result$highScoringIntraCSMCount, 0L)
  expect_identical(result$selected, "large")
})

test_that("automatic complexity includes annotated inter-protein classes", {
  input <- data.frame(
    Acc.1 = c("P1", "P3"),
    Acc.2 = c("P2", "P4"),
    Decoy = "Target",
    Score.Diff = 11,
    xlinkClass = c("interProtein, homomeric", "interProtein, heteromeric")
  )

  result <- touchstone:::resolveDatasetComplexity(
    input,
    complexity = "small"
  )

  expect_identical(result$proteinCount, 4L)
})

test_that("protein dominance counts participation once per CSM", {
  input <- data.frame(
    Acc.1 = c("P1", "P1"),
    Acc.2 = c("P1", "P2"),
    Decoy = "Target",
    Score.Diff = 11,
    xlinkClass = c("intraProtein", "interProtein")
  )

  result <- touchstone:::resolveDatasetComplexity(
    input,
    complexity = "small"
  )

  expect_identical(result$proteinCount, 2L)
  expect_identical(result$intraProteinCount, 1L)
  expect_equal(result$dominantProteinRatio, 2)
})

test_that("automatic complexity ignores unsupported background accessions", {
  input <- make_complexity_data(250)
  input$Score.Diff[-1] <- 5

  result <- touchstone:::resolveDatasetComplexity(input)

  expect_identical(result$rawProteinCount, 250L)
  expect_identical(result$proteinCount, 1L)
  expect_identical(result$highScoringCSMCount, 1L)
  expect_identical(result$selected, "small")
})

test_that("dominant plausible CSM evidence keeps a background-rich system small", {
  dominant <- make_complexity_data(30)
  dominant <- dplyr::bind_rows(
    dominant,
    dominant[rep(1, 100), , drop = FALSE]
  )

  result <- touchstone:::resolveDatasetComplexity(dominant)

  expect_identical(result$proteinCount, 30L)
  expect_true(result$dominanceOverride)
  expect_gt(result$dominantProteinRatio, 100)
  expect_identical(result$selected, "small")
})

test_that("explicit complexity overrides a dominance downgrade", {
  dominant <- make_complexity_data(30)
  dominant <- dplyr::bind_rows(
    dominant,
    dominant[rep(1, 100), , drop = FALSE]
  )

  result <- touchstone:::resolveDatasetComplexity(
    dominant,
    complexity = "medium"
  )

  expect_false(result$dominanceOverride)
  expect_identical(result$selected, "medium")
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
    "Score.Diff", "percMatched", "z", "CSMsupport", "URPsupport", "xlinkClass",
    "Perc.Bond.Cleavage.1", "Perc.Bond.Cleavage.2"
  )
  datTab <- as.data.frame(stats::setNames(rep(list(numeric()), length(columns)),
                                      columns))

  small <- touchstone:::complexityFeatureProfile("small", datTab)
  medium <- touchstone:::complexityFeatureProfile("medium", datTab)
  large <- touchstone:::complexityFeatureProfile("large", datTab)

  expect_true("CSMsupport" %in% small)
  expect_false(any(c("xlinkClass", "URPsupport") %in% small))
  expect_true("xlinkClass" %in% medium)
  expect_false("URPsupport" %in% medium)
  expect_true(all(c("xlinkClass", "URPsupport") %in% large))
  expect_true(all(c("Perc.Bond.Cleavage.1", "Perc.Bond.Cleavage.2") %in% small))
})

make_training_complexity_data <- function(n) {
  tibble::tibble(
    Acc.1 = paste0("P", seq_len(n)),
    Acc.2 = paste0("P", seq_len(n)),
    Decoy = "Target",
    Decoy2 = factor(
      rep(c("Target", "Decoy"), length.out = n),
      levels = c("Decoy", "Target")
    ),
    Score.Diff = 11 + seq_len(n),
    numCSM = 2,
    percMatched = 0.5,
    z = 3,
    CSMsupport = 1,
    URPsupport = 1,
    xlinkClass = "intraProtein"
  )
}

mock_tuning_result <- function(datTab, params) {
  datTab$SVM.score <- datTab$Score.Diff
  list(list(
    CSMs = datTab,
    URPs = datTab,
    thresh = list(intraThresh = 0, interThresh = 0),
    errorTable = tibble::tibble(
      fdr.inter = 0.01, inter = 0,
      fdr.intra = 0.01, intra = 1
    ),
    interInt = 1,
    interHits = 0,
    intraHits = 1,
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
  expect_identical(training$settings$complexity$intraProteinCount, 10L)
  expect_identical(training$settings$complexity$rawProteinCount, 10L)
  expect_identical(training$settings$featureSource, "complexity-profile")
  expect_identical(training$settings$features, observed$params)
  expect_false(any(c("xlinkClass", "URPsupport") %in% observed$params))
  expect_false(training$prefilter$applied)
  expect_null(training$prefilter$selectedScoreDiff)
  expect_true(all(c(
    "complexity", "requestedComplexity", "featureSource", "featureCount",
    "features", "scoreDiffPrefilterEvaluated", "scoreDiffPrefilter",
    "rowsBeforePrefilter", "rowsAfterPrefilter", "interThreshold",
    "intraThreshold", "interHits", "intraHits", "totalHits",
    "achievedFDR", "interCorrelation", "intraCorrelation",
    "interTailCorrelation", "intraTailCorrelation",
    "interCorrelationN", "intraCorrelationN", "worstClassCorrelation",
    "correlationAvailable", "minimumCorrelationRequired",
    "correlationCredible", "targetFDR",
    "scalingFactor", "validation", "selectionBasis", "recoveryHits",
    "bestRecoveryHits", "recoveryRelativeToBest", "nearBestRecovery",
    "selectionReason"
  ) %in% names(training$candidates)))
  expect_false("interInt" %in% names(training$candidates))
  expect_identical(training$models[[1]]$interInt, 1)
  expect_identical(training$candidates$complexity, "small")
  expect_identical(training$candidates$featureSource, "complexity-profile")
  expect_identical(training$candidates$featureCount, length(observed$params))
  expect_identical(
    training$candidates$features,
    paste(observed$params, collapse = ", ")
  )
  expect_false(training$candidates$scoreDiffPrefilterEvaluated)
  expect_true(is.na(training$candidates$scoreDiffPrefilter))
  expect_identical(training$candidates$interThreshold, 0)
  expect_identical(training$candidates$intraThreshold, 0)
  expect_identical(training$candidates$totalHits, 1)
  expect_identical(training$candidates$selectionBasis, "intra")
  expect_identical(
    training$candidates$selectionReason,
    "Recommended linear candidate"
  )
})

test_that("training stores compact candidates and materializes requested fits", {
  input <- make_training_complexity_data(10)
  score <- input$Score.Diff
  diagnostics <- touchstone:::summarizeScoreBehavior(
    transform(input, SVM.score = score)
  )
  compact <- list(
    score = score,
    scoreDiagnostics = diagnostics,
    achievedFDR = 0.01,
    thresh = list(intraThresh = 0, interThresh = Inf),
    errorTable = tibble::tibble(
      fdr.inter = 0, inter = 0,
      fdr.intra = 0.01, intra = 10
    ),
    interInt = NaN,
    interHits = 0,
    intraHits = 10,
    corScore = 100,
    cost = 0.001,
    gamma = NA_real_,
    kernel = "linear",
    params = "Score.Diff"
  )

  testthat::local_mocked_bindings(
    tuneSVM = function(...) list(compact),
    bestResPair = function(datTab, ...) datTab,
    .package = "touchstone"
  )

  training <- trainCrosslinkScore(input, params = "Score.Diff")

  expect_null(training$models[[1]]$CSMs)
  expect_null(training$models[[1]]$URPs)
  expect_identical(training$sourceCSMs, input)
  expect_equal(training$recommended$CSMs$SVM.score, score)
  resolved <- touchstone:::resolveCrosslinkFit(training, 1)
  expect_equal(resolved$fit$CSMs$SVM.score, score)
})

test_that("score diagnostics report concise within-class rank correlations", {
  monotonic <- data.frame(
    SVM.score = seq_len(100),
    Score.Diff = seq_len(100),
    Decoy = "Target",
    xlinkClass = "interProtein",
    CSMsupport = seq_len(100) / 100
  )
  reversing <- monotonic
  reversing$Score.Diff <- c(seq_len(50), 50:1)

  monotonic.result <- touchstone:::summarizeScoreBehavior(monotonic)
  reversing.result <- touchstone:::summarizeScoreBehavior(reversing)

  expect_equal(monotonic.result$interCorrelation, 1)
  expect_equal(monotonic.result$interTailCorrelation, 1)
  expect_equal(monotonic.result$interCorrelationN, 100)
  expect_true(is.na(monotonic.result$intraCorrelation))
  expect_equal(monotonic.result$worstClassCorrelation, 1)
  expect_true(monotonic.result$correlationAvailable)
  expect_lt(reversing.result$interCorrelation, 1)
})

test_that("upper-tail correlation detects a fold hidden by full correlation", {
  input <- data.frame(
    SVM.score = c(1:60, 100:61),
    Score.Diff = 1:100,
    Decoy = "Target",
    xlinkClass = "intraProtein"
  )

  result <- touchstone:::summarizeScoreBehavior(input)

  expect_gt(result$intraCorrelation, 0.8)
  expect_lte(result$intraTailCorrelation, 0)
  expect_true(result$correlationAvailable)
})

test_that("score diagnostics report intra- and inter-protein behavior separately", {
  input <- data.frame(
    SVM.score = rep(seq_len(20), 2),
    Score.Diff = c(seq_len(20), rev(seq_len(20))),
    Decoy = "Target",
    xlinkClass = rep(c("interProtein, heteromeric", "intraProtein"), each = 20)
  )

  result <- touchstone:::summarizeScoreBehavior(input)

  expect_equal(result$interCorrelation, 1)
  expect_equal(result$intraCorrelation, -1)
  expect_equal(result$interTailCorrelation, 1)
  expect_equal(result$intraTailCorrelation, -1)
  expect_equal(result$interCorrelationN, 20)
  expect_equal(result$intraCorrelationN, 20)
  expect_equal(result$worstClassCorrelation, -1)
  expect_true(result$correlationAvailable)
})

test_that("poor within-class correlation is excluded before recovery ranking", {
  input <- data.frame(
    Acc.1 = "P1",
    Acc.2 = "P1",
    Decoy = "Target",
    Score.Diff = rep(seq_len(20), 2),
    xlinkClass = rep(c("interProtein", "intraProtein"), each = 20)
  )
  make.candidate <- function(score, hits, cost) {
    scored <- input
    scored$SVM.score <- score
    list(
      CSMs = scored,
      URPs = scored,
      thresh = list(intraThresh = 0, interThresh = 0),
      errorTable = tibble::tibble(
        fdr.inter = 0.01, inter = hits,
        fdr.intra = 0.01, intra = 20
      ),
      interInt = hits,
      interHits = hits,
      intraHits = 20,
      corScore = 80,
      cost = cost,
      gamma = NA_real_,
      kernel = "linear",
      params = "Score.Diff"
    )
  }
  badly.fitted <- make.candidate(
    c(seq_len(20), rev(seq_len(20))),
    hits = 200,
    cost = 0.1
  )
  credible <- make.candidate(input$Score.Diff, hits = 100, cost = 1)

  testthat::local_mocked_bindings(
    tuneSVM = function(...) list(badly.fitted, credible),
    .package = "touchstone"
  )

  training <- trainCrosslinkScore(input, params = "Score.Diff")

  expect_false(training$candidates$eligible[1])
  expect_identical(training$candidates$bestRecoveryHits, c(100, 100))
  expect_identical(training$candidates$recommended, c(FALSE, TRUE))
})

test_that("radial candidates must meet the stronger correlation minimum", {
  input <- data.frame(
    Acc.1 = "P1",
    Acc.2 = "P1",
    Decoy = "Target",
    Score.Diff = 1:100,
    xlinkClass = "intraProtein"
  )
  folded.score <- c(1:65, 100:66)
  make.candidate <- function(kernel, gamma) {
    scored <- input
    scored$SVM.score <- folded.score
    list(
      CSMs = scored,
      URPs = scored,
      thresh = list(intraThresh = 0, interThresh = Inf),
      errorTable = tibble::tibble(
        fdr.inter = 0, inter = 0,
        fdr.intra = 0.01, intra = 100
      ),
      interInt = NaN,
      interHits = 0,
      intraHits = 100,
      corScore = 50,
      cost = 0.1,
      gamma = gamma,
      kernel = kernel,
      params = "Score.Diff"
    )
  }

  testthat::local_mocked_bindings(
    tuneSVM = function(...) list(
      make.candidate("linear", NA_real_),
      make.candidate("radial", 0.01)
    ),
    .package = "touchstone"
  )

  training <- trainCrosslinkScore(
    input,
    params = "Score.Diff",
    kernels = c("linear", "radial")
  )

  expect_gt(training$candidates$worstClassCorrelation[1], 0.2)
  expect_lt(training$candidates$worstClassCorrelation[2], 0.5)
  expect_identical(training$candidates$minimumCorrelationRequired, c(0.2, 0.5))
  expect_identical(training$candidates$eligible, c(TRUE, FALSE))
  expect_null(training$recommendedRadial)
})

test_that("large profile runs and records Score.Diff prefilter selection", {
  input <- make_training_complexity_data(501)
  observed <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    chooseScoreDiffPrefilter = function(datTab, ...) {
      observed$prefilter.called <- TRUE
      list(
        datTab = datTab[datTab$Score.Diff >= 413, , drop = FALSE],
        preFilter.summary = tibble::tibble(sd.thresh = c(0, 413)),
        bestPreFilter = 413
      )
    },
    tuneSVM = function(datTab, params, trainingRows, ...) {
      observed$scored.rows <- nrow(datTab)
      observed$training.rows <- sum(trainingRows)
      observed$params <- params
      mock_tuning_result(datTab, params)
    },
    .package = "touchstone"
  )

  training <- trainCrosslinkScore(input)

  expect_true(observed$prefilter.called)
  expect_identical(training$settings$complexity$selected, "large")
  expect_true(all(c("xlinkClass", "URPsupport") %in% observed$params))
  expect_true(training$prefilter$applied)
  expect_identical(training$prefilter$selectedScoreDiff, 413)
  expect_identical(training$prefilter$rowsBefore, 501L)
  expect_identical(training$prefilter$rowsAfter, 100L)
  expect_identical(training$prefilter$rowsScored, 501L)
  expect_true(training$prefilter$trainingOnly)
  expect_identical(observed$scored.rows, 501L)
  expect_identical(observed$training.rows, 100L)
  expect_true(training$candidates$scoreDiffPrefilterEvaluated)
  expect_identical(training$candidates$scoreDiffPrefilter, 413)
  expect_identical(training$candidates$rowsBeforePrefilter, 501L)
  expect_identical(training$candidates$rowsAfterPrefilter, 100L)
  expect_identical(training$candidates$rowsScored, 501L)
})

test_that("fixed Score.Diff prefilter bypasses automatic selection", {
  input <- make_training_complexity_data(100)
  observed <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    chooseScoreDiffPrefilter = function(...) {
      stop("automatic prefilter selection should not be called")
    },
    tuneSVM = function(datTab, params, trainingRows, ...) {
      observed$trainingRows <- trainingRows
      mock_tuning_result(datTab, params)
    },
    .package = "touchstone"
  )

  training <- trainCrosslinkScore(
    input,
    params = "Score.Diff",
    scoreDiffPrefilter = 61
  )

  expect_identical(observed$trainingRows, input$Score.Diff >= 61)
  expect_identical(training$prefilter$selectedScoreDiff, 61)
  expect_identical(training$prefilter$assessment$reason, "fixed by the user")
  expect_identical(training$settings$fixedScoreDiffPrefilter, 61)
})

test_that("prefilter assessment is independent of the feature profile", {
  input <- make_training_complexity_data(100)
  input <- dplyr::bind_rows(replicate(6, input, simplify = FALSE))
  observed <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    chooseScoreDiffPrefilter = function(datTab, ...) {
      observed$called <- TRUE
      list(
        datTab = datTab,
        preFilter.summary = tibble::tibble(sd.thresh = 10),
        bestPreFilter = 10
      )
    },
    tuneSVM = function(datTab, params, trainingRows, ...) {
      observed$training.rows <- sum(trainingRows)
      mock_tuning_result(datTab, params)
    },
    .package = "touchstone"
  )

  training <- trainCrosslinkScore(input)

  expect_identical(training$settings$complexity$selected, "medium")
  expect_true(observed$called)
  expect_true(training$prefilter$applied)
  expect_identical(observed$training.rows, 600L)
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
    intraHits = c(20, 20, 20, 20),
    interTargetFDRReached = TRUE,
    intraTargetFDRReached = TRUE,
    worstClassCorrelation = c(0.8, 0.8, 0.7, 0.75),
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
  expect_equal(
    selection$candidates$recoveryRelativeToBest,
    c(40 / 45, 1, 1, 50 / 55)
  )
})

test_that("radial recommendation remains separate when linear is ineligible", {
  candidates <- tibble::tibble(
    index = 1:3,
    kernel = c("linear", "radial", "radial"),
    cost = c(1, 1, 5),
    gamma = c(NA_real_, 0.1, 0.01),
    interInt = c(60, 40, 45),
    interHits = c(60, 40, 45),
    intraHits = c(20, 20, 20),
    interTargetFDRReached = TRUE,
    intraTargetFDRReached = TRUE,
    worstClassCorrelation = c(0.9, 0.7, 0.75),
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
    intraHits = c(20, 20, 20),
    interTargetFDRReached = TRUE,
    intraTargetFDRReached = TRUE,
    worstClassCorrelation = c(0.7, 0.8, 0.99),
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

test_that("ineligible candidates do not establish the recovery benchmark", {
  candidates <- tibble::tibble(
    index = 1:3,
    kernel = rep("linear", 3),
    cost = c(0.01, 0.1, 1),
    gamma = NA_real_,
    interHits = c(200, 100, 95),
    intraHits = c(50, 50, 50),
    interTargetFDRReached = TRUE,
    intraTargetFDRReached = TRUE,
    worstClassCorrelation = c(-0.2, 0.7, 0.8),
    eligible = c(FALSE, TRUE, TRUE)
  )

  selection <- touchstone:::selectSVMCandidates(candidates)

  expect_identical(selection$candidates$bestRecoveryHits, rep(100, 3))
  expect_identical(selection$recommended, 3L)
})

test_that("selection falls back to intraprotein recovery when inter hits are zero", {
  candidates <- tibble::tibble(
    index = 1:3,
    kernel = rep("linear", 3),
    cost = c(0.01, 0.1, 1),
    gamma = NA_real_,
    interHits = c(0, 0, 0),
    intraHits = c(100, 91, 89),
    interTargetFDRReached = FALSE,
    intraTargetFDRReached = TRUE,
    worstClassCorrelation = c(0.7, 0.8, 0.99),
    eligible = rep(TRUE, 3)
  )

  selection <- touchstone:::selectSVMCandidates(candidates)

  expect_identical(selection$candidates$selectionBasis, rep("intra", 3))
  expect_identical(selection$candidates$nearBestRecovery, c(TRUE, TRUE, FALSE))
  expect_identical(selection$recommended, 2L)
})

test_that("recovery ignores a class that did not reach its target FDR", {
  candidates <- tibble::tibble(
    index = 1:2,
    kernel = rep("linear", 2),
    cost = c(0.1, 1),
    gamma = NA_real_,
    interHits = c(10, 0),
    intraHits = c(90, 100),
    interTargetFDRReached = FALSE,
    intraTargetFDRReached = TRUE,
    worstClassCorrelation = c(0.8, 0.7),
    eligible = TRUE
  )

  selection <- touchstone:::selectSVMCandidates(candidates)

  expect_identical(selection$candidates$selectionBasis, rep("intra", 2))
  expect_identical(selection$candidates$bestRecoveryHits, rep(100, 2))
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
