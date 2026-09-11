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
    bestResPair = function(datTab) datTab,
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
