test_that("decoy labels distinguish target, mixed, and double-decoy matches", {
  input <- data.frame(
    Acc.1 = c("P1", "r1_P1", "decoy_P1", "P1"),
    Acc.2 = c("P2", "P2", "DECOY_P2", "DECOY_P2")
  )

  result <- calculateDecoys(input)

  expect_identical(
    as.character(result$Decoy),
    c("Target", "Decoy", "DoubleDecoy", "Decoy")
  )
  expect_identical(
    as.character(result$Decoy2),
    c("Target", "Decoy", "Decoy", "Decoy")
  )
  expect_identical(
    levels(result$Decoy),
    c("DoubleDecoy", "Decoy", "Target")
  )
})

test_that("decoy fractions and unseparated FDR match a hand calculation", {
  scored <- data.frame(
    Decoy = factor(
      c(rep("DoubleDecoy", 4), rep("Decoy", 10), rep("Target", 100)),
      levels = c("DoubleDecoy", "Decoy", "Target")
    ),
    SVM.score = 1
  )

  fractions <- touchstone:::calculateDecoyFractions(
    scored,
    scalingFactor = 2
  )

  # Four DD / 2^2 = 1; ten mixed / 2 - 2 * one DD-equivalent = 3.
  expect_equal(fractions, c(TT = 100, ftTT = 3, ffTT = 1))
  expect_equal(
    calculateFDR.unseparated(scored, threshold = 0, scalingFactor = 2),
    0.04
  )
})

test_that("negative mixed-decoy estimates use the documented correction", {
  scored <- data.frame(
    Decoy = factor(
      c(rep("DoubleDecoy", 4), rep("Decoy", 2), rep("Target", 20)),
      levels = c("DoubleDecoy", "Decoy", "Target")
    )
  )

  fractions <- touchstone:::calculateDecoyFractions(
    scored,
    scalingFactor = 2
  )

  expect_equal(fractions, c(TT = 20, ftTT = 0, ffTT = 0.5))
})

test_that("countDecoys uses an explicitly supplied scaling factor", {
  scored <- tibble::tibble(
    xlinkClass = c("interProtein", "interProtein"),
    Decoy = c("Target", "Decoy")
  )
  observed <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    deScaler = function(datTab, scalingFactor, ...) {
      observed$scalingFactor <- scalingFactor
      datTab
    },
    .package = "touchstone"
  )

  countDecoys(scored, scalingFactor = 5)

  expect_identical(observed$scalingFactor, 5)
})

test_that("deScaler is reproducible without changing the caller RNG state", {
  scored <- data.frame(
    row = seq_len(30),
    xlinkClass = rep(c("interProtein", "intraProtein"), 15),
    Decoy = factor(
      c(rep("Target", 10), rep("Decoy", 10), rep("DoubleDecoy", 10)),
      levels = c("DoubleDecoy", "Decoy", "Target")
    )
  )

  first <- deScaler(scored, scalingFactor = 2)
  stats::runif(5)
  second <- deScaler(scored, scalingFactor = 2)
  expect_identical(first, second)

  set.seed(917)
  expected.next <- stats::runif(1)
  set.seed(917)
  deScaler(scored, scalingFactor = 2)
  observed.next <- stats::runif(1)
  expect_identical(observed.next, expected.next)
})

test_that("deScaler supports an explicit seed", {
  scored <- data.frame(
    row = seq_len(100),
    Decoy = factor(
      rep(c("Decoy", "DoubleDecoy"), each = 50),
      levels = c("DoubleDecoy", "Decoy", "Target")
    )
  )

  first <- deScaler(scored, scalingFactor = 2, seed = 11)
  second <- deScaler(scored, scalingFactor = 2, seed = 12)

  expect_false(identical(first$row, second$row))
})

test_that("countDecoys accepts prepared results and their classifier", {
  scored <- data.frame(
    xlinkClass = c("interProtein", "interProtein", "intraProtein"),
    Decoy = factor(
      c("Target", "Decoy", "Target"),
      levels = c("DoubleDecoy", "Decoy", "Target")
    ),
    SVM.score = c(10, 10, 10),
    experimental.score = c(2, -2, 2)
  )
  prepared <- structure(
    list(
      data = scored,
      thresholds = list(intraThresh = 0, interThresh = 0),
      settings = list(scalingFactor = 1, classifier = "experimental.score")
    ),
    class = "touchstone_results"
  )

  expected <- countDecoys(
    scored,
    threshold = prepared$thresholds,
    classifier = "experimental.score",
    scalingFactor = 1
  )

  expect_identical(countDecoys(prepared), expected)
})

test_that("calculateFDR accepts prepared Touchstone results", {
  scored <- data.frame(
    Decoy = factor(
      c(rep("DoubleDecoy", 4), rep("Decoy", 10), rep("Target", 100)),
      levels = c("DoubleDecoy", "Decoy", "Target")
    ),
    xlinkClass = "interProtein",
    SVM.score = 1
  )
  prepared <- structure(
    list(
      data = scored,
      thresholds = list(intraThresh = 0, interThresh = 0),
      settings = list(scalingFactor = 2)
    ),
    class = "touchstone_results"
  )

  expect_equal(calculateFDR(prepared), 0.04)
  expect_equal(
    calculateFDR(
      prepared,
      threshold = list(intraThresh = 0, interThresh = 0),
      scalingFactor = 2
    ),
    0.04
  )
})

test_that("calculateFDR uses the classifier stored in prepared results", {
  scored <- data.frame(
    Decoy = factor(
      c(rep("DoubleDecoy", 4), rep("Decoy", 10), rep("Target", 100)),
      levels = c("DoubleDecoy", "Decoy", "Target")
    ),
    xlinkClass = "interProtein",
    experimental.score = 1
  )
  prepared <- structure(
    list(
      data = scored,
      thresholds = list(intraThresh = 0, interThresh = 0),
      settings = list(scalingFactor = 2, classifier = "experimental.score")
    ),
    class = "touchstone_results"
  )

  expect_equal(calculateFDR(prepared), 0.04)
})

test_that("classification accepts default, bare, and character score columns", {
  scored <- data.frame(
    xlinkClass = c("interProtein", "intraProtein"),
    SVM.score = c(2, -2),
    Score.Diff = c(-2, 2)
  )
  thresholds <- list(interThresh = 0, intraThresh = 0)
  score.column <- "Score.Diff"

  default.result <- classifyDataset(scored, thresholds)
  bare.result <- classifyDataset(scored, thresholds, classifier = Score.Diff)
  character.result <- classifyDataset(
    scored, thresholds, classifier = score.column
  )

  expect_identical(default.result$SVM.score, 2)
  expect_identical(bare.result$Score.Diff, 2)
  expect_identical(character.result$Score.Diff, 2)
})

test_that("separate modelled thresholds tolerate an absent crosslink class", {
  intra.only <- data.frame(
    xlinkClass = rep("intraProtein", 3),
    SVM.score = 1:3
  )

  testthat::local_mocked_bindings(
    findThreshold = function(...) list(globalThresh = 1),
    findThresholdModelled = function(...) {
      stop("The absent interprotein class should not be modelled")
    },
    .package = "touchstone"
  )

  thresholds <- findSeparateThresholdsModelled(intra.only, plot = FALSE)

  expect_identical(thresholds$intraThresh, 1)
  expect_identical(thresholds$interThresh, Inf)
})

test_that("logistic FDR modelling is constrained to a decreasing curve", {
  decoy.table <- tibble::tibble(
    thresh = seq(-4, 4, length.out = 100),
    fdr.exp = 0.2 / (1 + exp(1.25 * (thresh - 0.5)))
  )
  decoy.table$fdr.orig <- decoy.table$fdr.exp

  result <- touchstone:::modelFDRCurve(
    decoy.table,
    maxFDR = 0.2,
    targetER = 0.01
  )

  expect_identical(result$method, "logistic-model")
  expect_true(all(is.finite(result$fdr)))
  expect_true(all(diff(result$fdr) <= 0))
})

test_that("flat FDR data use a conservative monotonic empirical fallback", {
  decoy.table <- tibble::tibble(
    thresh = 1:5,
    fdr.exp = rep(0.1, 5),
    fdr.orig = c(0.2, 0.1, 0.15, -0.1, 0.02)
  )

  result <- touchstone:::modelFDRCurve(
    decoy.table,
    maxFDR = 0.2,
    targetER = 0.01
  )

  expect_identical(result$method, "monotonic-empirical-fallback")
  expect_equal(result$fdr, c(0.2, 0.15, 0.15, 0.02, 0.02))
  expect_true(all(result$fdr >= 0))
  expect_true(all(diff(result$fdr) <= 0))
})

test_that("modelled threshold records conservative fallback diagnostics", {
  decoy.table <- tibble::tibble(
    thresh = 1:4,
    fdr = c(0.2, 0.05, 0.01, 0.005)
  )
  attr(decoy.table, "thresholdMethod") <-
    "monotonic-empirical-fallback"
  attr(decoy.table, "thresholdMessage") <- "test fit failure"

  testthat::local_mocked_bindings(
    generateDecoyTable = function(...) decoy.table,
    .package = "touchstone"
  )

  expect_warning(
    result <- findThresholdModelled(
      data.frame(SVM.score = 1:4), targetER = 0.01
    ),
    "conservative monotonic empirical"
  )

  expect_identical(result$globalThresh, 3L)
  expect_identical(result$method, "monotonic-empirical-fallback")
  expect_true(result$targetFDRReached)
})

test_that("separate modelled thresholds retain per-class methods", {
  scored <- data.frame(
    xlinkClass = c("interProtein", "intraProtein"),
    SVM.score = c(1, 1)
  )

  testthat::local_mocked_bindings(
    findThresholdModelled = function(...) {
      list(
        globalThresh = 2,
        correspondingFDR = 0.01,
        method = "monotonic-empirical-fallback",
        targetFDRReached = TRUE
      )
    },
    findThreshold = function(...) {
      list(globalThresh = 1, correspondingFDR = 0.005)
    },
    .package = "touchstone"
  )

  thresholds <- findSeparateThresholdsModelled(scored, plot = FALSE)

  expect_equal(thresholds, list(intraThresh = 1, interThresh = 2),
               ignore_attr = TRUE)
  expect_identical(
    attr(thresholds, "thresholdMethods"),
    c(
      intraProtein = "empirical",
      interProtein = "monotonic-empirical-fallback"
    )
  )
  expect_identical(
    attr(thresholds, "targetFDRReached"),
    c(intraProtein = TRUE, interProtein = TRUE)
  )
})
