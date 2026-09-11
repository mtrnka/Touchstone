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
