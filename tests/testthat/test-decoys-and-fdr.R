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

