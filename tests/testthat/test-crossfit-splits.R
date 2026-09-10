test_that("cross-fitting keeps residue pairs together and apart from training", {
  input <- data.frame(
    xlinkedResPair = rep(LETTERS[1:8], each = 2),
    Decoy2 = factor(rep(c("Target", "Target", "Decoy", "Decoy"), each = 4))
  )

  split <- touchstone:::makeCrossfitSplit(
    input,
    sampleNo = 6,
    seed = 42
  )

  train.1.groups <- unique(split$group.id[split$train.1])
  train.2.groups <- unique(split$group.id[split$train.2])

  expect_identical(split$splitBy, "xlinkedResPair")
  expect_length(intersect(train.1.groups, train.2.groups), 0)
  expect_false(any(split$score.1[split$group.id %in% train.1.groups]))
  expect_false(any(split$score.2[split$group.id %in% train.2.groups]))
  expect_true(all(split$score.1[split$group.id %in% train.2.groups]))
  expect_true(all(split$score.2[split$group.id %in% train.1.groups]))
  expect_lte(length(split$train.1), 6)
  expect_lte(length(split$train.2), 6)
  expect_setequal(unique(input$Decoy2[split$train.1]), c("Decoy", "Target"))
  expect_setequal(unique(input$Decoy2[split$train.2]), c("Decoy", "Target"))
})

test_that("cross-fitting is reproducible without changing the caller RNG", {
  input <- data.frame(xlinkedResPair = rep(LETTERS[1:10], each = 2))

  set.seed(99)
  expected.next.draw <- runif(1)
  set.seed(99)

  split.1 <- touchstone:::makeCrossfitSplit(input, sampleNo = 6, seed = 7)
  observed.next.draw <- runif(1)
  split.2 <- touchstone:::makeCrossfitSplit(input, sampleNo = 6, seed = 7)

  expect_identical(split.1, split.2)
  expect_equal(observed.next.draw, expected.next.draw)
})

test_that("spectrum identity is used when residue pairs are unavailable", {
  input <- data.frame(
    Fraction = rep(c("fraction1", "fraction2"), each = 4),
    Spectrum = rep(rep(1:2, each = 2), 2)
  )

  split <- touchstone:::makeCrossfitSplit(input, sampleNo = 4, seed = 1)

  expect_identical(split$splitBy, c("Fraction", "Spectrum"))
  expect_equal(length(unique(split$group.id)), 4)
})

test_that("multi-column group identifiers do not collide on punctuation", {
  input <- data.frame(
    first = c("a.b", "a"),
    second = c("c", "b.c")
  )

  split <- touchstone:::makeCrossfitSplit(
    input,
    sampleNo = 1,
    splitBy = c("first", "second"),
    seed = 1
  )

  expect_equal(length(unique(split$group.id)), 2)
})

test_that("cross-fitting rejects outcome classes with only one independent group", {
  input <- data.frame(
    xlinkedResPair = c("target1", "target2", "decoy1"),
    Decoy2 = factor(c("Target", "Target", "Decoy"))
  )

  expect_error(
    touchstone:::makeCrossfitSplit(input, seed = 1),
    "at least two independent splitBy groups"
  )
})

test_that("buildSVM produces reproducible cross-fitted scores", {
  input <- data.frame(
    xlinkedResPair = rep(sprintf("pair%02d", 1:20), each = 2),
    Decoy2 = factor(
      rep(rep(c("Target", "Decoy"), each = 10), each = 2),
      levels = c("Decoy", "Target")
    ),
    Score.Diff = c(seq(12, 21.5, length.out = 20),
                   seq(1, 10.5, length.out = 20)),
    percMatched = seq(0.1, 0.9, length.out = 40),
    ppm = rep(seq(-2, 2, length.out = 20), 2)
  )

  expect_silent(
    score.1 <- buildSVM(
      input,
      params = c("Score.Diff", "percMatched", "massError"),
      sampleNo = 16,
      seed = 11,
      kernel = "linear",
      cost = 1
    )
  )
  expect_silent(
    score.2 <- buildSVM(
      input,
      params = c("Score.Diff", "percMatched", "massError"),
      sampleNo = 16,
      seed = 11,
      kernel = "linear",
      cost = 1
    )
  )

  expect_identical(score.1$SVM.score, score.2$SVM.score)
  expect_true(all(is.finite(score.1$SVM.score)))
})
