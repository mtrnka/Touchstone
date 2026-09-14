test_that("pair construction is orientation-independent and counts support", {
  input <- data.frame(
    Acc.1 = c("A", "A", "B", "A"),
    Acc.2 = c("B", "B", "A", "B"),
    Protein.1 = c("protein A", "protein A", "protein B", "protein A"),
    Protein.2 = c("protein B", "protein B", "protein A", "protein B"),
    XLink.AA.1 = c(10, 10, 20, 11),
    XLink.AA.2 = c(20, 20, 10, 21),
    DB.Peptide.1 = c("AAA", "AAA", "BBB", "AAC"),
    DB.Peptide.2 = c("BBB", "BBB", "AAA", "BBC"),
    Score.Diff = c(20, 10, 18, 16)
  )

  result <- calculatePairs(input)

  expect_identical(as.character(result$xlinkedProtPair), rep("A::B", 4))
  expect_identical(
    as.character(result$xlinkedResPair),
    c("10.A::20.B", "10.A::20.B", "10.A::20.B", "11.A::21.B")
  )
  expect_equal(result$numCSM, c(3, 3, 3, 1))
  expect_equal(result$wtCSM, log1p(c(2, 2, 2, 1)))
  csm.weights <- stats::plogis((c(20, 10, 18) - 15) / 2.5)
  expect_equal(
    result$CSMsupport,
    c(
      log1p(sum(csm.weights[-1])),
      log1p(sum(csm.weights[-2])),
      log1p(sum(csm.weights[-3])),
      0
    )
  )
  expect_equal(result$numURP, rep(2, 4))
  expect_equal(result$wtURP, rep(log1p(2), 4))
  urp.weights <- stats::plogis((c(20, 16) - 15) / 2.5)
  expect_equal(
    result$URPsupport,
    c(rep(log1p(urp.weights[[2]]), 3), log1p(urp.weights[[1]]))
  )
})

test_that("support features exclude the evidence being scored", {
  input <- data.frame(
    Acc.1 = c("A", "A"),
    Acc.2 = c("B", "B"),
    Protein.1 = "protein A",
    Protein.2 = "protein B",
    XLink.AA.1 = c(10, 11),
    XLink.AA.2 = c(20, 21),
    DB.Peptide.1 = c("AAA", "AAC"),
    DB.Peptide.2 = c("BBB", "BBC"),
    Score.Diff = c(20, 16)
  )

  result <- calculatePairs(input)

  expect_equal(result$CSMsupport, c(0, 0))
  expect_equal(
    result$URPsupport,
    log1p(rev(stats::plogis((input$Score.Diff - 15) / 2.5)))
  )
})

test_that("Score.Diff values below five contribute no support", {
  input <- data.frame(
    Acc.1 = c("A", "A"),
    Acc.2 = c("B", "B"),
    Protein.1 = "protein A",
    Protein.2 = "protein B",
    XLink.AA.1 = 10,
    XLink.AA.2 = 20,
    DB.Peptide.1 = "AAA",
    DB.Peptide.2 = "BBB",
    Score.Diff = c(4.99, 15)
  )

  result <- calculatePairs(input)

  expect_equal(result$CSMsupport, c(log1p(0.5), 0))
  expect_equal(result$URPsupport, c(0, 0))
})

test_that("support construction does not use target or decoy labels", {
  input <- data.frame(
    Acc.1 = c("A", "A"),
    Acc.2 = c("B", "B"),
    Protein.1 = "protein A",
    Protein.2 = "protein B",
    XLink.AA.1 = 10,
    XLink.AA.2 = 20,
    DB.Peptide.1 = "AAA",
    DB.Peptide.2 = "BBB",
    Score.Diff = c(10, 20),
    Decoy = c("Target", "Target")
  )
  relabeled <- input
  relabeled$Decoy <- c("DoubleDecoy", "Decoy")

  target.result <- calculatePairs(input)
  decoy.result <- calculatePairs(relabeled)

  expect_equal(target.result$CSMsupport, decoy.result$CSMsupport)
  expect_equal(target.result$URPsupport, decoy.result$URPsupport)
})

test_that("residue-pair summarization retains the highest-scoring CSM", {
  input <- data.frame(
    xlinkedResPair = factor(c("A::B", "A::B", "C::D")),
    SVM.score = c(0.2, 0.9, 0.3),
    row_id = c("low", "high", "only")
  )

  result <- bestResPair(input)

  expect_identical(result$row_id, c("high", "only"))
  expect_equal(result$SVM.score, c(0.9, 0.3))
})
