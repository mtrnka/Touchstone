test_that("score correlation plot reproduces the CSM diagnostic", {
  csms <- tibble::tibble(
    SVM.score = c(-1, 0, 1, 2),
    Score.Diff = c(1, 2, 4, 5),
    Decoy = c("Decoy", "Target", "Target", "Target"),
    xlinkClass = c(
      "interProtein", "interProtein", "intraProtein", "intraProtein"
    )
  )

  plot <- plotScoreCorrelation(csms)

  expect_s3_class(plot, "ggplot")
  expect_identical(plot$data, csms)
  expect_identical(plot$labels$x, "SVM.score")
  expect_identical(plot$labels$y, "Score.Diff")
  expect_length(plot$layers, 1)
  expect_s3_class(plot$layers[[1]]$geom, "GeomPoint")
  expect_equal(length(unique(ggplot2::ggplot_build(plot)$layout$layout$PANEL)), 2)
})

test_that("score correlation plot selects requested training model", {
  linear <- list(
    CSMs = data.frame(
      SVM.score = 1,
      Score.Diff = 2,
      Decoy = "Target",
      xlinkClass = "interProtein"
    ),
    kernel = "linear",
    cost = 1,
    gamma = NA_real_
  )
  radial <- linear
  radial$CSMs$SVM.score <- 10
  radial$kernel <- "radial"
  radial$gamma <- 0.01
  training <- structure(
    list(
      recommended = linear,
      recommendedRadial = radial,
      candidates = tibble::tibble(
        index = 1:2,
        recommended = c(TRUE, FALSE),
        recommendedRadial = c(FALSE, TRUE)
      ),
      models = list(linear, radial),
      settings = list(scoreName = "SVM.score")
    ),
    class = "touchstone_training"
  )

  plot <- plotScoreCorrelation(training, model = "radial")

  expect_identical(plot$data$SVM.score, 10)
})

test_that("score correlation plot supports alternate classifier columns", {
  scored <- data.frame(
    experimental.score = c(1, 2),
    Score.Diff = c(3, 4),
    Decoy = c("Target", "Decoy"),
    xlinkClass = c("interProtein", "intraProtein")
  )

  plot <- plotScoreCorrelation(scored, classifier = "experimental.score")

  expect_identical(plot$labels$x, "experimental.score")
})
