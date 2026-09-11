test_that("fdrPlots uses settings stored in prepared results", {
  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  scored <- tibble::tibble(
    experimental.score = c(-2, -1, 1, 2),
    Decoy = factor(
      c("Target", "Decoy", "Target", "DoubleDecoy"),
      levels = c("DoubleDecoy", "Decoy", "Target")
    ),
    xlinkClass = c(
      "interProtein", "interProtein", "intraProtein", "intraProtein"
    )
  )
  prepared <- structure(
    list(
      data = scored,
      thresholds = list(intraThresh = 0.5, interThresh = -0.5),
      fdr = 0.0123,
      summarizationLevel = "urp",
      settings = list(
        classifier = "experimental.score",
        scalingFactor = 5
      )
    ),
    class = "touchstone_results"
  )
  observed <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    deScaler = function(datTab, scalingFactor, ...) {
      observed$scalingFactor <- scalingFactor
      datTab
    },
    calculateFDR = function(datTab, threshold, classifier, scalingFactor) {
      observed$fdr.threshold <- threshold
      observed$fdr.classifier <- classifier
      observed$fdr.scalingFactor <- scalingFactor
      0.0456
    },
    .package = "touchstone"
  )

  plot <- fdrPlots(prepared)
  built <- suppressWarnings(ggplot2::ggplot_build(plot))
  reference.lines <- unlist(lapply(built$data[-1], function(layer) {
    if ("xintercept" %in% names(layer)) unique(layer$xintercept) else NULL
  }))

  expect_s3_class(plot, "ggplot")
  expect_identical(plot$data$experimental.score, scored$experimental.score)
  expect_identical(observed$scalingFactor, 5)
  expect_identical(observed$fdr.threshold, prepared$thresholds)
  expect_identical(observed$fdr.classifier, "experimental.score")
  expect_identical(observed$fdr.scalingFactor, 5)
  expect_setequal(reference.lines, c(0.5, -0.5))
  expect_identical(
    plot$labels$subtitle,
    "Summarization level: URP | Calculated FDR: 4.56%"
  )
})

test_that("fdrPlots recalculates subtitle FDR at an override threshold", {
  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  scored <- tibble::tibble(
    SVM.score = c(-2, -1, 1, 2),
    Decoy = c("Target", "Decoy", "Target", "DoubleDecoy"),
    xlinkClass = c(
      "interProtein", "interProtein", "intraProtein", "intraProtein"
    )
  )
  prepared <- structure(
    list(
      data = scored,
      thresholds = list(intraThresh = 0.5, interThresh = -0.5),
      fdr = 0.99,
      summarizationLevel = "urp",
      settings = list(classifier = "SVM.score", scalingFactor = 1)
    ),
    class = "touchstone_results"
  )
  override <- list(intraThresh = 1.5, interThresh = -1.5)
  observed <- new.env(parent = emptyenv())

  testthat::local_mocked_bindings(
    calculateFDR = function(datTab, threshold, classifier, scalingFactor) {
      observed$threshold <- threshold
      0.0789
    },
    .package = "touchstone"
  )

  plot <- fdrPlots(prepared, threshold = override)
  built <- suppressWarnings(ggplot2::ggplot_build(plot))
  reference.lines <- unlist(lapply(built$data[-1], function(layer) {
    if ("xintercept" %in% names(layer)) unique(layer$xintercept) else NULL
  }))

  expect_identical(observed$threshold, override)
  expect_setequal(reference.lines, c(1.5, -1.5))
  expect_identical(
    plot$labels$subtitle,
    "Summarization level: URP | Calculated FDR: 7.89%"
  )
})

test_that("fdrPlots retains its data-frame interface", {
  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  scored <- tibble::tibble(
    Score.Diff = c(-2, -1, 1, 2),
    Decoy = c("Target", "Decoy", "Target", "DoubleDecoy"),
    xlinkClass = c(
      "interProtein", "interProtein", "intraProtein", "intraProtein"
    )
  )

  plot <- fdrPlots(
    scored,
    threshold = 0,
    classifier = "Score.Diff",
    scalingFactor = 1
  )

  expect_s3_class(plot, "ggplot")
})

test_that("fdrPlots orders all match classes from largest to smallest", {
  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  scored <- tibble::tibble(
    SVM.score = seq_len(10),
    Decoy = c(rep("Decoy", 5), rep("Target", 3), rep("DoubleDecoy", 2)),
    xlinkClass = "interProtein"
  )

  plot <- fdrPlots(scored, scalingFactor = 1)

  expect_identical(
    levels(plot$data$Decoy),
    c("Decoy", "Target", "DoubleDecoy")
  )
})

test_that("fdrPlots supports x-axis zoom and dodged histograms", {
  grDevices::pdf(tempfile(fileext = ".pdf"))
  on.exit(grDevices::dev.off(), add = TRUE)
  scored <- tibble::tibble(
    SVM.score = c(-5, -2, 1, 5),
    Decoy = c("Target", "Decoy", "Target", "DoubleDecoy"),
    xlinkClass = c(
      "interProtein", "interProtein", "intraProtein", "intraProtein"
    )
  )

  plot <- fdrPlots(
    scored,
    scalingFactor = 1,
    xLimits = c(-2, 2),
    histogramPosition = "dodge"
  )

  expect_s3_class(plot$layers[[1]]$position, "PositionDodge2")
  expect_equal(plot$coordinates$limits$x, c(-2, 2))
  expect_setequal(plot$data$SVM.score, c(-2, 1))
  expect_error(fdrPlots(scored, xLimits = c(2, -2)), "increasing")
  expect_error(
    fdrPlots(scored, xLimits = c(10, 20)),
    "No observations"
  )
})
