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

  plot <- plotFDRHits(training)

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

  plot.without.raw <- plotFDRHits(training, showRaw = FALSE)
  expect_length(plot.without.raw$layers, 2)
})
