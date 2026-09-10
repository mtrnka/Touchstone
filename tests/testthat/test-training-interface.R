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
    make_model("radial", 0.1, 1, c(2, 1, 1, 1))
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
  expect_true(any(plot$data$hits == 1))
})

