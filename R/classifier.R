#' Train and select a crosslink-scoring SVM
#'
#' Performs feature selection, optional Score.Diff prefilter selection, and a
#' reproducible grid search. Linear SVMs are the conservative default. Candidate
#' models must reach the requested interprotein FDR with nonzero recovery and
#' remain positively correlated with Score.Diff. Eligible linear candidates are
#' ranked by recovery at the target FDR, then by average recovery over the
#' low-FDR region and lower cost.
#'
#' @param datTab Parsed CLMS search results
#' @param params Character vector specifying names of the features in `datTab` used to train model.
#' @param scoreName Name for the new scoring function.
#' @param scalingFactor An integer k. The multiple by which decoy DB is larger than target DB
#' @param targetER Desired FDR for classification of CSMs
#' @param sampleNo Size of the training dataset (integer).
#' @param cost_values Numeric vector of cost values used for hyperparameter tuning of the SVM model
#' @param gamma_values Numeric vector of gamma values used only when radial
#'   kernels are explicitly requested.
#' @param sd_values Numeric vector of Score Diff values to use for prefilitering optimiziation.
#' @param kernels Character vector of SVM kernels to evaluate. The conservative
#'   default is `"linear"`; include `"radial"` to evaluate experimental radial
#'   candidates without making them eligible for automatic recommendation.
#' @param seed Integer seed used to make cross-fitting reproducible.
#' @param splitBy Character vector naming columns whose rows must remain together
#'   during cross-fitting. The default uses residue pairs when available, then a
#'   spectrum identifier, and finally individual rows.
#' @param verbose Print progress and the candidate table.
#' @seealso [tuneSVM()], [tuneSVM.helper()], [buildSVM()]
#' @returns A `touchstone_training` object containing the recommended candidate,
#'   the candidate audit table, all fitted candidates, prefilter information,
#'   and training settings.
#' @export
trainCrosslinkScore <- function(datTab,
                                params = NULL,
                                scoreName="SVM.score",
                                scalingFactor = the$decoyScalingFactor,
                                targetER = 0.01,
                                sampleNo = 20000,
                                cost_values = c(1, 5, 10),
                                gamma_values = c(0.01, 0.05, 0.1),
                                sd_values = c(0,5,10,15,20),
                                kernels = "linear",
                                seed = 1,
                                splitBy = NULL,
                                verbose = FALSE) {
  datTab <- dplyr::ungroup(datTab)

  # feature selection
  plausibleHits <- datTab %>%
    filter(.data$Decoy == "Target",
           .data$Score.Diff > 10,
           .data$numCSM > 1,
           .data$xlinkClass == "intraProtein") %>%
    dplyr::count(.data$Acc.1, name = "n") %>%
    dplyr::arrange(desc(.data$n))

  if (is.null(params)) {
    params <- dplyr::case_when(
      nrow(plausibleHits) <= 50 ~ list(params.best.nop),
      nrow(plausibleHits) >= 2 &&
        plausibleHits$n[1] > 100 * plausibleHits$n[2] ~ list(params.best.nop),
      !("Perc.Bond.Cleavage.1" %in% names(datTab)) ~ list(params.noPercBond),
      TRUE ~ list(params.best)
    ) %>%
      unlist()
  }

  # prefiltering
  nProteinPairs <- removeDecoys(datTab) %>%
    dplyr::distinct(.data$Acc.1, .data$Acc.2) %>%
    nrow()

  preFilter.summary <- NULL
  bestPreFilter <- min(sd_values, na.rm = TRUE)

  if (nProteinPairs > 1000) {
    prefilter <- chooseScoreDiffPrefilter(
      datTab = datTab,
      sd_values = sd_values,
      targetER = targetER,
      params = params,
      scoreName = scoreName,
      scalingFactor = scalingFactor,
      sampleNo = sampleNo,
      cost = min(cost_values),
      gamma = NA_real_,
      kernel = "linear",
      seed = seed,
      splitBy = splitBy,
      verbose = verbose,
      fallback.threshold = min(sd_values, na.rm = TRUE)
    )

    datTab <- prefilter$datTab
    preFilter.summary <- prefilter$preFilter.summary
    bestPreFilter <- prefilter$bestPreFilter
  }

  # Hyperparameter optimization
  if (verbose) {
    message("Hyperparameter optimization...")
  }
  tuned <- tuneSVM(datTab,
                   params=params,
                   scoreName=scoreName,
                   scalingFactor = scalingFactor,
                   targetER = targetER,
                   sampleNo = sampleNo,
                   cost_values = cost_values,
                   gamma_values = gamma_values,
                   kernels = kernels,
                   seed = seed,
                   splitBy = splitBy,
                   verbose = verbose)
  tuned.parse <- tuned %>%
    purrr::imap_dfr(function(x,i) {
      data.frame(
        "index" = i,
        "kernel" = x$kernel,
        "cost" = x$cost,
        "gamma" = x$gamma,
        "interInt" = x$interInt,
        "interHits" = x$interHits,
        "scoreCorrelation" = x$corScore / 100,
        "targetFDRReached" = any(
          x$errorTable$fdr.inter <= targetER & x$errorTable$inter > 0,
          na.rm = TRUE
        )
      )
    }) %>%
    mutate(
      eligible = .data$kernel == "linear" &
        .data$targetFDRReached &
        is.finite(.data$interInt) &
        .data$interHits > 0 &
        is.finite(.data$scoreCorrelation) &
        .data$scoreCorrelation > 0,
      rejectionReason = dplyr::case_when(
        .data$kernel != "linear" ~
          "Radial candidates are inspection-only",
        !.data$targetFDRReached ~
          "Target interprotein FDR was not reached with nonzero hits",
        !is.finite(.data$interInt) ~
          "FDR-versus-hit summary was not finite",
        .data$interHits <= 0 ~
          "No interprotein hits at the target FDR",
        !is.finite(.data$scoreCorrelation) ~
          "Score correlation was not finite",
        .data$scoreCorrelation <= 0 ~
          "SVM score was not positively correlated with Score.Diff",
        TRUE ~ NA_character_
      )
    )

  eligible.models <- tuned.parse %>%
    filter(.data$eligible) %>%
    arrange(
      desc(.data$interHits),
      desc(.data$interInt),
      .data$cost,
      desc(.data$scoreCorrelation)
    )

  recommended.index <- if (nrow(eligible.models) > 0) {
    eligible.models$index[[1]]
  } else {
    NA_integer_
  }

  tuned.parse <- tuned.parse %>%
    mutate(
      recommended = if (is.na(recommended.index)) {
        FALSE
      } else {
        .data$index == recommended.index
      }
    ) %>%
    arrange(desc(.data$recommended), desc(.data$eligible), .data$kernel,
            .data$cost, .data$gamma)

  if (verbose) {
    print(tuned.parse)
  }

  structure(
    list(
      recommended = if (is.na(recommended.index)) NULL else tuned[[recommended.index]],
      recommendedIndex = recommended.index,
      candidates = tuned.parse,
      models = tuned,
      prefilter = list(
        selectedScoreDiff = bestPreFilter,
        candidates = preFilter.summary
      ),
      settings = list(
        targetER = targetER,
        kernels = kernels,
        seed = seed,
        splitBy = splitBy,
        features = params
      )
    ),
    class = "touchstone_training"
  )
}

#' Print a Touchstone training result
#'
#' @param x A result returned by [trainCrosslinkScore()].
#' @param ... Additional arguments passed to `print()` for the candidate table.
#' @return `x`, invisibly.
#' @export
print.touchstone_training <- function(x, ...) {
  if (is.null(x$recommended)) {
    cat("Touchstone training result: no eligible model was recommended.\n")
  } else {
    cat("Touchstone training result: recommended candidate ",
        x$recommendedIndex, ".\n", sep = "")
  }
  print(x$candidates, ...)
  invisible(x)
}

#' Plot FDR versus the number of recovered crosslinks
#'
#' Creates the faceted diagnostic plot formerly printed automatically by
#' `trainCrosslinkScore()`. Linear and radial/gamma model families occupy
#' separate facets, and color denotes SVM cost. Weak candidates are retained so
#' failed or unstable model families remain visible during inspection.
#'
#' @param training Result returned by `trainCrosslinkScore()`, or the model list
#'   returned by `tuneSVM()`.
#' @param targetER Desired FDR shown by the vertical reference line. By default,
#'   uses the value stored in a `trainCrosslinkScore()` result.
#' @param maxFDR Largest FDR value displayed.
#' @param linkClass Plot `"inter"` or `"intra"` protein crosslinks.
#' @return A `ggplot2` plot.
#' @export
plotFDRHits <- function(training,
                        targetER = NULL,
                        maxFDR = 0.05,
                        linkClass = c("inter", "intra")) {
  linkClass <- match.arg(linkClass)

  if (inherits(training, "touchstone_training")) {
    models <- training$models
    if (is.null(targetER)) {
      targetER <- training$settings$targetER
    }
  } else if (is.list(training)) {
    models <- training
  } else {
    stop(
      "training must be a trainCrosslinkScore() result or tuneSVM() model list.",
      call. = FALSE
    )
  }

  if (is.null(targetER)) {
    targetER <- 0.01
  }
  if (length(maxFDR) != 1 || !is.finite(maxFDR) || maxFDR <= 0) {
    stop("maxFDR must be one positive, finite number.", call. = FALSE)
  }

  fdr.column <- paste0("fdr.", linkClass)
  hit.column <- linkClass

  plot.data <- purrr::imap_dfr(models, function(model, index) {
    required <- c("errorTable", "kernel", "cost", "gamma")
    if (!all(required %in% names(model)) ||
        !all(c(fdr.column, hit.column) %in% names(model$errorTable))) {
      stop(
        "Every candidate must contain model settings and an FDR error table.",
        call. = FALSE
      )
    }

    model.label <- if (identical(model$kernel, "linear")) {
      "linear"
    } else {
      paste0("radial (gamma = ", format(model$gamma), ")")
    }
    model.cost <- model$cost

    model$errorTable %>%
      dplyr::transmute(
        candidate = index,
        model = model.label,
        cost = factor(model.cost),
        fdr = .data[[fdr.column]],
        hits = .data[[hit.column]]
      )
  }) %>%
    filter(is.finite(.data$fdr), is.finite(.data$hits),
           .data$fdr >= 0, .data$fdr <= maxFDR)

  if (nrow(plot.data) == 0) {
    stop("No finite FDR-versus-hit values fall within maxFDR.", call. = FALSE)
  }

  plot.data$model <- factor(plot.data$model, levels = unique(plot.data$model))

  ggplot2::ggplot(
    plot.data,
    ggplot2::aes(x = .data$fdr, y = .data$hits,
                 color = .data$cost, group = .data$candidate)
  ) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::geom_vline(xintercept = targetER, color = "red") +
    ggplot2::scale_color_viridis_d(option = "C", name = "Cost") +
    ggplot2::facet_grid(rows = ggplot2::vars(.data$model), scales = "free_y") +
    ggplot2::labs(
      x = paste0(toupper(linkClass), "-protein FDR"),
      y = paste0(toupper(linkClass), "-protein crosslinks")
    ) +
    ggplot2::theme_bw()
}

as_scalar_numeric <- function(x, default = NA_real_) {
  if (is.null(x) || length(x) == 0) {
    return(default)
  }

  x <- suppressWarnings(as.numeric(x[1]))

  if (is.na(x) || is.nan(x) || !is.finite(x)) {
    return(default)
  }

  x
}

chooseScoreDiffPrefilter <- function(datTab,
                                     sd_values = c(0, 5, 10, 15, 20),
                                     targetER = 0.01,
                                     params,
                                     scoreName = "SVM.score",
                                     scalingFactor = the$decoyScalingFactor,
                                     sampleNo = 20000,
                                     cost = 1,
                                     gamma = NA_real_,
                                     kernel = "linear",
                                     seed = 1,
                                     splitBy = NULL,
                                     verbose = FALSE,
                                     class.col = NULL,
                                     target.label = "Target",
                                     min.total = 500,
                                     min.target = 50,
                                     min.decoy = 50,
                                     fallback.threshold = NULL) {

  if (!"Score.Diff" %in% names(datTab)) {
    stop("datTab must contain a Score.Diff column.", call. = FALSE)
  }

  if (is.null(class.col)) {
    class.col <- dplyr::case_when(
      "Decoy2" %in% names(datTab) ~ "Decoy2",
      TRUE ~ NA_character_
    )
  }

  if (is.na(class.col) || !class.col %in% names(datTab)) {
    stop("Could not find a class column. Expected Decoy2", call. = FALSE)
  }

  if (is.null(fallback.threshold)) {
    fallback.threshold <- min(sd_values, na.rm = TRUE)
  }

  sd_values <- sort(unique(sd_values))
  sd_values <- sd_values[!is.na(sd_values)]

  if (length(sd_values) == 0) {
    warning("No valid sd_values supplied; using fallback threshold.", call. = FALSE)
    sd_values <- fallback.threshold
  }

  if (verbose) {
    message("Score.Diff prefilter optimization...")
  }

  prefilter.results <- purrr::map(sd_values, function(sd) {
    datTab.pre <- datTab %>%
      dplyr::filter(.data$Score.Diff >= sd)
    class.values <- datTab.pre[[class.col]]
    n.total <- nrow(datTab.pre)
    n.target <- sum(class.values == target.label, na.rm = TRUE)
    n.decoy <- sum(class.values != target.label & !is.na(class.values))
    n.inter.target <- if ("xlinkClass" %in% names(datTab.pre)) {
      sum(
        datTab.pre[[class.col]] == target.label &
          datTab.pre$xlinkClass == "interProtein",
        na.rm = TRUE
      )
    } else {
      NA_integer_
    }

    valid.training <- n.total >= min.total &&
      n.target >= min.target &&
      n.decoy >= min.decoy

    if (!valid.training) {
      return(list(
        result = NULL,
        summary = tibble::tibble(
          sd.thresh = sd,
          n.total = n.total,
          n.target = n.target,
          n.decoy = n.decoy,
          n.inter.target = n.inter.target,
          valid.training = FALSE,
          valid.result = FALSE,
          interInt = NA_real_,
          interHits = NA_real_,
          error = "Insufficient training data after Score.Diff prefilter"
        )
      ))
    }

    fit <- tryCatch(
      {
        tuneSVM.helper(
          datTab = datTab.pre,
          targetER = targetER,
          params = params,
          scoreName = scoreName,
          scalingFactor = scalingFactor,
          sampleNo = sampleNo,
          cost = cost,
          gamma = gamma,
          kernel = kernel,
          seed = seed,
          splitBy = splitBy,
          verbose = verbose
        )
      },
      error = function(e) e
    )

    if (inherits(fit, "error")) {
      return(list(
        result = NULL,
        summary = tibble::tibble(
          sd.thresh = sd,
          n.total = n.total,
          n.target = n.target,
          n.decoy = n.decoy,
          n.inter.target = n.inter.target,
          valid.training = TRUE,
          valid.result = FALSE,
          interInt = NA_real_,
          interHits = NA_real_,
          error = conditionMessage(fit)
        )
      ))
    }

    interInt <- as_scalar_numeric(fit$interInt)
    interHits <- as_scalar_numeric(fit$interHits)

    valid.result <- is.finite(interInt) && !is.na(interInt)

    list(
      result = fit,
      summary = tibble::tibble(
        sd.thresh = sd,
        n.total = n.total,
        n.target = n.target,
        n.decoy = n.decoy,
        n.inter.target = n.inter.target,
        valid.training = TRUE,
        valid.result = valid.result,
        interInt = interInt,
        interHits = interHits,
        error = NA_character_
      )
    )
  })

  preFilter.summary <- purrr::map_dfr(prefilter.results, "summary") %>%
    dplyr::arrange(.data$sd.thresh)

  valid.results <- preFilter.summary %>%
    dplyr::filter(.data$valid.training, .data$valid.result)

  positive.inter <- valid.results %>%
    dplyr::filter(.data$interInt > 0)

  if (nrow(positive.inter) > 0) {
    max.inter <- max(positive.inter$interInt, na.rm = TRUE)

    bestPreFilter <- positive.inter %>%
      dplyr::filter(.data$interInt >= 0.95 * max.inter) %>%
      dplyr::summarise(best = min(.data$sd.thresh, na.rm = TRUE)) %>%
      dplyr::pull(.data$best)

  } else if (nrow(valid.results) > 0) {
    warning(
      "No Score.Diff prefilter produced positive interprotein recovery. ",
      "Using the mildest valid Score.Diff prefilter.",
      call. = FALSE
    )

    bestPreFilter <- min(valid.results$sd.thresh, na.rm = TRUE)

  } else {
    warning(
      "No Score.Diff prefilter produced a valid model. ",
      "Using fallback Score.Diff threshold: ",
      fallback.threshold,
      call. = FALSE
    )

    bestPreFilter <- fallback.threshold
  }

  if (!is.finite(bestPreFilter)) {
    warning(
      "Selected Score.Diff prefilter was not finite. ",
      "Using fallback Score.Diff threshold: ",
      fallback.threshold,
      call. = FALSE
    )

    bestPreFilter <- fallback.threshold
  }

  datTab.filtered <- datTab %>%
    dplyr::filter(.data$Score.Diff >= bestPreFilter)

  if (nrow(datTab.filtered) == 0) {
    warning(
      "Selected Score.Diff prefilter produced an empty data set. ",
      "Using unfiltered data.",
      call. = FALSE
    )

    bestPreFilter <- -Inf
    datTab.filtered <- datTab
  }

  list(
    datTab = datTab.filtered,
    bestPreFilter = bestPreFilter,
    preFilter.summary = preFilter.summary,
    prefilter.results = prefilter.results
  )
}


#' Performs hyperparamter optimziation for SVM model building by calling `tuneSVM.helper()`
#' across the specified grid of cost and gamma values.
#'
#' @param datTab Parsed CLMS search results
#' @param params Character vector specifying names of the features in `datTab` used to train model.
#' @param scoreName Name for the new scoring function.
#' @param scalingFactor An integer k. The multiple by which decoy DB is larger than target DB
#' @param targetER Desired FDR for classification of CSMs
#' @param sampleNo Size of the training dataset (integer).
#' @param cost_values Numeric vector of cost values used for hyperparameter tuning of the SVM model
#' @param gamma_values Numeric vector of gamma values used only for radial kernels.
#' @param kernels Character vector of kernels to evaluate. Defaults to
#'   `"linear"`. Include `"radial"` for experimental comparison.
#' @param seed Integer seed used to make cross-fitting reproducible.
#' @param splitBy Character vector naming columns whose rows must remain together
#'   during cross-fitting.
#' @param verbose Print training diagnostics.
#' @seealso [trainCrosslinkScore()], [tuneSVM.helper()], [buildSVM()]
#' @returns A list containing all of the SVM models at different cost and gamma values.
#' @export
tuneSVM <- function(datTab,
                    params = params.best,
                    scoreName="SVM.score",
                    scalingFactor = the$decoyScalingFactor,
                    targetER = 0.01,
                    sampleNo = 20000,
                    cost_values = c(0.5, 1, 5, 10),
                    gamma_values = c(0.01, 0.05, 0.1, 0.5),
                    kernels = "linear",
                    seed = 1,
                    splitBy = NULL,
                    verbose = FALSE) {
  param_grid <- makeSVMParameterGrid(
    cost_values = cost_values,
    gamma_values = gamma_values,
    kernels = kernels
  )

  tuned <- purrr::pmap(param_grid, function(cost, gamma, kernel) {
    tuneSVM.helper(datTab=datTab,
                   targetER=targetER,
                   params=params,
                   scoreName=scoreName,
                   scalingFactor = scalingFactor,
                   sampleNo = sampleNo,
                   cost, gamma, kernel,
                   seed = seed,
                   splitBy = splitBy,
                   verbose = verbose)
  })
  return(tuned)
}

makeSVMParameterGrid <- function(cost_values,
                                 gamma_values,
                                 kernels = "linear") {
  kernels <- match.arg(kernels, c("linear", "radial"), several.ok = TRUE)

  param_grid <- tibble::tibble()
  if ("linear" %in% kernels) {
    param_grid <- dplyr::bind_rows(
      param_grid,
      tibble::tibble(
        cost = cost_values,
        gamma = NA_real_,
        kernel = "linear"
      )
    )
  }
  if ("radial" %in% kernels) {
    param_grid <- dplyr::bind_rows(
      param_grid,
      tidyr::expand_grid(
        cost = cost_values,
        gamma = gamma_values,
        kernel = "radial"
      )
    )
  }

  param_grid
}

#' Helper function called by `tuneSVM()` that in turn calls `buildSVM()` and
#' calculates intermediate values in determining the objective functions used
#' for optimization. Namely, the errorTable generated by `generateErrorTable.sep()`
#' that reports number of interProtein and intraProtein crosslinks by FDR threshold.
#'
#' @param datTab Parsed CLMS search results
#' @param params Character vector specifying names of the features in `datTab` used to train model.
#' @param scoreName Name for the new scoring function.
#' @param scalingFactor An integer k. The multiple by which decoy DB is larger than target DB
#' @param targetER Desired FDR for classification of URPss
#' @param sampleNo Size of the training dataset (integer).
#' @param cost Cost value passed to `e1071::svm()`
#' @param gamma Gamma value passed to `e1071::svm()`
#' @param kernel Kernel value passed to `e1071::svm()`
#' @param seed Integer seed used to make cross-fitting reproducible.
#' @param splitBy Character vector naming columns whose rows must remain together
#'   during cross-fitting.
#' @param verbose Print training diagnostics.
#' @seealso [trainCrosslinkScore()], [tuneSVM()], [buildSVM()]
#' @returns A list containing the trained data at CSM and URP levels, score thresholds
#' for the targetER, the error table and some other information used for tuning.
#' @export
tuneSVM.helper <- function(datTab,
                           targetER=0.01,
                           params=params.best,
                           scoreName="SVM.score",
                           scalingFactor = the$decoyScalingFactor,
                           sampleNo = 20000,
                           cost = 1,
                           gamma = NA_real_,
                           kernel = "linear",
                           seed = 1,
                           splitBy = NULL,
                           verbose = FALSE) {
  kernel <- match.arg(kernel, c("linear", "radial"))
  if (kernel == "radial" &&
      (length(gamma) != 1 || !is.finite(gamma) || gamma <= 0)) {
    stop("gamma must be positive and finite for a radial SVM.", call. = FALSE)
  }
  svm.args <- list(
    datTab = datTab,
    params = params,
    scoreName = scoreName,
    sampleNo = sampleNo,
    showTab = FALSE,
    seed = seed,
    splitBy = splitBy,
    verbose = verbose,
    cost = cost,
    kernel = kernel
  )
  if (kernel == "radial") {
    svm.args$gamma <- gamma
  }
  datTab.csm <- do.call(buildSVM, svm.args)
  datTab.urp <- bestResPair(datTab.csm)
  datTab.urp.thresh <- findSeparateThresholdsModelled(datTab.urp,
                                                      targetER = targetER,
                                                      scalingFactor = scalingFactor,
                                                      plot = F)
  numHits <- classifyDataset(datTab.urp, datTab.urp.thresh) %>%
    removeDecoys() %>%
    count(.data$xlinkClass)
  intraHits = numHits[numHits$xlinkClass=="intraProtein", "n"][[1]]
  interHits = numHits[numHits$xlinkClass=="interProtein", "n"][[1]]
  if (length(intraHits)==0) {intraHits <- 0}
  if (length(interHits)==0) {interHits <- 0}
  errorTable <- generateErrorTable.sep(datTab.urp)
  inter.integral <- errorTable %>%
    filter(dplyr::between(.data$fdr.inter, 0.01, 0.05)) %>%
    summarize(inter.sum = sum(.data$inter), n= n(), inter.int = .data$inter.sum / n) %>%
    pull(.data$inter.int)
  top.inter.csms <- datTab.csm %>%
    filter(.data$Decoy == "Target",
           .data$xlinkClass == "interProtein") %>%
    arrange(desc(.data$Score.Diff))
  top.inter.csms <- dplyr::slice_head(
    top.inter.csms,
    n = floor(nrow(top.inter.csms) / 2)
  )
  correlation_score <- if (nrow(top.inter.csms) >= 2) {
    100 * stats::cor(
      top.inter.csms$Score.Diff,
      top.inter.csms$SVM.score,
      method = "spearman"
    )
  } else {
    NA_real_
  }

  list("CSMs" = datTab.csm,
       "URPs" = datTab.urp,
       "thresh" = datTab.urp.thresh,
       "intraHits" = intraHits,
       "interHits" = interHits,
       "errorTable" = errorTable,
       "interInt" = inter.integral,
       "corScore" = correlation_score,
       "cost" = cost,
       "gamma" = gamma,
       "kernel" = kernel,
       "sd.thresh" = min(datTab.csm$Score.Diff),
       "params" = params)
}

#' Basic function to build a new SVM classifier.  Doesn't do any feature selection or
#' hyperparamter tuning. Builds two separate SVM models on non-overlapping groups
#' of the data and averages their out-of-training-group predictions.
#'
#' @param datTab Parsed CLMS search results.
#' @param params Character vector specifying names of the features in `datTab` used to train model.
#' @param scoreName Name for the new scoring function.
#' @param sampleNo Target maximum size of each training subset (integer).
#'   A subset can exceed this target when necessary to keep a group intact or
#'   retain both outcome classes.
#' @param showTab print classificaiton table?
#' @param seed Integer seed used to make cross-fitting reproducible. Use `NULL`
#'   to use R's current random-number state.
#' @param splitBy Character vector naming columns whose rows must remain together
#'   during cross-fitting. The default uses `xlinkedResPair` when present, then
#'   a spectrum identifier, and finally individual rows.
#' @param verbose Print training-data diagnostics.
#' @param ... paramters passed to `e1071:svm()` function
#' @seealso [trainCrosslinkScore()], [tuneSVM.helper()], [tuneSVM()]
#' @return A data frame, one column larger than the input containing the new score.
#' @export
buildSVM <- function(datTab,
                     params=params.best,
                     scoreName="SVM.score",
                     sampleNo = 20000,
                     showTab = F,
                     seed = 1,
                     splitBy = NULL,
                     verbose = FALSE,
                     ...) {
  datTab$massError <- abs(datTab$ppm - mean(datTab$ppm))
  split <- makeCrossfitSplit(
    datTab,
    sampleNo = sampleNo,
    splitBy = splitBy,
    seed = seed
  )
  ind.1 <- split$train.1
  ind.2 <- split$train.2
  train.1 <- datTab[ind.1,]
  train.2 <- datTab[ind.2,]
  wghts.1 <- numeric(0)
  wghts.2 <- numeric(0)
  wghts.1["Target"] <- table(train.1$Decoy2)["Decoy"] / sum(table(train.1$Decoy2),na.rm=T)
  wghts.1["Decoy"] <- table(train.1$Decoy2)["Target"] / sum(table(train.1$Decoy2),na.rm=T)
  wghts.2["Target"] <- table(train.2$Decoy2)["Decoy"] / sum(table(train.2$Decoy2),na.rm=T)
  wghts.2["Decoy"] <- table(train.2$Decoy2)["Target"] / sum(table(train.2$Decoy2),na.rm=T)

  if (verbose) {
    diagnoseSVMdata(
      train.df = train.1,
      response.col = "Decoy2",
      feature.cols = params
    )

    diagnoseSVMdata(
      train.df = train.2,
      response.col = "Decoy2",
      feature.cols = params
    )
  }


  fit.1 <- e1071::svm(train.1$Decoy2 ~.,
                      subset(train.1, select=params),
                      class.weights=wghts.1,
                      ...
  )
  fit.2 <- e1071::svm(train.2$Decoy2 ~.,
                      subset(train.2, select=params),
                      class.weights=wghts.2,
                      ...
  )
  p.1 <- stats::predict(fit.1, subset(datTab, select=params),decision.values=T)
  p.2 <- stats::predict(fit.2, subset(datTab, select=params),decision.values=T)
  datTab$score.1 = as.numeric(attr(p.1, "decision.values"))
  datTab$score.2 = as.numeric(attr(p.2, "decision.values"))
  datTab[!split$score.1, "score.1"] <- NA
  datTab[!split$score.2, "score.2"] <- NA
  if (stats::cor(datTab$Score.Diff, datTab$score.1, use = "complete.obs") < 0) {
    datTab$score.1 <- -1 * datTab$score.1
  }
  if (stats::cor(datTab$Score.Diff, datTab$score.2, use = "complete.obs") < 0) {
    datTab$score.2 <- -1 * datTab$score.2
  }
  datTab[[scoreName]] <- purrr::map2_dbl(datTab$score.1, datTab$score.2, function(x, y) mean(c(x, y), na.rm=T))
  if (showTab) {
    tab <- table(datTab$Decoy2, datTab[[scoreName]] > 0)
    print(tab)
    print(paste("specificity:", round(tab[1]/(tab[1]+tab[3]),2)))
  }
  return(datTab)
}

makeCrossfitSplit <- function(datTab,
                              sampleNo = 20000,
                              splitBy = NULL,
                              seed = 1) {
  n <- nrow(datTab)

  if (n < 2) {
    stop("At least two rows are required for cross-fitting.", call. = FALSE)
  }

  if (length(sampleNo) != 1 || is.na(sampleNo) || sampleNo < 1) {
    stop("sampleNo must be one positive number.", call. = FALSE)
  }

  sampleNo <- as.integer(sampleNo)

  if (is.null(splitBy)) {
    splitBy <- dplyr::case_when(
      "xlinkedResPair" %in% names(datTab) ~ list("xlinkedResPair"),
      all(c("Fraction", "Spectrum") %in% names(datTab)) ~
        list(c("Fraction", "Spectrum")),
      all(c("Fraction", "MSMS.Info") %in% names(datTab)) ~
        list(c("Fraction", "MSMS.Info")),
      "Spectrum" %in% names(datTab) ~ list("Spectrum"),
      "MSMS.Info" %in% names(datTab) ~ list("MSMS.Info"),
      TRUE ~ list(character())
    )[[1]]
  }

  missing.split.columns <- setdiff(splitBy, names(datTab))
  if (length(missing.split.columns) > 0) {
    stop(
      "Cross-fitting group column(s) not found: ",
      paste(missing.split.columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (length(splitBy) == 0) {
    group.id <- as.character(seq_len(n))
  } else {
    group.parts <- lapply(datTab[splitBy], function(x) {
      x <- as.character(x)
      encoded <- paste0(nchar(enc2utf8(x), type = "bytes"), ":", x)
      encoded[is.na(x)] <- "-1:"
      encoded
    })
    group.id <- do.call(paste, c(group.parts, sep = "|"))
  }

  group.rows <- split(seq_len(n), group.id)

  if (length(group.rows) < 2) {
    stop(
      "Cross-fitting requires at least two distinct groups in splitBy.",
      call. = FALSE
    )
  }

  had.seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had.seed) {
    old.seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (!is.null(seed)) {
      if (had.seed) {
        assign(".Random.seed", old.seed, envir = .GlobalEnv)
      } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        rm(".Random.seed", envir = .GlobalEnv)
      }
    }
  }, add = TRUE)

  if (!is.null(seed)) {
    set.seed(seed)
  }

  fold.groups <- list(character(), character())
  fold.sizes <- c(0L, 0L)

  group.strata <- if ("Decoy2" %in% names(datTab)) {
    vapply(group.rows, function(rows) {
      paste(sort(unique(as.character(datTab$Decoy2[rows]))), collapse = "|")
    }, character(1))
  } else {
    stats::setNames(rep("all", length(group.rows)), names(group.rows))
  }

  strata <- split(names(group.rows), group.strata)

  if ("Decoy2" %in% names(datTab) && any(lengths(strata) < 2)) {
    sparse.strata <- names(strata)[lengths(strata) < 2]
    stop(
      "Cross-fitting requires at least two independent splitBy groups for ",
      "each outcome. Insufficient groups for: ",
      paste(sparse.strata, collapse = ", "),
      ".",
      call. = FALSE
    )
  }

  for (stratum.groups in strata) {
    shuffled.groups <- sample(stratum.groups, length(stratum.groups))
    stratum.fold.sizes <- c(0L, 0L)

    for (group in shuffled.groups) {
      smallest.stratum.folds <- which(
        stratum.fold.sizes == min(stratum.fold.sizes)
      )
      destination <- smallest.stratum.folds[
        which.min(fold.sizes[smallest.stratum.folds])
      ]
      fold.groups[[destination]] <- c(fold.groups[[destination]], group)
      group.size <- length(group.rows[[group]])
      fold.sizes[destination] <- fold.sizes[destination] + group.size
      stratum.fold.sizes[destination] <-
        stratum.fold.sizes[destination] + group.size
    }
  }

  limit.fold <- function(groups) {
    if (sum(lengths(group.rows[groups])) <= sampleNo) {
      return(groups)
    }

    groups.by.stratum <- split(groups, group.strata[groups])
    groups.by.stratum <- lapply(groups.by.stratum, function(x) {
      sample(x, length(x))
    })

    # Retain at least one independent group from every outcome class. This can
    # exceed sampleNo when a single group is unusually large, but avoids
    # creating an SVM training subset with a missing class.
    selected <- vapply(groups.by.stratum, `[[`, character(1), 1)
    selected.size <- sum(lengths(group.rows[selected]))
    remaining <- unlist(lapply(groups.by.stratum, function(x) x[-1]),
                        use.names = FALSE)

    if (length(remaining) > 1) {
      remaining <- sample(remaining, length(remaining))
    }

    for (group in remaining) {
      proposed.size <- selected.size + length(group.rows[[group]])
      if (abs(sampleNo - proposed.size) <= abs(sampleNo - selected.size)) {
        selected <- c(selected, group)
        selected.size <- proposed.size
      }
    }

    selected
  }

  fold.groups <- lapply(fold.groups, limit.fold)
  train.1 <- unlist(group.rows[fold.groups[[1]]], use.names = FALSE)
  train.2 <- unlist(group.rows[fold.groups[[2]]], use.names = FALSE)

  score.1 <- !group.id %in% fold.groups[[1]]
  score.2 <- !group.id %in% fold.groups[[2]]

  list(
    train.1 = sort(train.1),
    train.2 = sort(train.2),
    score.1 = score.1,
    score.2 = score.2,
    group.id = group.id,
    splitBy = splitBy
  )
}

#' Automated function to select features, build SVM score, and perform hyper-parameter tuning
#'
#' replaced by buildCrosslinkScore but keeping this function around for backward capability.
#' `trainClassifier()` is the main function that most users will use to automatically
#' used for model training and then calls `buildClassifier()` at different pre-filter
#' values and selects the best performing model.
#'
#' @param datTab Parsed CLMS search results
#' @param params Character vector specifying names of the features in `datTab` used to train model.
#' @param scoreName Name for the new scoring function.
#' @param scalingFactor An integer k. The multiple by which decoy DB is larger than target DB
#' @param targetER Desired FDR for classification of CSMs
#' @param preFilterER.values Numeric vector of preFilter FDR values for hyperparamter tuning
#' @return A list object containing the preFiltered data with the SVM.score as well as other information about the training procedure.
#' @seealso [trainClassifier_parallel()], [buildClassifier]
#' @export
#'
trainClassifier <- function(datTab, params=NULL, scoreName="SVM.score",
                            scalingFactor = the$decoyScalingFactor, targetER = 0.01,
                            preFilterER.values = c(0.45, 0.35, 0.25)) {
  #  start.time <- Sys.time()

  datTab <- ungroup(datTab)
  # feature selection
  plausibleHits <- datTab %>%
    filter(.data$Decoy == "Target", .data$Score.Diff > 10, .data$numCSM > 1, .data$xlinkClass == "intraProtein") %>%
    group_by(.data$Acc.1) %>%
    count %>%
    arrange(desc(.data$n))
  if (is.null(params)) {
    params <- case_when(
      nrow(plausibleHits) <= 50 ~ list(params.best.nop),
      plausibleHits$n[1] > 100 * plausibleHits$n[2] ~ list(params.best.nop),
      TRUE ~ list(params.best)
    ) %>% unlist()
  }
  # test model at different pre-filter error rates:
  baseFDR <- calculateFDR(datTab, classifier = "Score.Diff", scalingFactor = scalingFactor)
  preFilterER.values <- c(baseFDR, preFilterER.values[baseFDR > preFilterER.values])
  resultsList <- preFilterER.values %>%
    map(function(x) buildClassifier(datTab=datTab, params=params,
                                    scoreName=scoreName, scalingFactor=scalingFactor,
                                    preFilterER = x))
  resultsProts <- resultsList %>%
    map(bestProtPair)
  thresholds <- map(resultsProts, function(datTab)
    findSeparateThresholdsModelled(datTab, targetER = targetER,
                                   scalingFactor = scalingFactor))
  classifiedList <- purrr::map2(resultsProts, thresholds, function(datTab, thresh)
    classifySeparateThresholds(datTab, separateThresh = thresh))
  numInterHits <- purrr::map_dbl(classifiedList, function(classedTable)
    classedTable %>%
      filter(.data$Decoy == "Target", .data$xlinkClass == "interProtein") %>%
      count %>% pull)
  bestPreFilter <- which.max(numInterHits)
  scoreThreshold.csm <- findSeparateThresholdsModelled(resultsList[[bestPreFilter]],
                                                       targetER = targetER,
                                                       scalingFactor = scalingFactor)
  # end.time <- Sys.time()
  # time.taken <- end.time - start.time
  # print(time.taken)

  return(list(
    "params" = params,
    "preFilterER.best" = preFilterER.values[[bestPreFilter]],
    "preFilter.SD.min" = min(resultsList[[bestPreFilter]]$Score.Diff, na.rm=T),
    "dataTable" = resultsList[[bestPreFilter]],
    "scoreThreshold" = scoreThreshold.csm,
    "preFilterERs.all" = preFilterER.values,
    "numInterHits.all" = numInterHits
  ))
}

#' Automated function to select features, build SVM score, and perform hyper-paramter tuning
#'
#' `trainClassifier_parallel()` is the parallelized version `trainClassifier()`.
#' Requires that `furrr` package to be installed and user must select an appropriate
#' `furrr::plan()`. In testing with large datasets, the parallelized computation
#' does not results in significant improvements in processing time. This function is
#' still under development.
#'
#' @param datTab Parsed CLMS search results
#' @param params Character vector specifying names of the features in `datTab` used to train model.
#' @param scoreName Name for the new scoring function.
#' @param scalingFactor An integer k. The multiple by which decoy DB is larger than target DB
#' @param targetER Desired FDR for classification of CSMs
#' @param preFilterER.values Numeric vector of preFilter FDR values for hyperparamter tuning
#' @return A list object containing the preFiltered data with the SVM.score as well as other information about the training procedure.
#' @seealso [trainClassifier_parallel()], [buildClassifier]
#' @export
#'
trainClassifier_parallel <- function(datTab, params=NA, scoreName="SVM.score",
                                     scalingFactor = the$decoyScalingFactor, targetER = 0.01,
                                     preFilterER.values = c(0.45, 0.35, 0.25)) {
  if (!requireNamespace("furrr", quietly = TRUE) ||
      !requireNamespace("future", quietly = TRUE)) {
    stop(
      "Packages 'furrr' and 'future' are required for parallel training.",
      call. = FALSE
    )
  }
  # start.time = Sys.time()

  oopts <- options(future.globals.maxSize = 8000 * 1024^2)
  on.exit(options(oopts))

  datTab <- ungroup(datTab)
  # feature selection
  plausibleHits <- datTab %>%
    filter(.data$Decoy == "Target", .data$Score.Diff > 10, .data$numCSM > 1, .data$xlinkClass == "intraProtein") %>%
    group_by(.data$Acc.1) %>%
    count %>%
    arrange(desc(.data$n))
  if (is.na(params)) {
    params <- case_when(
      nrow(plausibleHits) <= 1 ~ list(params.best.nop),
      plausibleHits$n[1] > 100 * plausibleHits$n[2] ~ list(params.best.nop),
      TRUE ~ list(params.best)
    ) %>% unlist()
  }
  # test model at different pre-filter error rates:
  baseFDR <- calculateFDR(datTab, classifier = "Score.Diff", scalingFactor = scalingFactor)
  preFilterER.values <- c(baseFDR, preFilterER.values[baseFDR > preFilterER.values])

  # time.1 <- Sys.time()

  resultsList <- preFilterER.values %>%
    furrr::future_map(function(x) buildClassifier(datTab=datTab, params=params,
                                                  scoreName=scoreName, scalingFactor=scalingFactor,
                                                  preFilterER = x),
                      .progress = T, .options = furrr::furrr_options(seed = T)
    )
  resultsProts <- resultsList %>%
    furrr::future_map(bestProtPair)
  thresholds <- furrr::future_map(resultsProts, function(datTab)
    findSeparateThresholdsModelled(datTab, targetER = targetER,
                                   scalingFactor = scalingFactor),
    .progress=T, .options = furrr::furrr_options(seed = T))

  # time.2 <- Sys.time()

  classifiedList <- purrr::map2(resultsProts, thresholds, function(datTab, thresh)
    classifySeparateThresholds(datTab, separateThresh = thresh))
  numInterHits <- purrr::map_dbl(classifiedList, function(classedTable)
    classedTable %>%
      filter(.data$Decoy == "Target", .data$xlinkClass == "interProtein") %>%
      count %>% pull)
  bestPreFilter <- which.max(numInterHits)
  scoreThreshold.csm <- findSeparateThresholdsModelled(resultsList[[bestPreFilter]],
                                                       targetER = targetER,
                                                       scalingFactor = scalingFactor)
  # end.time <- Sys.time()
  # time.taken <- c(time.1, time.2, end.time) - start.time
  # print(time.taken)

  return(list(
    "params" = params,
    "preFilterER.best" = preFilterER.values[[bestPreFilter]],
    "preFilter.SD.min" = min(resultsList[[bestPreFilter]]$Score.Diff, na.rm=T),
    "dataTable" = resultsList[[bestPreFilter]],
    "scoreThreshold" = scoreThreshold.csm,
    "preFilterERs.all" = preFilterER.values,
    "numInterHits.all" = numInterHits
  ))
}

#' Train an SVM scoring function for CLMS data
#'
#' replaced by buildSVM, but keeping it around for backward capability.
#' `buildClassifier()` is the basic function used to create an SVM scoring model
#' for CLMS results. Useful for testing different sets of parameters. Results
#' are pre-filtered to a specified FDR (based on Score.Diff) prior to training.
#'
#' @param datTab Parsed CLMS search results.
#' @param params Character vector specifying names of the features in `datTab` used to train model.
#' @param preFilterER Desired error rate (FDR) for pre-filtering.
#' @param scoreName Name for the new scoring function.
#' @param scalingFactor An integer k. he multiple by which decoy DB is larger than target DB
#' @param ... paramters passed down to `e1071::svm()`
#' @param sampleNo Size of the training dataset (integer).
#' @return A data frame, one column larger than the input containing the new score.
#' @seealso [trainClassifier()]
#' @export
#'

buildClassifier <- function(datTab, params=params.best, preFilterER = NA,
                            scoreName="SVM.score", scalingFactor = the$decoyScalingFactor,
                            sampleNo = 20000, ...) {
  baseFDR <- calculateFDR(datTab, classifier = "Score.Diff", scalingFactor = scalingFactor)
  if (!is.na(preFilterER) & baseFDR > preFilterER) {
    if (nrow(datTab) > sampleNo) {
      datTab.pre <- slice(datTab, sample(nrow(datTab), sampleNo))
    } else {
      datTab.pre <- datTab
    }
    preFilter.thresh <- findSeparateThresholds(datTab.pre,
                                               targetER = preFilterER,
                                               minThreshold = -5,
                                               classifier = "Score.Diff",
                                               errorFUN = calculateFDR.unseparated,
                                               scalingFactor = scalingFactor)
    datTab <- classifySeparateThresholds(datTab, preFilter.thresh, classifier="Score.Diff")
  }
  buildSVM(datTab, params, scoreName, sampleNo, ...)
}


diagnoseSVMdata <- function(train.df, response.col, feature.cols) {
  if (!response.col %in% colnames(train.df)) {
    stop("response.col not found: ", response.col, call. = FALSE)
  }

  missing.features <- setdiff(feature.cols, colnames(train.df))
  if (length(missing.features) > 0) {
    stop(
      "Feature column(s) not found: ",
      paste(missing.features, collapse = ", "),
      call. = FALSE
    )
  }

  message("\nOutcome counts:")
  print(table(train.df[[response.col]], useNA = "ifany"))

  y.n <- length(unique(stats::na.omit(train.df[[response.col]])))
  if (y.n < 2) {
    message("PROBLEM: response variable has fewer than 2 non-NA classes.")
  }

  feature.summary <- purrr::map_dfr(feature.cols, function(col) {
    x <- train.df[[col]]

    tibble::tibble(
      feature = col,
      class = paste(class(x), collapse = "/"),
      n_non_missing = sum(!is.na(x)),
      n_unique_non_missing = length(unique(stats::na.omit(x))),
      n_levels_after_factor = if (is.factor(x) || is.character(x)) {
        nlevels(factor(stats::na.omit(x)))
      } else {
        NA_integer_
      },
      example_values = paste(utils::head(unique(stats::na.omit(x)), 5), collapse = ", ")
    )
  })

  message("\nFeature summary:")
  print(feature.summary, n = Inf)

  bad.categorical <- feature.summary |>
    dplyr::filter(
      .data$class %in% c("factor", "character") |
        stringr::str_detect(.data$class, "factor|character")
    ) |>
    dplyr::filter(.data$n_levels_after_factor < 2)

  bad.constant <- feature.summary |>
    dplyr::filter(.data$n_unique_non_missing < 2)

  if (nrow(bad.categorical) > 0) {
    message("\nCategorical predictors with fewer than 2 levels:")
    print(bad.categorical, n = Inf)
  }

  if (nrow(bad.constant) > 0) {
    message("\nConstant / all-NA predictors:")
    print(bad.constant, n = Inf)
  }

  invisible(feature.summary)
}


check_training_df <- function(df, stage, response_col, feature_cols = NULL) {
  message("\n--- ", stage, " ---")
  message("nrow: ", nrow(df))
  message("ncol: ", ncol(df))

  if (response_col %in% colnames(df)) {
    message("response counts:")
    print(table(df[[response_col]], useNA = "ifany"))
  } else {
    message("response column missing: ", response_col)
  }

  if (!is.null(feature_cols)) {
    missing_features <- setdiff(feature_cols, colnames(df))
    if (length(missing_features) > 0) {
      message("missing features: ", paste(missing_features, collapse = ", "))
    }

    present_features <- intersect(feature_cols, colnames(df))

    feature_summary <- purrr::map_dfr(present_features, function(col) {
      x <- df[[col]]
      tibble::tibble(
        feature = col,
        class = paste(class(x), collapse = "/"),
        n_non_missing = sum(!is.na(x)),
        n_unique = length(unique(stats::na.omit(x))),
        values = paste(utils::head(unique(stats::na.omit(x)), 5), collapse = ", ")
      )
    })

    print(feature_summary, n = Inf)
  }

  invisible(df)
}
