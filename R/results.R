#' Prepare crosslink results at a reporting level
#'
#' Starts from a scored CSM table, summarizes it at the requested proteomics
#' level, calculates level-specific inter- and intra-protein thresholds, and
#' generates a scaling-adjusted classification summary with [countDecoys()]. The
#' complete summarized target-and-decoy table and its thresholds are retained as
#' the canonical result. Classified and target-only tables can be generated
#' without storing duplicate copies; see Examples.
#'
#' When `x` is a result from [trainCrosslinkScore()], `model = "selected"` uses
#' the linear recommendation. `model = "radial"` uses `recommendedRadial`, and
#' a numeric value uses that candidate index. When
#' the selected candidate's URP table and thresholds already match the requested
#' settings, those cached results are reused.
#'
#' @param x A `touchstone_training` result or a CSM-level data frame containing
#'   the requested `classifier` column.
#' @param summarizationLevel Reporting level. One of `"csm"`, `"urp"`,
#'   `"peptide-pair"`, `"protein-pair"`, or `"module-pair"`.
#' @param model For a training result, `"selected"` or `"linear"` for the
#'   recommended linear model, `"radial"` for the recommended radial model, or
#'   a numeric candidate index. Ignored for a CSM data frame.
#' @param targetER Desired FDR. Defaults to the value stored in a training result,
#'   or 0.01 for a data frame.
#' @param scalingFactor Multiple by which the decoy database is larger than the
#'   target database. Defaults to the value stored in a training result, or the
#'   value established by [setDecoyScalingFactor()] for a data frame.
#' @param retainGroups Retain existing data-frame groups during pair
#'   summarization; passed to the existing `best*Pair()` functions.
#' @param classifier Column in the CSM table used to rank, summarize, threshold,
#'   and classify hits. Defaults to `"SVM.score"`; alternatives such as
#'   `"Score.Diff"` or an experimental score column are supported.
#' @return A `touchstone_results` object containing the complete summarized
#'   `data`, level-specific `thresholds`, the resulting `fdr`, a compact
#'   `classificationSummary`, and the settings used. Use [classifyDataset()] and
#'   [removeDecoys()] to generate filtered views.
#' @examples
#' \dontrun{
#' urp.results <- prepareCrosslinkResults(training, "urp")
#' accepted <- classifyDataset(urp.results$data, urp.results$thresholds)
#' reported <- removeDecoys(accepted)
#'
#' radial.pp <- prepareCrosslinkResults(
#'   training,
#'   summarizationLevel = "protein-pair",
#'   model = "radial"
#' )
#' }
#' @export
prepareCrosslinkResults <- function(x,
                                    summarizationLevel = "urp",
                                    model = "selected",
                                    targetER = NULL,
                                    scalingFactor = NULL,
                                    retainGroups = TRUE,
                                    classifier = "SVM.score") {
  level.aliases <- c(
    "csm" = "csm",
    "spectra" = "csm",
    "urp" = "urp",
    "residue-pair" = "urp",
    "peptide-pair" = "peptide-pair",
    "protein-pair" = "protein-pair",
    "module-pair" = "module-pair"
  )

  if (length(summarizationLevel) != 1 ||
      !summarizationLevel %in% names(level.aliases)) {
    stop(
      "summarizationLevel must be one of: csm, urp, peptide-pair, ",
      "protein-pair, or module-pair.",
      call. = FALSE
    )
  }
  summarizationLevel <- unname(level.aliases[[summarizationLevel]])

  classifier <- rlang::as_name(rlang::ensym(classifier))

  resolved <- resolveCrosslinkFit(x, model)
  scored.csms <- resolved$fit$CSMs
  if (is.null(scored.csms)) {
    scored.csms <- resolved$fit$scoredCSMs
  }
  required <- c(classifier, "xlinkClass", "Decoy")
  missing.columns <- setdiff(required, names(scored.csms))
  if (length(missing.columns) > 0) {
    if (classifier %in% missing.columns) {
      stop(
        "The scored CSM data have no classifier column named '", classifier,
        "'. Supply an existing score column with classifier, or pass data ",
        "scored with that classifier.",
        call. = FALSE
      )
    }
    stop(
      "Scored CSM data are missing required column(s): ",
      paste(missing.columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(targetER)) {
    targetER <- resolved$settings$targetER
  }
  if (is.null(targetER)) {
    targetER <- 0.01
  }
  if (length(targetER) != 1 || !is.finite(targetER) ||
      targetER <= 0 || targetER >= 1) {
    stop("targetER must be one finite number between 0 and 1.", call. = FALSE)
  }

  if (is.null(scalingFactor)) {
    scalingFactor <- resolved$settings$scalingFactor
  }
  if (is.null(scalingFactor)) {
    scalingFactor <- the$decoyScalingFactor
  }
  if (length(scalingFactor) != 1 || !is.finite(scalingFactor) ||
      scalingFactor <= 0) {
    stop("scalingFactor must be one positive, finite number.", call. = FALSE)
  }

  cached.classifier <- resolved$settings$scoreName
  if (is.null(cached.classifier)) {
    cached.classifier <- "SVM.score"
  }
  use.cached.urp <- summarizationLevel == "urp" && retainGroups &&
    identical(classifier, cached.classifier) &&
    !is.null(resolved$fit$URPs)

  summarized <- if (use.cached.urp) {
    resolved$fit$URPs
  } else {
    summarize.with <- function(fun) {
      do.call(
        fun,
        list(
          datTab = scored.csms,
          classifier = classifier,
          retainGroups = retainGroups
        )
      )
    }
    switch(
      summarizationLevel,
      csm = scored.csms,
      urp = summarize.with(bestResPair),
      `peptide-pair` = summarize.with(bestPepPair),
      `protein-pair` = summarize.with(bestProtPair),
      `module-pair` = summarize.with(bestModPair)
    )
  }

  use.cached.thresholds <- use.cached.urp &&
    !is.null(resolved$fit$thresh) &&
    identical(as.numeric(targetER), as.numeric(resolved$settings$targetER)) &&
    identical(
      as.numeric(scalingFactor),
      as.numeric(resolved$settings$scalingFactor)
    )

  thresholds <- if (use.cached.thresholds) {
    resolved$fit$thresh
  } else {
    findSeparateThresholdsModelled(
      summarized,
      targetER = targetER,
      scalingFactor = scalingFactor,
      plot = FALSE,
      classifier = classifier
    )
  }

  classification.summary <- countDecoys(
    summarized,
    threshold = thresholds,
    classifier = classifier,
    scalingFactor = scalingFactor
  )
  if (!is.data.frame(classification.summary)) {
    stop("countDecoys() did not return a data frame.", call. = FALSE)
  }
  fdr <- calculateFDR(
    summarized,
    threshold = thresholds,
    classifier = classifier,
    scalingFactor = scalingFactor
  )

  structure(
    list(
      data = summarized,
      thresholds = thresholds,
      fdr = fdr,
      classificationSummary = classification.summary,
      summarizationLevel = summarizationLevel,
      model = resolved$model,
      settings = list(
        targetER = targetER,
        scalingFactor = scalingFactor,
        retainGroups = retainGroups,
        classifier = classifier
      )
    ),
    class = "touchstone_results"
  )
}

resolveCrosslinkFit <- function(x, model = "selected") {
  if (is.data.frame(x)) {
    return(list(
      fit = list(CSMs = x),
      model = list(requested = "scored-csms"),
      settings = list()
    ))
  }
  if (!inherits(x, "touchstone_training")) {
    stop("x must be a trainCrosslinkScore() result or scored CSM data frame.",
         call. = FALSE)
  }

  if (is.numeric(model)) {
    if (length(model) != 1 || !is.finite(model) || model != as.integer(model) ||
        model < 1 || model > length(x$models)) {
      stop("Numeric model must be a valid candidate index.", call. = FALSE)
    }
    index <- as.integer(model)
    fit <- x$models[[index]]
    requested <- paste0("candidate-", index)
  } else {
    if (length(model) != 1 || !model %in% c("selected", "linear", "radial")) {
      stop("model must be selected, linear, radial, or a candidate index.",
           call. = FALSE)
    }
    requested <- model
    if (model %in% c("selected", "linear")) {
      fit <- x$recommended
      index <- x$candidates$index[x$candidates$recommended]
    } else {
      fit <- x$recommendedRadial
      index <- x$candidates$index[x$candidates$recommendedRadial]
    }
  }

  if (is.null(fit) || length(index) == 0 || is.na(index)) {
    stop("The requested model was not trained or had no eligible candidate.",
         call. = FALSE)
  }

  list(
    fit = fit,
    model = list(
      requested = requested,
      index = index,
      kernel = fit$kernel,
      cost = fit$cost,
      gamma = fit$gamma
    ),
    settings = x$settings
  )
}

#' Print prepared crosslink results
#'
#' @param x A result returned by [prepareCrosslinkResults()].
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.touchstone_results <- function(x, ...) {
  cat(
    "Touchstone results: ", nrow(x$data), " ", x$summarizationLevel,
    " rows; target FDR ", format(x$settings$targetER),
    "; calculated FDR ", format(x$fdr), ".\n",
    sep = ""
  )
  print(x$classificationSummary)
  invisible(x)
}
