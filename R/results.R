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
#' @param thresholds Optional manually selected threshold: one numeric global
#'   threshold, a list containing `globalThresh`, or a list containing both
#'   `interThresh` and `intraThresh`. When omitted, thresholds are estimated
#'   with [findSeparateThresholdsModelled()].
#' @return A `touchstone_results` object containing the complete summarized
#'   `data`, its complete scored CSM source, level-specific `thresholds`, the
#'   resulting `fdr`, a compact `classificationSummary`, and the settings used.
#'   Use [classifyCrosslinkResults()] to generate a classified and optionally
#'   polished reporting table with recalculated support counts.
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
                                    classifier = "SVM.score",
                                    thresholds = NULL) {
  summarizationLevel <- normalizeSummarizationLevel(summarizationLevel)

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
    summarizeCrosslinkData(
      scored.csms,
      summarizationLevel = summarizationLevel,
      classifier = classifier,
      retainGroups = retainGroups
    )
  }

  use.cached.thresholds <- use.cached.urp &&
    !is.null(resolved$fit$thresh) &&
    identical(as.numeric(targetER), as.numeric(resolved$settings$targetER)) &&
    identical(
      as.numeric(scalingFactor),
      as.numeric(resolved$settings$scalingFactor)
    )

  threshold.source <- if (!is.null(thresholds)) "manual" else "modelled"
  thresholds <- if (!is.null(thresholds)) {
    validateCrosslinkThresholds(thresholds)
  } else if (use.cached.thresholds) {
    threshold.source <- "training-cache"
    validateCrosslinkThresholds(resolved$fit$thresh)
  } else {
    validateCrosslinkThresholds(findSeparateThresholdsModelled(
      summarized,
      targetER = targetER,
      scalingFactor = scalingFactor,
      plot = FALSE,
      classifier = classifier
    ))
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
      sourceCSMs = scored.csms,
      thresholds = thresholds,
      fdr = fdr,
      classificationSummary = classification.summary,
      summarizationLevel = summarizationLevel,
      stage = "prepared",
      model = resolved$model,
      settings = list(
        targetER = targetER,
        scalingFactor = scalingFactor,
        retainGroups = retainGroups,
        classifier = classifier,
        thresholdSource = threshold.source,
        polishing = NULL
      )
    ),
    class = "touchstone_results"
  )
}

normalizeSummarizationLevel <- function(summarizationLevel) {
  level.aliases <- c(
    "csm" = "csm",
    "spectra" = "csm",
    "urp" = "urp",
    "residue-pair" = "urp",
    "peptide-pair" = "peptide-pair",
    "protein-pair" = "protein-pair",
    "module-pair" = "module-pair"
  )
  if (length(summarizationLevel) != 1 || is.na(summarizationLevel) ||
      !summarizationLevel %in% names(level.aliases)) {
    stop(
      "summarizationLevel must be one of: csm, urp, peptide-pair, ",
      "protein-pair, or module-pair.",
      call. = FALSE
    )
  }
  unname(level.aliases[[summarizationLevel]])
}

summarizeCrosslinkData <- function(datTab,
                                   summarizationLevel,
                                   classifier,
                                   retainGroups) {
  summarize.with <- function(fun) {
    do.call(
      fun,
      list(
        datTab = datTab,
        classifier = classifier,
        retainGroups = retainGroups
      )
    )
  }
  switch(
    summarizationLevel,
    csm = datTab,
    urp = summarize.with(bestResPair),
    `peptide-pair` = summarize.with(bestPepPair),
    `protein-pair` = summarize.with(bestProtPair),
    `module-pair` = summarize.with(bestModPair)
  )
}

validateCrosslinkThresholds <- function(thresholds) {
  valid.number <- function(x) {
    is.numeric(x) && length(x) == 1 && !is.na(x) && is.finite(x)
  }
  if (valid.number(thresholds)) {
    return(thresholds)
  }
  if (is.list(thresholds) && valid.number(thresholds$globalThresh)) {
    return(thresholds)
  }
  if (is.list(thresholds) && valid.number(thresholds$interThresh) &&
      valid.number(thresholds$intraThresh)) {
    return(thresholds)
  }
  stop(
    "thresholds must be one finite numeric value, a list containing ",
    "globalThresh, or a list containing finite interThresh and intraThresh ",
    "values.",
    call. = FALSE
  )
}

#' Classify and optionally polish prepared crosslink results
#'
#' Applies a reporting threshold to the complete scored CSM source retained by
#' [prepareCrosslinkResults()], optionally applies transparent evidence filters,
#' recalculates support counts with [calculatePairs()], and only then summarizes
#' to the requested reporting level. This keeps `numCSM` and `numURP` aligned
#' with the evidence that actually passes the reporting policy.
#'
#' Polishing is supplied as a named list. Its names intentionally follow
#' [readProspectorXLOutput()]: `minPepLen`, `minPepScore`, `minScoreDiff`, and
#' `minIons`. The additional `minLadderCoverage` rule requires both peptides'
#' sequential product-ion ladder lengths to be at least that fraction of their
#' peptide lengths. Only explicitly supplied rules are considered. When the
#' columns needed by a requested rule were not selected in the Search Compare
#' output, that rule is skipped with a warning and recorded as unavailable in
#' `polishingAudit`.
#'
#' @param x A `touchstone_results` object returned by
#'   [prepareCrosslinkResults()].
#' @param thresholds Optional manual threshold. Defaults to the threshold stored
#'   in `x` and accepts the same forms as the `thresholds` argument to
#'   [prepareCrosslinkResults()].
#' @param polishing `NULL` for no evidence polishing, or a named list containing
#'   any of `minPepLen`, `minPepScore`, `minScoreDiff`, `minIons`, and
#'   `minLadderCoverage`.
#' @return A `touchstone_results` object whose `data` contain classified and
#'   optionally polished target and decoy results. The original complete scored
#'   CSM source is retained in `sourceCSMs`; `polishingAudit` reports whether
#'   each rule was applied and the number of CSMs it removed.
#' @examples
#' \dontrun{
#' prepared <- prepareCrosslinkResults(training, "urp")
#' classified <- classifyCrosslinkResults(prepared)
#' polished <- classifyCrosslinkResults(
#'   prepared,
#'   polishing = list(minIons = 3, minLadderCoverage = 0.25)
#' )
#' }
#' @export
classifyCrosslinkResults <- function(x, thresholds = NULL, polishing = NULL) {
  if (!inherits(x, "touchstone_results")) {
    stop("x must be a result returned by prepareCrosslinkResults().",
         call. = FALSE)
  }
  if (is.null(x$sourceCSMs) || !is.data.frame(x$sourceCSMs)) {
    stop(
      "This results object does not retain its scored CSM source. Recreate it ",
      "with prepareCrosslinkResults() before classification.",
      call. = FALSE
    )
  }
  if (is.null(thresholds)) {
    thresholds <- x$thresholds
    threshold.source <- x$settings$thresholdSource
    if (is.null(threshold.source)) threshold.source <- "prepared-result"
  } else {
    threshold.source <- "manual"
  }
  thresholds <- validateCrosslinkThresholds(thresholds)
  polishing <- normalizePolishingOptions(polishing)
  classifier <- x$settings$classifier
  if (is.null(classifier)) classifier <- "SVM.score"

  csms <- classifyDataset(
    x$sourceCSMs,
    threshold = thresholds,
    classifier = classifier
  )
  if (!is.data.frame(csms)) {
    stop("The supplied threshold could not be applied to the scored CSMs.",
         call. = FALSE)
  }

  audit <- tibble::tibble(
    rule = "threshold",
    value = formatThresholdForAudit(thresholds),
    before = nrow(x$sourceCSMs),
    after = nrow(csms),
    removed = nrow(x$sourceCSMs) - nrow(csms),
    applied = TRUE,
    reason = NA_character_
  )
  polished <- applyCrosslinkPolishing(csms, polishing)
  csms <- polished$data
  audit <- dplyr::bind_rows(audit, polished$audit)
  if (nrow(csms) == 0) {
    stop("No CSMs remain after thresholding and polishing.", call. = FALSE)
  }

  csms <- csms %>%
    dplyr::select(-dplyr::any_of(c("numCSM", "numURP", "wtCSM", "wtURP"))) %>%
    calculatePairs(scalingFactor = x$settings$scalingFactor) %>%
    dplyr::select(-dplyr::any_of(c("wtCSM", "wtURP")))
  summarized <- summarizeCrosslinkData(
    csms,
    summarizationLevel = x$summarizationLevel,
    classifier = classifier,
    retainGroups = x$settings$retainGroups
  ) %>%
    dplyr::select(-dplyr::any_of(c("wtCSM", "wtURP")))

  classification.summary <- countDecoys(
    summarized,
    threshold = thresholds,
    classifier = classifier,
    scalingFactor = x$settings$scalingFactor
  )
  fdr <- calculateFDR(
    summarized,
    threshold = thresholds,
    classifier = classifier,
    scalingFactor = x$settings$scalingFactor
  )
  result.settings <- x$settings
  result.settings$thresholdSource <- threshold.source
  result.settings$polishing <- polishing
  polishing.applied <- nrow(audit) > 1 && any(audit$applied[-1])

  structure(
    list(
      data = summarized,
      sourceCSMs = x$sourceCSMs,
      thresholds = thresholds,
      fdr = fdr,
      classificationSummary = classification.summary,
      summarizationLevel = x$summarizationLevel,
      stage = if (polishing.applied) "polished" else "classified",
      model = x$model,
      polishingAudit = audit,
      settings = result.settings
    ),
    class = "touchstone_results"
  )
}

normalizePolishingOptions <- function(polishing) {
  if (is.null(polishing)) return(list())
  if (!is.list(polishing) || is.null(names(polishing)) ||
      any(!nzchar(names(polishing)))) {
    stop("polishing must be NULL or a named list.", call. = FALSE)
  }
  allowed <- c(
    "minPepLen", "minPepScore", "minScoreDiff", "minIons",
    "minLadderCoverage"
  )
  unknown <- setdiff(names(polishing), allowed)
  if (length(unknown) > 0) {
    stop(
      "Unknown polishing option(s): ", paste(unknown, collapse = ", "),
      ". Allowed options are: ", paste(allowed, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (anyDuplicated(names(polishing))) {
    stop("Each polishing option may be supplied only once.", call. = FALSE)
  }
  for (name in names(polishing)) {
    value <- polishing[[name]]
    if (!is.numeric(value) || length(value) != 1 || is.na(value) ||
        !is.finite(value) || value < 0) {
      stop(name, " must be one non-negative, finite number.", call. = FALSE)
    }
  }
  if (!is.null(polishing$minLadderCoverage) &&
      polishing$minLadderCoverage > 1) {
    stop("minLadderCoverage must be between 0 and 1.", call. = FALSE)
  }
  polishing
}

applyCrosslinkPolishing <- function(datTab, polishing) {
  audit <- tibble::tibble(
    rule = character(), value = character(), before = integer(),
    after = integer(), removed = integer(), applied = logical(),
    reason = character()
  )
  apply.rule <- function(data, name, value, keep) {
    before <- nrow(data)
    keep[is.na(keep)] <- FALSE
    data <- data[keep, , drop = FALSE]
    audit <<- dplyr::bind_rows(
      audit,
      tibble::tibble(
        rule = name,
        value = format(value),
        before = before,
        after = nrow(data),
        removed = before - nrow(data),
        applied = TRUE,
        reason = NA_character_
      )
    )
    data
  }
  skip.rule <- function(data, name, value, missing.columns) {
    reason <- paste0(
      "Missing columns: ", paste(missing.columns, collapse = ", ")
    )
    warning(
      "Polishing rule '", name, "' was skipped because its required ",
      "column(s) are unavailable: ",
      paste(missing.columns, collapse = ", "), ".",
      call. = FALSE
    )
    audit <<- dplyr::bind_rows(
      audit,
      tibble::tibble(
        rule = name,
        value = format(value),
        before = nrow(data),
        after = nrow(data),
        removed = 0L,
        applied = FALSE,
        reason = reason
      )
    )
    data
  }
  apply.if.available <- function(data, name, value, columns, keep) {
    missing.columns <- setdiff(columns, names(data))
    if (length(missing.columns) > 0) {
      return(skip.rule(data, name, value, missing.columns))
    }
    apply.rule(data, name, value, keep)
  }

  if (!is.null(polishing$minPepLen)) {
    columns <- c("Len.Pep.1", "Len.Pep.2")
    keep <- if (all(columns %in% names(datTab))) {
      datTab$Len.Pep.1 >= polishing$minPepLen &
        datTab$Len.Pep.2 >= polishing$minPepLen
    } else logical()
    datTab <- apply.if.available(
      datTab, "minPepLen", polishing$minPepLen, columns, keep
    )
  }
  if (!is.null(polishing$minPepScore)) {
    columns <- c("Sc.1", "Sc.2")
    keep <- if (all(columns %in% names(datTab))) {
      datTab$Sc.1 >= polishing$minPepScore &
        datTab$Sc.2 >= polishing$minPepScore
    } else logical()
    datTab <- apply.if.available(
      datTab, "minPepScore", polishing$minPepScore, columns, keep
    )
  }
  if (!is.null(polishing$minScoreDiff)) {
    columns <- "Score.Diff"
    keep <- if (all(columns %in% names(datTab))) {
      datTab$Score.Diff >= polishing$minScoreDiff
    } else logical()
    datTab <- apply.if.available(
      datTab, "minScoreDiff", polishing$minScoreDiff, columns, keep
    )
  }
  if (!is.null(polishing$minIons)) {
    columns <- c("numProdIons.1", "numProdIons.2")
    keep <- if (all(columns %in% names(datTab))) {
      datTab$numProdIons.1 >= polishing$minIons &
        datTab$numProdIons.2 >= polishing$minIons
    } else logical()
    datTab <- apply.if.available(
      datTab, "minIons", polishing$minIons, columns, keep
    )
  }
  if (!is.null(polishing$minLadderCoverage)) {
    columns <- c("Len.Pep.1", "Len.Pep.2", "ladderLen.1", "ladderLen.2")
    missing.columns <- setdiff(columns, names(datTab))
    if (length(missing.columns) > 0) {
      datTab <- skip.rule(
        datTab, "minLadderCoverage", polishing$minLadderCoverage,
        missing.columns
      )
    } else {
      datTab <- apply.rule(
        datTab, "minLadderCoverage", polishing$minLadderCoverage,
        datTab$ladderLen.1 >=
          polishing$minLadderCoverage * datTab$Len.Pep.1 &
          datTab$ladderLen.2 >=
          polishing$minLadderCoverage * datTab$Len.Pep.2
      )
    }
  }
  list(data = datTab, audit = audit)
}

formatThresholdForAudit <- function(thresholds) {
  if (is.numeric(thresholds)) return(format(thresholds))
  paste(
    paste0(names(thresholds), "=", vapply(thresholds, format, character(1))),
    collapse = ", "
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
  stage <- x$stage
  if (is.null(stage)) stage <- "prepared"
  cat(
    "Touchstone ", stage, " results: ", nrow(x$data), " ",
    x$summarizationLevel,
    " rows; target FDR ", format(x$settings$targetER),
    "; calculated FDR ", format(x$fdr), ".\n",
    sep = ""
  )
  print(x$classificationSummary)
  invisible(x)
}
