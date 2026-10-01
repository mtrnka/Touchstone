#' Calibrate and classify context-aware protein-pair results
#'
#' Estimates a local posterior error probability (`contextPEP`) separately for
#' each PPI evidence-context group. Target-decoy and double-decoy densities are
#' combined using Touchstone's decoy-scaling calculation, sparse group curves
#' are shrunk toward the pooled curve, and weighted isotonic regression enforces
#' the common-sense requirement that error probability cannot increase as the
#' primary classifier improves.
#'
#' The fitted `contextPEP` ranks PPIs, while direct cumulative target-decoy
#' q-values determine the classification boundary. Complete tied context-PEP
#' plateaus are accepted or rejected together, using Touchstone's target-decoy,
#' double-decoy, and scaling-factor arithmetic. Ordinary PPI calls are not
#' protected from re-evaluation, but remain annotated in the complete table.
#'
#' Bootstrap selection frequency is a stability annotation, not an additional
#' error estimate. It distinguishes stable and unstable context-enhanced calls
#' in the reported classification status.
#'
#' @param context A `touchstone_ppi_context` object returned by
#'   [annotatePPIContext()].
#' @param targetER Requested protein-pair error rate.
#' @param bootstrapReplicates Number of stratified bootstrap fits. Use `0` to
#'   omit stability estimation.
#' @param seed Random seed used only for bootstrap resampling.
#' @param stabilityThreshold Minimum bootstrap selection frequency for a
#'   `Stably context-enhanced` classification status.
#' @param bandwidthAdjust Multiplier applied to the score-density bandwidth.
#' @param priorDecoys Number of decoys controlling shrinkage of each context
#'   curve toward the pooled curve.
#' @return A `touchstone_ppi_results` object. `PPIs` contains one row per PPI
#'   candidate with context PEP, empirical FDR and q-value, bootstrap stability,
#'   ordinary-core status, context-selection status, and one classification
#'   label. `URPs` and `CSMs` retain keyed evidence; model, thresholds, FDR, and
#'   settings preserve classification provenance.
#' @export
classifyPPIContext <- function(context,
                               targetER = 0.02,
                               bootstrapReplicates = 200,
                               seed = 1,
                               stabilityThreshold = 0.90,
                               bandwidthAdjust = 1,
                               priorDecoys = 20) {
  if (!inherits(context, "touchstone_ppi_context")) {
    stop(
      "context must be a touchstone_ppi_context object from ",
      "annotatePPIContext().",
      call. = FALSE
    )
  }
  if (length(targetER) != 1 || !is.finite(targetER) ||
      targetER <= 0 || targetER >= 1) {
    stop("targetER must be one number between 0 and 1.", call. = FALSE)
  }
  if (length(bootstrapReplicates) != 1 ||
      !is.finite(bootstrapReplicates) || bootstrapReplicates < 0 ||
      bootstrapReplicates != as.integer(bootstrapReplicates)) {
    stop("bootstrapReplicates must be one non-negative integer.", call. = FALSE)
  }
  if (length(seed) != 1 || !is.finite(seed) || seed < 0 ||
      seed > .Machine$integer.max || seed != as.integer(seed)) {
    stop("seed must be one non-negative integer.", call. = FALSE)
  }
  if (length(stabilityThreshold) != 1 ||
      !is.finite(stabilityThreshold) || stabilityThreshold < 0 ||
      stabilityThreshold > 1) {
    stop("stabilityThreshold must be between 0 and 1.", call. = FALSE)
  }
  if (length(bandwidthAdjust) != 1 || !is.finite(bandwidthAdjust) ||
      bandwidthAdjust <= 0) {
    stop("bandwidthAdjust must be positive and finite.", call. = FALSE)
  }
  if (length(priorDecoys) != 1 || !is.finite(priorDecoys) ||
      priorDecoys < 0) {
    stop("priorDecoys must be non-negative and finite.", call. = FALSE)
  }

  ppis <- context$PPIs
  required <- c(
    context$settings$classifier, "contextGroup", "coreSupported", "Decoy",
    "xlinkClass"
  )
  missing.columns <- setdiff(required, names(ppis))
  if (length(missing.columns) > 0) {
    stop(
      "PPI context data are missing column(s): ",
      paste(missing.columns, collapse = ", "), ".",
      call. = FALSE
    )
  }

  model <- fitPPIContextPEP(
    ppis,
    scoreColumn = context$settings$classifier,
    scalingFactor = context$settings$scalingFactor,
    bandwidthAdjust = bandwidthAdjust,
    priorDecoys = priorDecoys
  )
  ppis$contextPEP <- predictPPIContextPEP(model, ppis)
  bootstrap <- if (bootstrapReplicates > 0) {
    bootstrapPPIContextPEP(
      ppis,
      targetER = targetER,
      replicates = bootstrapReplicates,
      seed = seed,
      scoreColumn = context$settings$classifier,
      scalingFactor = context$settings$scalingFactor,
      bandwidthAdjust = bandwidthAdjust,
      priorDecoys = priorDecoys
    )
  } else {
    NULL
  }
  classification <- classifyPPIContextTable(
    ppis,
    bootstrap = bootstrap,
    targetER = targetER,
    stabilityThreshold = stabilityThreshold,
    scoreColumn = context$settings$classifier,
    scalingFactor = context$settings$scalingFactor
  )
  selected.ppis <- classification$PPIs %>%
    dplyr::filter(.data$contextSelected)
  classification.summary <- countDecoys(
    selected.ppis,
    scalingFactor = context$settings$scalingFactor
  )
  reported.ppis <- classification$PPIs
  target.fdr.reached <- is.finite(classification$fdr) &&
    classification$fdr <= targetER

  structure(
    list(
      PPIs = reported.ppis,
      URPs = context$URPs,
      CSMs = context$CSMs,
      model = model,
      bootstrap = bootstrap,
      targetDecoyCurve = classification$targetDecoyCurve,
      thresholds = list(
        coreThreshold = context$settings$coreThreshold,
        contextPEPThreshold = classification$threshold
      ),
      fdr = list(
        requested = targetER,
        estimated = classification$fdr,
        counts = classification.summary,
        targetFDRReached = target.fdr.reached
      ),
      settings = c(
        context$settings,
        list(
          targetER = targetER,
          bootstrapReplicates = as.integer(bootstrapReplicates),
          seed = seed,
          stabilityThreshold = stabilityThreshold,
          bandwidthAdjust = bandwidthAdjust,
          priorDecoys = priorDecoys,
          calibration = "weighted-isotonic-context-PEP",
          classification = "direct-target-decoy-q-value",
          protectedCore = FALSE,
          targetFDRReached = target.fdr.reached
        )
      )
    ),
    class = "touchstone_ppi_results"
  )
}

#' Retrieve protein pairs from contextual PPI results
#'
#' Returns a view of the single canonical PPI table without storing separate
#' classified and target-only copies in the results object.
#'
#' @param x A `touchstone_ppi_results` object.
#' @param view One of `"classified"`, `"clean"`, or `"all"`. The clean view
#'   contains classified target PPIs only.
#' @param tiers Optional classification-status labels to retain.
#' @return A data frame containing the requested PPI view.
#' @export
getPPIs <- function(x,
                    view = c("classified", "clean", "all"),
                    tiers = NULL) {
  if (!inherits(x, "touchstone_ppi_results")) {
    stop("x must be a touchstone_ppi_results object.", call. = FALSE)
  }
  view <- match.arg(view)
  result <- x$PPIs
  if (view != "all") {
    result <- dplyr::filter(result, .data$contextSelected)
  }
  if (view == "clean") {
    result <- dplyr::filter(result, .data$Decoy == "Target")
  }
  if (!is.null(tiers)) {
    allowed <- levels(x$PPIs$classificationStatus)
    unknown <- setdiff(tiers, allowed)
    if (length(unknown) > 0) {
      stop(
        "Unknown classification status label(s): ",
        paste(unknown, collapse = ", "), ".",
        call. = FALSE
      )
    }
    result <- dplyr::filter(
      result, as.character(.data$classificationStatus) %in% tiers
    )
  }
  result
}

#' Retrieve the evidence supporting one protein pair
#'
#' @param x A `touchstone_ppi_results` object.
#' @param ppiID Internal PPI identifier from `x$PPIs$ppiID`.
#' @param proteinPair Alternatively, one reported `xlinkedProtPair` value.
#' @return A list containing the one-row `PPI` summary and its keyed `URPs` and
#'   `CSMs` evidence tables.
#' @export
getPPIEvidence <- function(x, ppiID = NULL, proteinPair = NULL) {
  if (!inherits(x, "touchstone_ppi_results")) {
    stop("x must be a touchstone_ppi_results object.", call. = FALSE)
  }
  supplied <- c(!is.null(ppiID), !is.null(proteinPair))
  if (sum(supplied) != 1) {
    stop("Supply exactly one of ppiID or proteinPair.", call. = FALSE)
  }
  if (!is.null(ppiID)) {
    if (length(ppiID) != 1 || is.na(ppiID)) {
      stop("ppiID must be one non-missing value.", call. = FALSE)
    }
    ppi <- dplyr::filter(x$PPIs, .data$ppiID == !!ppiID)
  } else {
    if (length(proteinPair) != 1 || is.na(proteinPair)) {
      stop("proteinPair must be one non-missing value.", call. = FALSE)
    }
    ppi <- dplyr::filter(
      x$PPIs, as.character(.data$xlinkedProtPair) == !!proteinPair
    )
  }
  if (nrow(ppi) == 0) {
    stop("No matching PPI was found.", call. = FALSE)
  }
  if (nrow(ppi) > 1) {
    stop(
      "proteinPair matched more than one PPI; select one by ppiID.",
      call. = FALSE
    )
  }
  selected.id <- ppi$ppiID[[1]]
  urps <- if (is.null(x$URPs)) {
    tibble::tibble()
  } else {
    dplyr::filter(x$URPs, .data$ppiID == selected.id)
  }
  csms <- if (is.null(x$CSMs)) {
    tibble::tibble()
  } else {
    dplyr::filter(x$CSMs, .data$ppiID == selected.id)
  }
  list(PPI = ppi, URPs = urps, CSMs = csms)
}

#' Print contextual protein-pair results
#'
#' @param x A `touchstone_ppi_results` object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.touchstone_ppi_results <- function(x, ...) {
  classified <- getPPIs(x, view = "classified")
  clean <- getPPIs(x, view = "clean")
  cat(
    "Touchstone contextual PPI results: ", nrow(x$PPIs), " candidates; ",
    nrow(classified), " classified rows; ", nrow(clean),
    " classified target PPIs.\n",
    sep = ""
  )
  cat(
    "Requested FDR ", format(x$fdr$requested),
    "; calculated target-decoy FDR ", format(x$fdr$estimated), ".\n",
    sep = ""
  )
  cat("Classification status:\n")
  print(table(classified$classificationStatus, useNA = "no"))
  invisible(x)
}

fitPPIContextPEP <- function(data,
                             scoreColumn,
                             scalingFactor,
                             bandwidthAdjust,
                             priorDecoys,
                             gridLength = 1024) {
  score <- data[[scoreColumn]]
  finite <- is.finite(score)
  data <- data[finite, , drop = FALSE]
  score <- score[finite]
  if (length(score) < 2) {
    stop("At least two finite PPI classifier values are required.", call. = FALSE)
  }
  limits <- range(score)
  grid <- seq(limits[[1]], limits[[2]], length.out = gridLength)
  bandwidth <- stats::bw.nrd0(score) * bandwidthAdjust
  if (!is.finite(bandwidth) || bandwidth <= 0) bandwidth <- 0.1

  countDensity <- function(values) {
    if (length(values) == 0) return(rep(0, length(grid)))
    if (length(values) == 1) {
      return(stats::dnorm(grid, mean = values[[1]], sd = bandwidth))
    }
    length(values) * stats::density(
      values,
      bw = bandwidth,
      from = limits[[1]],
      to = limits[[2]],
      n = length(grid)
    )$y
  }
  calculateCurve <- function(d) {
    target <- countDensity(d[[scoreColumn]][d$Decoy == "Target"])
    targetDecoy <- countDensity(d[[scoreColumn]][d$Decoy == "Decoy"])
    doubleDecoy <- countDensity(
      d[[scoreColumn]][d$Decoy == "DoubleDecoy"]
    )
    scaled <- .scaleDecoyEvidence(
      targetDecoy, doubleDecoy, scalingFactor
    )
    false <- scaled$ftTT + scaled$ffTT
    list(
      target = target,
      false = false,
      pep = pmin(pmax(false / pmax(target, 1e-12), 0), 1),
      decoys = sum(d$Decoy != "Target")
    )
  }

  pooled <- calculateCurve(data)
  groups <- split(data, as.character(data$contextGroup), drop = TRUE)
  curves <- lapply(groups, function(d) {
    curve <- calculateCurve(d)
    shrinkageWeight <- curve$decoys / (curve$decoys + priorDecoys)
    curve$rawPEP <- curve$pep
    curve$shrunkPEP <- shrinkageWeight * curve$pep +
      (1 - shrinkageWeight) * pooled$pep
    curve$pep <- weightedDecreasingPAVA(
      curve$shrunkPEP,
      curve$target + curve$false
    )
    curve$shrinkageWeight <- shrinkageWeight
    curve
  })
  structure(
    list(
      grid = grid,
      curves = curves,
      pooled = pooled,
      settings = list(
        scoreColumn = scoreColumn,
        groupColumn = "contextGroup",
        scalingFactor = scalingFactor,
        bandwidthAdjust = bandwidthAdjust,
        priorDecoys = priorDecoys,
        monotonicMethod = "weighted-isotonic"
      )
    ),
    class = "touchstone_ppi_context_pep_fit"
  )
}

predictPPIContextPEP <- function(model, newdata) {
  score <- newdata[[model$settings$scoreColumn]]
  group <- as.character(newdata[[model$settings$groupColumn]])
  prediction <- rep(NA_real_, length(score))
  for (groupName in unique(group)) {
    index <- which(group == groupName & is.finite(score))
    curve <- model$curves[[groupName]]
    values <- if (is.null(curve)) model$pooled$pep else curve$pep
    prediction[index] <- stats::approx(
      model$grid,
      values,
      xout = score[index],
      rule = 2,
      ties = "ordered"
    )$y
  }
  prediction
}

weightedDecreasingPAVA <- function(values, weights) {
  if (length(values) != length(weights)) {
    stop("values and weights must have equal lengths.", call. = FALSE)
  }
  if (length(values) < 2) return(values)
  weights[!is.finite(weights) | weights <= 0] <- 0
  weightFloor <- max(weights, na.rm = TRUE) * 1e-8
  if (!is.finite(weightFloor) || weightFloor <= 0) weightFloor <- 1
  weights <- pmax(weights, weightFloor)

  starts <- ends <- integer(length(values))
  blockValues <- blockWeights <- numeric(length(values))
  blocks <- 0L
  for (index in seq_along(values)) {
    blocks <- blocks + 1L
    starts[[blocks]] <- index
    ends[[blocks]] <- index
    blockValues[[blocks]] <- values[[index]]
    blockWeights[[blocks]] <- weights[[index]]
    while (blocks > 1L &&
           blockValues[[blocks - 1L]] < blockValues[[blocks]]) {
      combinedWeight <- blockWeights[[blocks - 1L]] +
        blockWeights[[blocks]]
      blockValues[[blocks - 1L]] <-
        (blockValues[[blocks - 1L]] * blockWeights[[blocks - 1L]] +
           blockValues[[blocks]] * blockWeights[[blocks]]) /
        combinedWeight
      blockWeights[[blocks - 1L]] <- combinedWeight
      ends[[blocks - 1L]] <- ends[[blocks]]
      blocks <- blocks - 1L
    }
  }

  fitted <- numeric(length(values))
  for (block in seq_len(blocks)) {
    fitted[starts[[block]]:ends[[block]]] <- blockValues[[block]]
  }
  pmin(pmax(fitted, 0), 1)
}

contextTargetDecoyCurve <- function(contextPEP,
                                    decoyClass,
                                    scalingFactor = 1) {
  if (length(contextPEP) != length(decoyClass)) {
    stop("contextPEP and decoyClass must have equal lengths.", call. = FALSE)
  }
  curve <- tibble::tibble(
    contextPEP = contextPEP,
    Decoy = as.character(decoyClass)
  ) %>%
    dplyr::filter(is.finite(.data$contextPEP)) %>%
    dplyr::group_by(.data$contextPEP) %>%
    dplyr::summarize(
      targets = sum(.data$Decoy == "Target"),
      targetDecoys = sum(.data$Decoy == "Decoy"),
      doubleDecoys = sum(.data$Decoy == "DoubleDecoy"),
      .groups = "drop"
    ) %>%
    dplyr::arrange(.data$contextPEP) %>%
    dplyr::mutate(
      targets = cumsum(.data$targets),
      targetDecoys = cumsum(.data$targetDecoys),
      doubleDecoys = cumsum(.data$doubleDecoys)
    )
  if (nrow(curve) == 0) {
    curve$contextFDR <- curve$contextQValue <- numeric()
    return(curve)
  }
  scaled <- .scaleDecoyEvidence(
    curve$targetDecoys, curve$doubleDecoys, scalingFactor
  )
  false.targets <- scaled$ftTT + scaled$ffTT
  curve$contextFDR <- ifelse(
    curve$targets > 0, false.targets / curve$targets, Inf
  )
  curve$contextQValue <- rev(cummin(rev(curve$contextFDR)))
  curve
}

contextTargetDecoyValues <- function(contextPEP,
                                     decoyClass,
                                     scalingFactor = 1) {
  curve <- contextTargetDecoyCurve(
    contextPEP, decoyClass, scalingFactor = scalingFactor
  )
  fdr <- q.value <- rep(NA_real_, length(contextPEP))
  finite <- is.finite(contextPEP)
  if (nrow(curve) > 0 && any(finite)) {
    index <- match(contextPEP[finite], curve$contextPEP)
    fdr[finite] <- curve$contextFDR[index]
    q.value[finite] <- curve$contextQValue[index]
  }
  list(fdr = fdr, qValue = q.value, curve = curve)
}

contextPEPQValues <- function(contextPEP,
                              decoyClass,
                              scalingFactor = 1) {
  contextTargetDecoyValues(
    contextPEP, decoyClass, scalingFactor = scalingFactor
  )$qValue
}

contextPEPThreshold <- function(contextPEP,
                                decoyClass,
                                targetER,
                                scalingFactor = 1) {
  curve <- contextTargetDecoyCurve(
    contextPEP, decoyClass, scalingFactor = scalingFactor
  )
  eligible <- is.finite(curve$contextQValue) &
    curve$contextQValue <= targetER
  if (!any(eligible)) -Inf else max(curve$contextPEP[eligible])
}

bootstrapPPIContextPEP <- function(data,
                                   targetER,
                                   replicates,
                                   seed,
                                   scoreColumn,
                                   scalingFactor,
                                   bandwidthAdjust,
                                   priorDecoys) {
  had.seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had.seed) {
    old.seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had.seed) {
      assign(".Random.seed", old.seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  groupIndices <- split(seq_len(nrow(data)), as.character(data$contextGroup))
  pep <- matrix(NA_real_, nrow = nrow(data), ncol = replicates)
  selected <- matrix(FALSE, nrow = nrow(data), ncol = replicates)
  thresholds <- targetCounts <- rep(NA_real_, replicates)

  for (iteration in seq_len(replicates)) {
    bootstrapIndex <- unlist(lapply(
      groupIndices,
      function(index) sample(index, length(index), replace = TRUE)
    ), use.names = FALSE)
    model <- fitPPIContextPEP(
      data[bootstrapIndex, , drop = FALSE],
      scoreColumn = scoreColumn,
      scalingFactor = scalingFactor,
      bandwidthAdjust = bandwidthAdjust,
      priorDecoys = priorDecoys
    )
    pep[, iteration] <- predictPPIContextPEP(model, data)
    thresholds[[iteration]] <- contextPEPThreshold(
      pep[, iteration], data$Decoy, targetER,
      scalingFactor = scalingFactor
    )
    selected[, iteration] <- pep[, iteration] <= thresholds[[iteration]]
    targetCounts[[iteration]] <- sum(
      selected[, iteration] & data$Decoy == "Target",
      na.rm = TRUE
    )
  }

  list(
    perPPI = data.frame(
      rowIndex = seq_len(nrow(data)),
      medianContextPEP = apply(pep, 1, stats::median, na.rm = TRUE),
      upperContextPEP = apply(
        pep, 1, stats::quantile, probs = 0.9, na.rm = TRUE
      ),
      selectionFrequency = rowMeans(selected, na.rm = TRUE)
    ),
    thresholds = thresholds,
    targetCounts = targetCounts,
    settings = list(
      targetER = targetER,
      replicates = replicates,
      seed = seed
    )
  )
}

classifyPPIContextTable <- function(data,
                                    bootstrap,
                                    targetER,
                                    stabilityThreshold,
                                    scoreColumn,
                                    scalingFactor) {
  context.values <- contextTargetDecoyValues(
    data$contextPEP,
    data$Decoy,
    scalingFactor = scalingFactor
  )
  threshold <- contextPEPThreshold(
    data$contextPEP,
    data$Decoy,
    targetER,
    scalingFactor = scalingFactor
  )
  selectionFrequency <- if (is.null(bootstrap)) {
    rep(NA_real_, nrow(data))
  } else {
    bootstrap$perPPI$selectionFrequency
  }
  contextSelected <- is.finite(context.values$qValue) &
    context.values$qValue <= targetER
  selectedFDR <- function(selected) {
    if (!any(data$Decoy[selected] == "Target")) return(Inf)
    as.numeric(calculateFDR(
      data[selected, , drop = FALSE],
      threshold = -Inf,
      classifier = scoreColumn,
      scalingFactor = scalingFactor
    ))
  }
  estimatedFDR <- selectedFDR(contextSelected)
  classificationStatus <- dplyr::case_when(
    data$coreSupported & contextSelected ~ "Ordinary core retained",
    !data$coreSupported & contextSelected &
      !is.na(selectionFrequency) &
      selectionFrequency >= stabilityThreshold ~ "Stably context-enhanced",
    !data$coreSupported & contextSelected &
      !is.na(selectionFrequency) ~ "Unstably context-enhanced",
    !data$coreSupported & contextSelected ~
      "Context-enhanced (stability not estimated)",
    data$coreSupported ~ "Ordinary core not context-selected",
    TRUE ~ "Unclassified candidate"
  )
  annotated <- dplyr::mutate(
    data,
    contextFDR = context.values$fdr,
    contextQValue = context.values$qValue,
    contextSelected = contextSelected,
    selectionFrequency = selectionFrequency,
    classificationStatus = factor(
      classificationStatus,
      levels = c(
        "Ordinary core retained",
        "Stably context-enhanced",
        "Unstably context-enhanced",
        "Context-enhanced (stability not estimated)",
        "Ordinary core not context-selected",
        "Unclassified candidate"
      )
    )
  )
  list(
    PPIs = annotated,
    threshold = threshold,
    fdr = estimatedFDR,
    targetDecoyCurve = context.values$curve
  )
}
