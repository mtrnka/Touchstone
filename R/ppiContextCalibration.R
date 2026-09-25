#' Calibrate and classify context-aware protein-pair results
#'
#' Estimates a local posterior error probability (`contextPEP`) separately for
#' each PPI evidence-context group. Target-decoy and double-decoy densities are
#' combined using Touchstone's decoy-scaling calculation, sparse group curves
#' are shrunk toward the pooled curve, and weighted isotonic regression enforces
#' the common-sense requirement that error probability cannot increase as the
#' primary classifier improves.
#'
#' The reported `contextQValue` is the cumulative mean context PEP through all
#' PPIs with equal or better context PEP. Complete tied plateaus are accepted or
#' rejected together. The ordinary PPI classification defined by
#' `coreThreshold` is retained; context classification can add candidates but
#' never remove core-supported PPIs. If their union exceeds `targetER` by
#' Touchstone's target-decoy FDR calculation, contextual additions are trimmed
#' from least to most confident.
#'
#' Bootstrap selection frequency is a stability annotation, not an additional
#' error estimate. It separates contextual additions into `stably-enhanced`
#' and `unstably-enhanced` reporting tiers.
#'
#' @param context A `touchstone_ppi_context` object returned by
#'   [annotatePPIContext()].
#' @param targetER Requested protein-pair error rate.
#' @param bootstrapReplicates Number of stratified bootstrap fits. Use `0` to
#'   omit stability estimation.
#' @param seed Random seed used only for bootstrap resampling.
#' @param stabilityThreshold Minimum bootstrap selection frequency for the
#'   `stably-enhanced` tier.
#' @param bandwidthAdjust Multiplier applied to the score-density bandwidth.
#' @param priorDecoys Number of decoys controlling shrinkage of each context
#'   curve toward the pooled curve.
#' @return A `touchstone_ppi_results` object. `PPIs` contains one row per PPI
#'   candidate with context PEP, context q-value, bootstrap stability, and one
#'   nonredundant classification tier. `URPs` and `CSMs` retain keyed evidence;
#'   model, thresholds, FDR, and settings preserve classification provenance.
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
    dplyr::filter(.data$classified)
  classification.summary <- countDecoys(
    selected.ppis,
    scalingFactor = context$settings$scalingFactor
  )
  reported.ppis <- classification$PPIs %>%
    dplyr::select(-dplyr::any_of(c(
      "coreSupported", "contextQualified", "contextSelected", "classified"
    )))
  target.fdr.reached <- is.finite(classification$fdr) &&
    classification$fdr <= targetER

  structure(
    list(
      PPIs = reported.ppis,
      URPs = context$URPs,
      CSMs = context$CSMs,
      model = model,
      bootstrap = bootstrap,
      thresholds = list(
        coreThreshold = context$settings$coreThreshold,
        initialContextPEPThreshold = classification$initialThreshold,
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
          classification = "ordinary-or-context",
          contextTrimmed = classification$contextTrimmed,
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
#' @param tiers Optional classification tiers to retain.
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
    result <- dplyr::filter(result, !is.na(.data$classificationTier))
  }
  if (view == "clean") {
    result <- dplyr::filter(result, .data$Decoy == "Target")
  }
  if (!is.null(tiers)) {
    allowed <- levels(x$PPIs$classificationTier)
    unknown <- setdiff(tiers, allowed)
    if (length(unknown) > 0) {
      stop(
        "Unknown classification tier(s): ",
        paste(unknown, collapse = ", "), ".",
        call. = FALSE
      )
    }
    result <- dplyr::filter(
      result, as.character(.data$classificationTier) %in% tiers
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
  cat("Classification tiers:\n")
  print(table(classified$classificationTier, useNA = "no"))
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

contextPEPQValues <- function(contextPEP, decoyClass) {
  targetPEP <- sort(
    contextPEP[decoyClass == "Target" & is.finite(contextPEP)]
  )
  result <- rep(NA_real_, length(contextPEP))
  if (length(targetPEP) == 0) return(result)

  tied <- rle(targetPEP)
  groupEnds <- cumsum(tied$lengths)
  qValues <- cumsum(tied$values * tied$lengths) / groupEnds
  finite <- is.finite(contextPEP)
  positions <- findInterval(contextPEP[finite], tied$values)
  positions <- pmin(pmax(positions, 1L), length(tied$values))
  result[finite] <- qValues[positions]
  result
}

contextPEPThreshold <- function(contextPEP, decoyClass, targetER) {
  qValues <- contextPEPQValues(contextPEP, decoyClass)
  eligible <- decoyClass == "Target" & is.finite(contextPEP) &
    is.finite(qValues) & qValues <= targetER
  if (!any(eligible)) -Inf else max(contextPEP[eligible])
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
      pep[, iteration], data$Decoy, targetER
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
  initialThreshold <- contextPEPThreshold(
    data$contextPEP, data$Decoy, targetER
  )
  contextQValue <- contextPEPQValues(data$contextPEP, data$Decoy)
  selectionFrequency <- if (is.null(bootstrap)) {
    rep(NA_real_, nrow(data))
  } else {
    bootstrap$perPPI$selectionFrequency
  }
  contextQualified <- is.finite(contextQValue) & contextQValue <= targetER
  selectedFDR <- function(selected) {
    if (!any(data$Decoy[selected] == "Target")) return(Inf)
    as.numeric(calculateFDR(
      data[selected, , drop = FALSE],
      threshold = -Inf,
      classifier = scoreColumn,
      scalingFactor = scalingFactor
    ))
  }

  finalThreshold <- initialThreshold
  classified <- data$coreSupported | contextQualified
  combinedFDR <- selectedFDR(classified)
  if (!is.finite(combinedFDR) || combinedFDR > targetER) {
    acceptable <- FALSE
    cutoffs <- c(sort(unique(data$contextPEP[
      !data$coreSupported & contextQualified
    ]), decreasing = TRUE), -Inf)
    for (cutoff in cutoffs) {
      trial <- data$coreSupported |
        (!data$coreSupported & is.finite(data$contextPEP) &
           data$contextPEP <= cutoff)
      trialFDR <- selectedFDR(trial)
      if (is.finite(trialFDR) && trialFDR <= targetER) {
        finalThreshold <- cutoff
        classified <- trial
        combinedFDR <- trialFDR
        acceptable <- TRUE
        break
      }
    }
    if (!acceptable) {
      finalThreshold <- -Inf
      classified <- data$coreSupported
      combinedFDR <- selectedFDR(classified)
    }
  }

  contextSelected <- !data$coreSupported & classified
  classificationTier <- dplyr::case_when(
    data$coreSupported ~ "core-supported",
    contextSelected & !is.na(selectionFrequency) &
      selectionFrequency >= stabilityThreshold ~ "stably-enhanced",
    contextSelected ~ "unstably-enhanced",
    TRUE ~ NA_character_
  )
  annotated <- dplyr::mutate(
    data,
    contextQValue = contextQValue,
    contextQualified = contextQualified,
    contextSelected = contextSelected,
    selectionFrequency = selectionFrequency,
    classificationTier = factor(
      classificationTier,
      levels = c(
        "core-supported", "stably-enhanced", "unstably-enhanced"
      )
    ),
    classified = classified
  )
  list(
    PPIs = annotated,
    threshold = finalThreshold,
    initialThreshold = initialThreshold,
    fdr = combinedFDR,
    contextTrimmed = finalThreshold < initialThreshold
  )
}
