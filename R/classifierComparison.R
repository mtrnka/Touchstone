#' Compare classifier recovery across the FDR range
#'
#' Builds on [generateErrorTable.sep()] to evaluate several score columns on the
#' same CSM evidence. Each classifier independently chooses the representative
#' match at the requested proteomics summarization level, after which target
#' recovery and decoy-estimated FDR are calculated over its score range. This
#' does not fit or retrain a model.
#'
#' The returned long table can be inspected directly or passed to
#' [plotClassifierComparison()]. Its `threshold` column remains on each
#' classifier's native scale; `fdr` and `hits` provide the scale-independent
#' comparison.
#'
#' @param x A CSM data frame, a `touchstone_training` result, or a
#'   `touchstone_results` object. Results objects use their complete retained
#'   CSM source when available.
#' @param classifiers Character vector naming numeric score columns to compare.
#' @param summarizationLevel Proteomics level at which recovery is evaluated.
#'   One of `"csm"`, `"urp"`, `"peptide-pair"`, `"protein-pair"`, or
#'   `"module-pair"`.
#' @param model Model selected when `x` is a training result; passed to the same
#'   resolver used by [prepareCrosslinkResults()].
#' @param scalingFactor Decoy-database scaling factor. By default this is read
#'   from `x` when available, otherwise from [setDecoyScalingFactor()].
#' @param retainGroups Retain existing groups during summarization.
#' @param polishing Optional named polishing rules accepted by
#'   [classifyCrosslinkResults()]. These filters are applied once, before the
#'   classifier-specific summarization, without applying an FDR threshold.
#' @return A tibble with classifier, summarization level, crosslink class,
#'   native score threshold, calculated FDR, and recovered target hits. The
#'   polishing audit is retained as a `polishingAudit` attribute.
#' @examples
#' \dontrun{
#' comparison <- compareClassifiers(
#'   scored.csms,
#'   classifiers = c("SVM.score", "Score.Diff", "experimental.score"),
#'   summarizationLevel = "urp",
#'   polishing = list(minIons = 3)
#' )
#' plotClassifierComparison(comparison, targetER = 0.01)
#' }
#' @export
compareClassifiers <- function(x,
                               classifiers,
                               summarizationLevel = "urp",
                               model = "selected",
                               scalingFactor = NULL,
                               retainGroups = TRUE,
                               polishing = NULL) {
  if (!is.character(classifiers) || length(classifiers) < 1 ||
      anyNA(classifiers) || any(!nzchar(classifiers))) {
    stop("classifiers must be a non-empty character vector.", call. = FALSE)
  }
  classifiers <- unique(classifiers)
  summarizationLevel <- normalizeSummarizationLevel(summarizationLevel)

  stored.scaling <- NULL
  if (inherits(x, "touchstone_results")) {
    stored.scaling <- x$settings$scalingFactor
    csms <- if (is.data.frame(x$sourceCSMs)) x$sourceCSMs else x$data
  } else {
    resolved <- resolveCrosslinkFit(x, model)
    stored.scaling <- resolved$settings$scalingFactor
    csms <- resolved$fit$CSMs
    if (is.null(csms)) csms <- resolved$fit$scoredCSMs
  }
  if (!is.data.frame(csms)) {
    stop("x does not contain a usable CSM data frame.", call. = FALSE)
  }
  missing.required <- setdiff(c("xlinkClass", "Decoy"), names(csms))
  if (length(missing.required) > 0) {
    stop(
      "CSM data are missing required column(s): ",
      paste(missing.required, collapse = ", "),
      call. = FALSE
    )
  }

  missing.classifiers <- setdiff(classifiers, names(csms))
  if (length(missing.classifiers) > 0) {
    stop(
      "CSM data are missing classifier column(s): ",
      paste(missing.classifiers, collapse = ", "),
      call. = FALSE
    )
  }
  nonnumeric <- classifiers[!vapply(csms[classifiers], is.numeric, logical(1))]
  if (length(nonnumeric) > 0) {
    stop(
      "Classifier columns must be numeric: ",
      paste(nonnumeric, collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(scalingFactor)) scalingFactor <- stored.scaling
  if (is.null(scalingFactor)) scalingFactor <- the$decoyScalingFactor
  if (length(scalingFactor) != 1 || !is.finite(scalingFactor) ||
      scalingFactor <= 0) {
    stop("scalingFactor must be one positive, finite number.", call. = FALSE)
  }

  polishing.result <- applyCrosslinkPolishing(
    csms,
    normalizePolishingOptions(polishing)
  )
  csms <- polishing.result$data
  if (nrow(csms) == 0) {
    stop("No CSMs remain after polishing.", call. = FALSE)
  }

  comparison <- purrr::map_dfr(classifiers, function(classifier) {
    score <- csms[[classifier]]
    if (!any(is.finite(score))) {
      stop(
        "Classifier '", classifier, "' has no finite values after polishing.",
        call. = FALSE
      )
    }
    classifier.csms <- csms[is.finite(score), , drop = FALSE]
    finite.range <- range(classifier.csms[[classifier]])
    if (diff(finite.range) == 0) {
      stop(
        "Classifier '", classifier,
        "' is constant after polishing and cannot define a threshold curve.",
        call. = FALSE
      )
    }
    summarized <- summarizeCrosslinkData(
      classifier.csms,
      summarizationLevel = summarizationLevel,
      classifier = classifier,
      retainGroups = retainGroups
    )
    error.table <- generateErrorTable.sep(
      summarized,
      classifier = classifier,
      scalingFactor = scalingFactor
    )
    dplyr::bind_rows(
      error.table %>%
        dplyr::transmute(
          classifier = classifier,
          summarizationLevel = summarizationLevel,
          xlinkClass = "interProtein",
          threshold = .data$thresh,
          fdr = .data$fdr.inter,
          hits = .data$inter
        ),
      error.table %>%
        dplyr::transmute(
          classifier = classifier,
          summarizationLevel = summarizationLevel,
          xlinkClass = "intraProtein",
          threshold = .data$thresh,
          fdr = .data$fdr.intra,
          hits = .data$intra
        )
    )
  })
  attr(comparison, "polishingAudit") <- polishing.result$audit
  attr(comparison, "scalingFactor") <- scalingFactor
  comparison
}

#' Plot a comparison of classifier recovery across the FDR range
#'
#' Displays target hits recovered by each classifier as the allowed FDR
#' increases. The emphasized step curve is the best-attainable envelope at each
#' FDR allowance; the raw threshold path can be retained in the background to
#' reveal jagged or non-monotonic behavior.
#'
#' @param comparison Table returned by [compareClassifiers()].
#' @param targetER Desired FDR shown with a vertical reference line.
#' @param maxFDR Largest FDR value displayed.
#' @param linkClass Crosslink classes to display: `"both"`, `"inter"`, or
#'   `"intra"`.
#' @param showRaw Show each classifier's raw threshold path faintly behind its
#'   best-attainable envelope.
#' @return A `ggplot2` object.
#' @export
plotClassifierComparison <- function(comparison,
                                     targetER = 0.01,
                                     maxFDR = 0.05,
                                     linkClass = c("both", "inter", "intra"),
                                     showRaw = TRUE) {
  linkClass <- match.arg(linkClass)
  required <- c("classifier", "xlinkClass", "threshold", "fdr", "hits")
  missing.columns <- setdiff(required, names(comparison))
  if (length(missing.columns) > 0) {
    stop(
      "comparison is missing required column(s): ",
      paste(missing.columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (length(targetER) != 1 || !is.finite(targetER) || targetER < 0) {
    stop("targetER must be one non-negative, finite number.", call. = FALSE)
  }
  if (length(maxFDR) != 1 || !is.finite(maxFDR) || maxFDR <= 0) {
    stop("maxFDR must be one positive, finite number.", call. = FALSE)
  }

  selected.classes <- switch(
    linkClass,
    both = c("interProtein", "intraProtein"),
    inter = "interProtein",
    intra = "intraProtein"
  )
  plot.data <- comparison %>%
    dplyr::filter(
      .data$xlinkClass %in% selected.classes,
      is.finite(.data$fdr), is.finite(.data$hits),
      .data$fdr >= 0, .data$fdr <= maxFDR
    )
  if (nrow(plot.data) == 0) {
    stop("No finite comparison values fall within maxFDR.", call. = FALSE)
  }

  frontier <- plot.data %>%
    dplyr::group_by(.data$classifier, .data$xlinkClass, .data$fdr) %>%
    dplyr::summarize(hits = max(.data$hits), .groups = "drop") %>%
    dplyr::arrange(.data$classifier, .data$xlinkClass, .data$fdr) %>%
    dplyr::group_by(.data$classifier, .data$xlinkClass) %>%
    dplyr::mutate(hits = cummax(.data$hits)) %>%
    dplyr::ungroup()

  result <- ggplot2::ggplot(
    frontier,
    ggplot2::aes(
      x = .data$fdr,
      y = .data$hits,
      color = .data$classifier,
      group = .data$classifier
    )
  )
  if (isTRUE(showRaw)) {
    result <- result + ggplot2::geom_path(
      data = plot.data,
      ggplot2::aes(
        x = .data$fdr,
        y = .data$hits,
        color = .data$classifier,
        group = .data$classifier
      ),
      linewidth = 0.45,
      alpha = 0.25
    )
  }
  result +
    ggplot2::geom_step(linewidth = 1.1, direction = "hv") +
    ggplot2::geom_vline(xintercept = targetER, color = "red") +
    ggplot2::facet_grid(
      rows = ggplot2::vars(.data$xlinkClass),
      scales = "free_y"
    ) +
    ggplot2::labs(
      x = "Estimated FDR",
      y = "Target hits",
      color = "Classifier"
    ) +
    ggplot2::theme_bw()
}
