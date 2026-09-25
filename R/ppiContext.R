#' Annotate protein-pair candidates with fused contextual evidence
#'
#' Constructs target-decoy-symmetric protein identities from the reported
#' species and protein names, summarizes unique residue-pair evidence, and adds
#' protein-presence and conservative-network annotations. Target/decoy origin
#' is retained in `Decoy`; fusion is used only to calculate context.
#'
#' Context annotation accepts Touchstone's decoy scaling factor, but contextual
#' PPI error calibration with scaled decoys remains experimental. Network
#' features are calculated leave-one-edge-out from protein pairs passing
#' `coreThreshold`, so a candidate interaction cannot create its own support.
#'
#' @param datTab Scored CSM-level cross-link data.
#' @param coreThreshold Numeric classifier threshold defining the conservative
#'   protein-pair network.
#' @param candidateThreshold Minimum classifier value for candidate URPs and
#'   protein pairs. Use `-Inf` to use all scored candidates.
#' @param supportThreshold Minimum classifier value for URPs allowed to provide
#'   distinct-residue or intra-protein contextual support. Defaults to
#'   `candidateThreshold` for backward compatibility.
#' @param classifier Score column used for summarization and thresholds.
#' @param scalingFactor Decoy database scaling factor. A value other than 1
#'   produces an experimental-calibration warning.
#' @param retainGroups Retain existing groups during URP and PPI summarization.
#' @return A `touchstone_ppi_context` object containing annotated `PPIs`,
#'   annotated `URPs`, and analysis settings.
#' @export
annotatePPIContext <- function(datTab,
                               coreThreshold,
                               candidateThreshold = -Inf,
                               supportThreshold = candidateThreshold,
                               classifier = "SVM.score",
                               scalingFactor = the$decoyScalingFactor,
                               retainGroups = FALSE) {
  classifier <- .classifierName(rlang::enquo(classifier))
  if (length(scalingFactor) != 1 || !is.finite(scalingFactor) ||
      scalingFactor <= 0) {
    stop("scalingFactor must be one positive, finite value.", call. = FALSE)
  }
  if (!isTRUE(all.equal(scalingFactor, 1))) {
    warning(
      "PPI context calibration with scaled decoys is experimental; ",
      "1x decoys are currently recommended.",
      call. = FALSE
    )
  }
  if (length(coreThreshold) != 1 || !is.finite(coreThreshold)) {
    stop("coreThreshold must be one finite numeric value.", call. = FALSE)
  }
  if (length(candidateThreshold) != 1 || is.na(candidateThreshold)) {
    stop("candidateThreshold must be one numeric value.", call. = FALSE)
  }
  if (length(supportThreshold) != 1 || is.na(supportThreshold)) {
    stop("supportThreshold must be one numeric value.", call. = FALSE)
  }
  required <- c(
    classifier, "Acc.1", "Acc.2", "Protein.1", "Protein.2",
    "Species.1", "Species.2", "XLink.AA.1", "XLink.AA.2",
    "xlinkClass", "Decoy"
  )
  missing.columns <- setdiff(required, names(datTab))
  if (length(missing.columns) > 0) {
    stop(
      "PPI context annotation requires column(s): ",
      paste(missing.columns, collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (!all(c("xlinkedResPair", "xlinkedProtPair") %in% names(datTab))) {
    datTab <- calculatePairs(datTab, scalingFactor = scalingFactor)
  }

  fused <- addFusedContextIdentity(datTab)
  scored <- dplyr::bind_cols(datTab, fused)
  urps <- bestResPair(
    scored,
    classifier = !!rlang::sym(classifier),
    retainGroups = retainGroups
  )
  candidate.urps <- urps %>%
    dplyr::filter(.data[[classifier]] >= candidateThreshold)
  support.urps <- urps %>%
    dplyr::filter(.data[[classifier]] >= supportThreshold)

  intra.urps <- support.urps %>%
    dplyr::filter(.data$xlinkClass == "intraProtein") %>%
    dplyr::distinct(.data$fusedProtein1, .data$fusedResiduePair, .keep_all = TRUE)
  intra.support <- if (nrow(intra.urps) == 0) {
    tibble::tibble(
      fusedProtein1 = character(), intraURPs = integer(),
      intraMaxScore = double()
    )
  } else {
    intra.urps %>%
      dplyr::group_by(.data$fusedProtein1) %>%
      dplyr::summarize(
        intraURPs = dplyr::n(),
        intraMaxScore = max(.data[[classifier]]),
        .groups = "drop"
      )
  }

  inter.urps <- support.urps %>%
    dplyr::filter(.data$xlinkClass == "interProtein") %>%
    dplyr::arrange(.data$fusedProteinPair, dplyr::desc(.data[[classifier]])) %>%
    dplyr::group_by(.data$fusedProteinPair) %>%
    dplyr::mutate(
      urpRank = dplyr::row_number(),
      bestPositionA = dplyr::first(.data$fusedPositionA),
      bestPositionB = dplyr::first(.data$fusedPositionB),
      fullyDistinctFromBest = .data$urpRank > 1L &
        .data$fusedPositionA != .data$bestPositionA &
        .data$fusedPositionB != .data$bestPositionB
    ) %>%
    dplyr::ungroup()

  pair.evidence <- inter.urps %>%
    dplyr::group_by(.data$fusedProteinPair) %>%
    dplyr::summarize(
      fusedProteinA = dplyr::first(.data$fusedProteinA),
      fusedProteinB = dplyr::first(.data$fusedProteinB),
      distinctURPs = dplyr::n_distinct(.data$fusedResiduePair),
      fullyDistinctURPs = sum(.data$fullyDistinctFromBest, na.rm = TRUE),
      distinctURPContext = hasFullyDistinctContext(
        .data$fusedPositionA,
        .data$fusedPositionB
      ),
      .groups = "drop"
    ) %>%
    dplyr::left_join(
      intra.support,
      by = c("fusedProteinA" = "fusedProtein1")
    ) %>%
    dplyr::rename(intraURPsA = "intraURPs",
                  intraMaxScoreA = "intraMaxScore") %>%
    dplyr::left_join(
      intra.support,
      by = c("fusedProteinB" = "fusedProtein1")
    ) %>%
    dplyr::rename(intraURPsB = "intraURPs",
                  intraMaxScoreB = "intraMaxScore") %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(c("intraURPsA", "intraURPsB")),
        ~ dplyr::coalesce(.x, 0L)
      ),
      bothIntraSupported = .data$intraURPsA > 0 & .data$intraURPsB > 0
    )

  ppis <- bestProtPair(
    scored,
    classifier = !!rlang::sym(classifier),
    retainGroups = retainGroups
  ) %>%
    dplyr::filter(
      .data$xlinkClass == "interProtein",
      .data[[classifier]] >= candidateThreshold
    ) %>%
    dplyr::left_join(
      dplyr::select(pair.evidence, -"fusedProteinA", -"fusedProteinB"),
      by = "fusedProteinPair"
    ) %>%
    dplyr::mutate(
      distinctURPs = dplyr::coalesce(.data$distinctURPs, 0L),
      fullyDistinctURPs = dplyr::coalesce(.data$fullyDistinctURPs, 0L),
      distinctURPContext = dplyr::coalesce(.data$distinctURPContext, FALSE),
      intraURPsA = dplyr::coalesce(.data$intraURPsA, 0L),
      intraURPsB = dplyr::coalesce(.data$intraURPsB, 0L),
      bothIntraSupported = dplyr::coalesce(.data$bothIntraSupported, FALSE)
    )

  core.edges <- ppis %>%
    dplyr::filter(.data[[classifier]] >= coreThreshold) %>%
    dplyr::distinct(.data$fusedProteinPair, .keep_all = TRUE) %>%
    dplyr::select("fusedProteinA", "fusedProteinB")
  network <- calculateLeaveOneEdgeOutContext(
    ppis$fusedProteinA,
    ppis$fusedProteinB,
    core.edges
  )
  ppis <- dplyr::bind_cols(ppis, network) %>%
    dplyr::mutate(
      coreSupported = .data[[classifier]] >= coreThreshold,
      bothCoreConnected = .data$coreDegreeA > 0 & .data$coreDegreeB > 0,
      networkEmbedded = .data$commonCoreNeighbors > 0,
      contextGroup = dplyr::case_when(
        .data$distinctURPContext ~ "distinct",
        .data$bothCoreConnected ~ "core-connected",
        .data$bothIntraSupported ~ "intra-supported",
        TRUE ~ "context-poor"
      )
    )

  structure(
    list(
      PPIs = ppis,
      URPs = candidate.urps,
      settings = list(
        classifier = classifier,
        coreThreshold = coreThreshold,
        candidateThreshold = candidateThreshold,
        supportThreshold = supportThreshold,
        scalingFactor = scalingFactor,
        identity = "species-and-protein-name",
        scaledContextCalibration = if (scalingFactor == 1) {
          "validated-prototype"
        } else {
          "experimental"
        }
      )
    ),
    class = "touchstone_ppi_context"
  )
}

addFusedContextIdentity <- function(datTab) {
  protein1 <- stringr::str_c(datTab$Species.1, "::", datTab$Protein.1)
  protein2 <- stringr::str_c(datTab$Species.2, "::", datTab$Protein.2)
  swap <- protein1 > protein2
  protein.a <- ifelse(swap, protein2, protein1)
  protein.b <- ifelse(swap, protein1, protein2)
  position.a <- ifelse(swap, datTab$XLink.AA.2, datTab$XLink.AA.1)
  position.b <- ifelse(swap, datTab$XLink.AA.1, datTab$XLink.AA.2)
  tibble::tibble(
    fusedProtein1 = protein1,
    fusedProtein2 = protein2,
    fusedProteinA = protein.a,
    fusedProteinB = protein.b,
    fusedPositionA = position.a,
    fusedPositionB = position.b,
    fusedProteinPair = stringr::str_c(protein.a, protein.b, sep = "::"),
    fusedResiduePair = stringr::str_c(
      protein.a, "@", position.a, "::", protein.b, "@", position.b
    )
  )
}

hasFullyDistinctContext <- function(position.a, position.b) {
  if (length(position.a) < 2) return(FALSE)
  any(outer(position.a, position.a, "!=") &
        outer(position.b, position.b, "!="))
}

calculateLeaveOneEdgeOutContext <- function(protein.a, protein.b, core.edges) {
  adjacency <- list()
  if (nrow(core.edges) > 0) {
    for (index in seq_len(nrow(core.edges))) {
      a <- core.edges$fusedProteinA[[index]]
      b <- core.edges$fusedProteinB[[index]]
      adjacency[[a]] <- unique(c(adjacency[[a]], b))
      adjacency[[b]] <- unique(c(adjacency[[b]], a))
    }
  }
  purrr::map2_dfr(protein.a, protein.b, function(a, b) {
    neighbors.a <- setdiff(adjacency[[a]] %||% character(), b)
    neighbors.b <- setdiff(adjacency[[b]] %||% character(), a)
    tibble::tibble(
      coreDegreeA = length(neighbors.a),
      coreDegreeB = length(neighbors.b),
      commonCoreNeighbors = length(intersect(neighbors.a, neighbors.b))
    )
  })
}
