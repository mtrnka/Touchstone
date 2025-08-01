#' Replace accession number in search results.
#'
#' @param datTab Parsed CLMS search results.
#' @param oldName Accession number to replace.
#' @param newName New accession number.
#' @return a data frame
#' @export
#'
renameProtein <- function(datTab, oldName, newName) {
  datTab$Acc.1 <- gsub(oldName, newName, datTab$Acc.1)
  datTab$Acc.2 <- gsub(oldName, newName, datTab$Acc.2)
  datTab <- calculatePairs(datTab)
  return(datTab)
}

#' Renumber protein in search results.
#'
#' This would typically be used when the sequence of a protein construct used in
#' an experiment differs from the UniProt reference, due to an epitope tag or deletion,
#'and one would like to harmonize the results reporting.
#'
#' @param datTab Parsed CLMS search results.
#' @param protein Accession number to renumber.
#' @param shift Difference between old sequence position and new sequence position
#' @return a data frame
#' @export
#'
renumberProtein <-function(datTab, protein, shift) {
  datTab[datTab$Acc.1==protein,"XLink.AA.1"] <-
    datTab[datTab$Acc.1==protein,"XLink.AA.1"] + shift
  datTab[datTab$Acc.2==protein,"XLink.AA.2"] <-
    datTab[datTab$Acc.2==protein,"XLink.AA.2"] + shift
  datTab <- calculatePairs(datTab)
  return(datTab)
}

#' Contact map style dot plot showing crosslinked residue pairs
#'
#' Function takes a datTab which is reduced to Residue Pairs and plot
#' Dot for each unique crosslink with size equal to num of instances. The
#' color can reflect the source of the crosslink to facilitate comparisons.
#' Basic skeleton of the function that will get built up later. Ultimately
#' would like this to work on the XLSearchOutput object, which probably needs
#' flag to show if it is residue pair reduced etc.
#'
#' @param datTab Parsed CLMS search results.
#' @param modTab A parsed Touchstone module file.
#' @param color Named or hex coded color for dots.
#' @param accOrder Character vector of Accession numbers in their desired plot order.
#' @param removeMods Modules to remove from plot.
#' @param displayEmpty Draw module boxes even when they contain no crosslinks.#'
#' @returns the ggplot object
#' @export
#'
pairPlot <- function(datTab, modTab=NULL, color="lightseagreen", accOrder=NULL,
                      removeMods=NA_character_, displayEmpty=T) {
  #    datTab <- removeModule(datTab, removeMods)
  datTabR <- datTab
  tempXLinkAA1 <- datTab$XLink.AA.1
  tempAcc1 <- datTab$Acc.1
  datTabR$XLink.AA.1 <- datTab$XLink.AA.2
  datTabR$Acc.1 <- datTab$Acc.2
  datTabR$XLink.AA.2 <- tempXLinkAA1
  datTabR$Acc.2 <- tempAcc1
  datTab <- rbind(datTab, datTabR)
  if (is.null(accOrder)) {
    accOrder = forcats::lvls_union(list(datTab$Acc.1, datTab$Acc.2))
  }
  datTab[c("Acc.1", "Acc.2")] <- datTab %>%
    select(.data$Acc.1, .data$Acc.2) %>%
    as.list() %>%
    forcats::fct_unify(levels = accOrder) %>%
    bind_cols
  gg <- ggplot(datTab, aes(.data$XLink.AA.1, .data$XLink.AA.2)) +
    geom_point(aes(size=.data$numCSM)) +
    facet_grid(.data$Acc.2 ~ .data$Acc.1, space="free", scales="free",
               as.table=F, switch="both") +
    theme(panel.background = element_rect(fill=rgb(0.98,0.98,0.98)),
          panel.spacing = grid::unit(0.5, units="line"),
          panel.border = element_rect(color = "grey", fill = NA, size = 1),
          panel.grid.major = element_line(color = "grey", size=.5),
          panel.grid.minor = element_line(color = "grey", size=.25),
          aspect.ratio=1) +#legend.position="none") +
    scale_x_continuous(breaks=seq(0,max(datTab$XLink.AA.1),by=50),
                       minor_breaks=seq(0,max(datTab$XLink.AA.1),by=25)) +
    scale_y_continuous(breaks=seq(0,max(datTab$XLink.AA.2),by=50),
                       minor_breaks=seq(0,max(datTab$XLink.AA.2),by=25)) +
    scale_size_area(limits=c(1,70)) +
    scale_alpha(limits=c(10,55), range=c(0.8,1)) +
    theme_bw()

  if (color=="byDistance") {
    gg <- gg + geom_point(aes(size=.data$numCSM, col=.data$distance <35))
  }   else if (color=="byQuant") {
    gg <- gg + geom_point(aes(size=.data$numCSM, col=.data$numCSM)) +
      scale_color_viridis_c()
  }   else {
    gg <- gg + geom_point(aes(size=.data$numCSM), col=color)
  }

  if (!is.null(modTab)) {
    modTab <- modTab %>%
      tidyr::replace_na(list(first.AA=0))
    modTabX <- modTab %>%
      select(.data$Acc, .data$first.AA) %>%
      rename(Acc.1 = .data$Acc, XLink.AA.1 = .data$first.AA) %>%
      mutate(Acc.1 = factor(.data$Acc.1, levels=accOrder))
    modTabY <- modTabX %>%
      rename(Acc.2 = .data$Acc.1, XLink.AA.2 = .data$XLink.AA.1) %>%
      mutate(Acc.2 = factor(.data$Acc.2, levels=accOrder))
    gg <- gg +
      geom_vline(data=modTabX, aes(xintercept = .data$XLink.AA.1), color = "palevioletred2", linetype="21", size=0.75, alpha=0.5) +
      geom_hline(data=modTabY, aes(yintercept = .data$XLink.AA.2), color = "palevioletred2", linetype="21", size=0.75, alpha=0.5)
  }

  plot(gg)
}

#' Removes all crosslinks from a given module from the dataset.
#'
#' @param datTab Parsed CLMS search results.
#' @param modules character vector of modules to remove.
#'
#' @returns A data frame
#' @export
#'
removeModule <- function(datTab, modules) {
  datTab <- datTab[!(datTab$Module.1 %in% modules | datTab$Module.2 %in% modules),]
  return(datTab)
}

#' Set the decoy scaling factor
#'
#' The decoy scaling factor is an integer k that describes how how much larger the
#' decoy database is to the target database.  Should be set globally once for each
#' analysis.  Defaults to 1.
#' @param dsf an integer, the decoy scaling factor k.
#'
#' @returns No reutrn value.
#' @export
#'
setDecoyScalingFactor <- function(dsf) {
  the$decoyScalingFactor <- dsf
  }

#' Skeleton function to parse heavy / light crosslink pairs.
#'
#' Currently only works for D12/H12 BS3
#' @param datTab Parsed CLMS search results.
#'
#' @returns A data frame
#' @export
#'
parseCrosslinker <- function(datTab) {
  datTab$Xlinker <- "BS3"
  datTab[grepl("\\:2H\\(12\\)",datTab$Peptide.1) |
           grepl("\\:2H\\(12\\)",datTab$Peptide.2),"Xlinker"] <- "BS3hvy"
  datTab$Xlinker <- as.factor(datTab$Xlinker)
  return(datTab)
}

#' Internal function to generate url links to MS-Product spectra. For MS2 data.
#'
#' @param path Base url to MS-Product spectra (typically a path on MS-Viewer for instance)
#' @param fraction Fraction identifier
#' @param z Precursor charge
#' @param peptide.1 Peptide 1 sequence
#' @param peptide.2 Peptide 2 sequence
#' @param spectrum spectrum identifier
#' @param instrumentType Prospector instrument type parameter
#' @param linkType Prospector link type parameter
#' @param outputType MS-Product output type (typically "text" or "HTML")
#' @keywords internal
#'
#' @returns A character containing a uri.
#' @export
#'
generateMSViewerLink <- function(path, fraction, z, peptide.1, peptide.2, spectrum,
                                 instrumentType = NA, linkType = NA, outputType = "HTML") {
  if(!stringr::str_detect(fraction, "\\.[[a-z]]+$")) {
    fraction <- paste0(fraction, ".mgf")
  }
  if (is.na(instrumentType)) {
    instrumentType <- switch(
      stringr::str_extract(fraction, "[[a-z]]+(?=\\.[[a-z]]+$)"),
      "ESI-Q-high-res",
      ethcd = "ESI-EThcD-high-res",
      etd = "ESI-ETD-high-res",
      hcd = "ESI-Q-high-res",
      cid = "ESI-ION-TRAP-low-res"
      #Add other instrument types
    )
  }
  if (is.na(linkType)) {
    linkType <- stringr::str_extract(peptide.1, "(?<=\\(\\+).+?(?=\\)([[A-Z]]|$|-))")
  } else if (linkType == "DSSOms3") {
    peptide.1 <- stringr::str_replace(peptide.1,
                                      "XL:A-[[A-Z]][[a-z]]+(\\(Unsaturated\\))?",
                                      stringr::fixed("+DSSO*"))
    peptide.2 <- stringr::str_replace(peptide.2,
                                      "XL:A-[[A-Z]][[a-z]]+(\\(Unsaturated\\))?",
                                      stringr::fixed("+DSSO*"))
    linkType <- "DSSO*"
  }
  templateVals.ms2["output_type"] <- outputType
  templateVals.ms2["data_filename"] <- file.path(path, fraction)
  templateVals.ms2["instrument_name"] <- instrumentType
  templateVals.ms2["scan_number"] <- spectrum
  templateVals.ms2["max_charge"] <- z
  templateVals.ms2["msms_precursor_charge"] <- z
  templateVals.ms2["sequence"] <- peptide.1
  templateVals.ms2["sequence2"] <- peptide.2
  templateVals.ms2["link_search_type"] <- linkType
  templateNames.ms2 <- names(templateVals.ms2)
  templateVals.ms2 <- urltools::url_encode(templateVals.ms2)
  names(templateVals.ms2) <- templateNames.ms2
  zipped <- stringr::str_c(names(templateVals.ms2), templateVals.ms2, sep="=", collapse="&")
  if (outputType == "HTML") {
    return(stringr::str_c('<a href=\"', zipped, '\" target=\"_blank\">Spectrum</a>'))
  } else {
    return(zipped)
  }
}

#' Internal function to generate url links to MS-Product spectra. For MS3 data.
#'
#' @param path Base url to MS-Product spectra (typically a path on MS-Viewer for instance)
#' @param fraction Fraction identifier
#' @param z Precursor charge
#' @param peptide Peptide sequence
#' @param spectrum spectrum identifier
#' @param instrumentType Prospector instrument type parameter
#' @param outputType MS-Product output type (typically "text" or "HTML")
#' @keywords internal
#'
#' @returns A character containing a uri.
#' @export
#'
generateMSViewerLink.ms3 <- function(path, fraction, z, peptide, spectrum,
                                     instrumentType = NA, outputType = "HTML") {
  if(!stringr::str_detect(fraction, "\\.[[a-z]]+$")) {
    fraction <- paste0(fraction, ".mgf")
  }
  if (is.na(instrumentType)) {
    instrumentType <- switch(
      stringr::str_extract(fraction, "[[a-z]]+(?=\\.[[a-z]]+$)"),
      "ESI-Q-high-res",
      ethcd = "ESI-EThcD-high-res",
      etd = "ESI-ETD-high-res",
      hcd = "ESI-Q-high-res",
      cid = "ESI-ION-TRAP-low-res"
      #Add other instrument types
    )
  }
  templateVals.ms3["output_type"] <- outputType
  templateVals.ms3["data_filename"] <- file.path(path, fraction)
  templateVals.ms3["instrument_name"] <- instrumentType
  templateVals.ms3["scan_number"] <- spectrum
  templateVals.ms3["max_charge"] <- z
  templateVals.ms3["msms_precursor_charge"] <- z
  templateVals.ms3["sequence"] <- peptide
  templateNames.ms3 <- names(templateVals.ms3)
  templateVals.ms3 <- urltools::url_encode(templateVals.ms3)
  names(templateVals.ms3) <- templateNames.ms3
  zipped <- stringr::str_c(names(templateVals.ms3), templateVals.ms3, sep="=", collapse="&")
  if (outputType == "HTML") {
    return(stringr::str_c('<a href=\"', zipped, '\" target=\"_blank\">Spectrum</a>'))
  } else {
    return(zipped)
  }
}

#' Formats the crosslink results for exporting / viewing
#'
#' @param datTab Parsed CLMS search results.
#' @param msviewer Formats column names to be compatible with MS-Viewer
#' @param extraCols Character vector of extra column names to be included in the report
#'
#' @returns A data frame
#' @export
#'
formatXLTable <- function(datTab, msviewer=F, extraCols=NULL) {
  annoyingColumns <- stringr::str_which(names(datTab), "(Int|Dec)[a-z]{2}\\.[[1-2]]")
  if (length(annoyingColumns) > 0) {
    datTab <- datTab[, -annoyingColumns]
  }
  datTab <- datTab %>%
    select(-starts_with("Res"),
           -starts_with("Num\\."),
           -any_of("massError"))
  if (sum(!is.na(datTab$distance)) == 0) {
    datTab <- datTab %>%
      select(-any_of("distance"))
  }
  if (sum(stringr::str_detect(names(datTab), "SVM.score")) > 0) datTab <- datTab[order(datTab$SVM.score, decreasing = T),]
  columnsToReport <-c(
    "keep", "specMS2", "specMS3.1", "specMS3.2", "Decoy", "groundTruth",
    "xlinkedResPair", "xlinkedProtPair", "xlinkedModulPair",
    "SVM.score", "SVM.new", "distance", "m/z", "z", "ppm",
    "DB.Peptide.1", "DB.Peptide.2", "Score", "Score.Diff", "percMatched",
    "Sc.1", "Rk.1", "Sc.2", "Rk.2", "Acc.1",
    "XLink.AA.1", "Protein.1", "Module.1", "Species.1",
    "Acc.2", "XLink.AA.2", "Protein.2", "Module.2", "Species.2",
    "xlinkClass", "Len.Pep.1", "Len.Pep.2",
    "Peptide.1", "Peptide.2", "numCSM", "numURP",
    "Fraction", "RT", "MSMS.Info", "Instrument", "id", "experiment", "Manual.Inspection"
  )
  if (!is.null(extraCols)) columnsToReport <- c(columnsToReport, extraCols)
  datTab <- datTab %>% select(any_of(columnsToReport))
  if ("Protein.1" %in% names(datTab) & "Protein.2" %in% names(datTab)) {
    datTab <- datTab %>%
      mutate(Protein.1 = factor(.data$Protein.1),
             Protein.2 = factor(.data$Protein.2))
  }
  if ("Module.1" %in% names(datTab) & "Module.2" %in% names(datTab)) {
    datTab <- datTab %>%
      mutate(Module.1 = factor(.data$Module.1),
             Module.2 = factor(.data$Module.2))
  }
  if ("Species.1" %in% names(datTab) & "Species.2" %in% names(datTab)) {
    datTab <- datTab %>%
      mutate(Species.1 = factor(.data$Species.1),
             Species.2 = factor(.data$Species.2))
  }
  if ("Acc.1" %in% names(datTab) & "Acc.2" %in% names(datTab)) {
    datTab <- datTab %>%
      mutate(Acc.1 = factor(.data$Acc.1),
             Acc.2 = factor(.data$Acc.2))
  }
  if ("xlinkClass" %in% names(datTab)) {
    datTab <- datTab %>%
      mutate(xlinkClass = factor(.data$xlinkClass))
  }
  if ("Fraction" %in% names(datTab)) {
    datTab$Fraction <- gsub("(.*)\\.[^.]+$", "\\1", datTab$Fraction)
    datTab <- datTab %>%
      mutate(Fraction = factor(.data$Fraction))
  }
  if ("percMatched" %in% names(datTab)) {
    datTab$percMatched <- round(datTab$percMatched, 2)
  }
  if ("Rk.1" %in% names(datTab) & "Rk.2" %in% names(datTab)) {
    datTab <- datTab %>%
      mutate(Rk.1 = as.integer(.data$Rk.1),
             Rk.2 = as.integer(.data$Rk.2))
  }
  if ("XLink.AA.1" %in% names(datTab) & "XLink.AA.2" %in% names(datTab)) {
    datTab <- datTab %>%
      mutate(XLink.AA.1 = as.integer(.data$XLink.AA.1),
             XLink.AA.2 = as.integer(.data$XLink.AA.2))
  }
  if ("Len.Pep.1" %in% names(datTab) & "Len.Pep.2" %in% names(datTab)) {
    datTab <- datTab %>%
      mutate(Len.Pep.1 = as.integer(.data$Len.Pep.1),
             Len.Pep.2 = as.integer(.data$Len.Pep.2))
  }
  if ("z" %in% names(datTab)) {
    datTab <- datTab %>%
      mutate(z = as.integer(.data$z))
  }
  #   if ("MSMS.Info" %in% names(datTab)) {
  #    datTab <- datTab %>%
  #   mutate(MSMS.Info = as.integer(.data$MSMS.Info))
  #   }
  if ("SVM.score" %in% names(datTab)) {
    datTab <- datTab %>%
      mutate(SVM.score = round(.data$SVM.score, 2))
  }
  if (msviewer) {datTab <- datTab %>% rename(
    `MSMS Info` = .data$MSMS.Info,
    `Peptide 1` = .data$Peptide.1,
    `Peptide 2` = .data$Peptide.2)
  }
  datTab <- datTab %>%
    mutate(xlinkedResPair = forcats::fct_drop(.data$xlinkedResPair))
  return(datTab)
}

#' Rescales data to plot and tabulate estimate decoy hits
#'
#' Called by the FDRplots function, but also useful when tabulating results.
#' For cases where a larger decoy database is searched to help estimate incorrect hits.
#' The decoy hits are divided by the scaling factor so that we tabulate the estimated
#' error rather than the artifically large number of decoys.
#' @param datTab Parsed CLMS search results.
#' @param scalingFactor The decoy scaling factor
#'
#' @returns A data frame
#' @export
#'
deScaler <- function(datTab, scalingFactor = the$decoyScalingFactor) {
  datTab.t <- datTab %>% filter(.data$Decoy == "Target")
  datTab.d <- datTab %>% filter(.data$Decoy == "Decoy")
  datTab.dd <- datTab %>% filter(.data$Decoy == "DoubleDecoy")
  dt.adjustment <- (nrow(datTab.d) / scalingFactor) - 2*(nrow(datTab.dd) / scalingFactor**2)
  dd.adjustment <- nrow(datTab.dd) / scalingFactor**2
  if (dt.adjustment < 0) {dt.adjustment <- 0}
  datTab.d <- datTab.d %>% slice(sample(nrow(datTab.d), dt.adjustment))
  datTab.dd <- datTab.dd %>% slice(sample(nrow(datTab.dd), dd.adjustment))
  bind_rows(datTab.t, datTab.d, datTab.dd)
}

#' Plots the score distributions of decoy and target crosslinked hits.
#'
#' @param datTab Parsed CLMS search results.
#' @param threshold Score threshold for classifying data.
#' @param classifier Classifier to use.
#' @param scalingFactor The decoy scaling factor.
#' @param separateFacets Whether to plot intraProtein and interProtein hits in separate facets
#' @param addLegend Whether to display the legend.
#' @param title Title to plot.
#'
#' @returns A ggplot object.
#' @export
#'
fdrPlots <- function(datTab,
                     threshold = 0,
                     classifier=.data$SVM.score,
                     scalingFactor = the$decoyScalingFactor,
                     separateFacets = T,
                     addLegend = T,
                     title = "FDR plot") {
  datTab <- ungroup(datTab)
  minValue = datTab %>% pull({{ classifier }}) %>% min(na.rm=T)
  minValue = floor(minValue)
  maxValue = datTab %>% pull({{ classifier }}) %>% max(na.rm=T)
  maxValue = ceiling(maxValue)
  stepSize = mmax((maxValue - minValue) / 100, 0.25)
  datTab <- deScaler(datTab, scalingFactor = scalingFactor)
  decCounts <- datTab %>% tally(.data$Decoy == "Decoy") %>% pull(.data$n)
  doubleCounts <- datTab %>% tally(.data$Decoy == "DoubleDecoy") %>% pull(.data$n)
  datTab <- datTab %>%
    mutate(Decoy = case_when(
      doubleCounts > decCounts ~ factor(Decoy, levels = c("Target", "DoubleDecoy", "Decoy")),
      doubleCounts <= decCounts ~ factor(Decoy, levels = c("Target", "Decoy", "DoubleDecoy"))
    ))
  fdr.plot <- datTab %>%
    ggplot(aes(x= {{ classifier }}, fill=.data$Decoy, alpha=.data$xlinkClass)) +
    geom_histogram(col="black", binwidth = stepSize, position="identity")
  if (is(threshold, "list")) {
    if (!is.null(threshold$interThresh) & !is.null(threshold$intraThresh)) {
      fdr.plot <- fdr.plot +
        geom_vline(data=filter(datTab, .data$xlinkClass=="intraProtein"), aes(xintercept=threshold$intraThresh), linetype = "dashed", col="green", size = 1.2) +
        geom_vline(data=filter(datTab, .data$xlinkClass=="interProtein"), aes(xintercept=threshold$interThresh), linetype = "dashed", col="green", size = 1.2)
    } else if (!is.null(threshold$globalThresh)) {
      fdr.plot <- fdr.plot +
        geom_vline(aes(xintercept=threshold$globalThresh), linetype = "dashed", col="green", size = 1.2)
    }
  } else {
    fdr.plot <- fdr.plot +
      geom_vline(xintercept = as.numeric(threshold), linetype = "dashed", col="green", size = 1.2)
  }
  if (separateFacets) {
    fdr.plot <- fdr.plot +
      facet_grid(rows=ggplot2::vars(.data$xlinkClass), scales="free_y")
  }
  fdr.plot <- fdr.plot +
    theme_bw() +
    xlim(minValue, maxValue) +
    scale_fill_manual(values=c("Target" = "lightblue",
                               "Decoy" = "salmon",
                               "DoubleDecoy" = "goldenrod1")) +
    scale_alpha_manual(values=c("interProtein" = 0.9, "intraProtein" = 0.4))
  if (!addLegend) {fdr.plot <- fdr.plot +
    theme(legend.position = "none") }
  fdr.plot <- fdr.plot +
    ggtitle(title)
  suppressWarnings(plot(fdr.plot))
}

#' Histogram of precursor mass errors.
#' Internal plotting function used by the Shiny interface to Touchstone
#'
#' @param massErrors precursor mass accuracy error (in ppm)
#' @param lowThresh mass accuracy threshold, low (in ppm)
#' @param highThresh mass accuracy threshold, high (in ppm)
#' @param lowPlotRange x-axis limit, low (in ppm)
#' @param highPlotRange x-axis limit, high (in ppm)
#' @keywords internal
#'
#' @returns No return value
#' @export
#'
massErrorPlot <- function(massErrors, lowThresh, highThresh, lowPlotRange, highPlotRange) {
  if (abs(lowPlotRange) <= 10 & abs(highPlotRange) <= 10) {
    binw <- 0.5
  } else {
    binw <- 1
  }
  hist(massErrors, col="slategray3",
       xlim=c(lowPlotRange, highPlotRange),
       breaks=seq(lowPlotRange, highPlotRange, binw),
       xlab = "Mass Error (ppm)",
       main = NA)
  title("Precursor Mass Deviation", adj=0)
  abline(v = c(lowThresh, highThresh), lwd = 2, lt = 1, col = "red")
}

#' Histogram of crosslinked distance distributions.
#' Internal plotting function used by the Shiny interface to Touchstone
#' Seems to be deprecated
#'
#' @param targetDists vector of c-alpha c-alpha distances
#' @param randomDists vector of random c-alpha c-alpha distances
#' @param threshold the distance violation threshold.
#' @keywords internal
#' @returns No return value
#' @export
#'
distancePlot <- function(targetDists, randomDists, threshold) {
  targetDists.size <- length(!is.na(targetDists))
  randomDists.size <- length(!is.na(randomDists))
  maxDist <- ceiling(max(c(randomDists, targetDists), na.rm=T)) + 5
  dists <- hist(targetDists, breaks=seq(-5, maxDist, 5), plot=F)
  r.dists <- hist(randomDists, breaks=seq(-5, maxDist, 5), plot=F)
  r.dists$counts <- r.dists$counts * (targetDists.size / randomDists.size)
  col1 <- rgb(141, 90, 151, alpha = 150, maxColorValue = 255)
  col2 <- rgb(184, 235, 208, alpha = 150, maxColorValue = 255)
  ylim.value <- ceiling(max(c(dists$counts, r.dists$counts)))
  plot(r.dists, col=col1, ylim=c(0, ylim.value),
       xlab = "C\u03B1-C\u03B1 Distance (\u00c5)",
       main = NA)
  title("Crosslink Violations", adj=0)
  plot(dists, col=col2, add=T)
  abline(v=threshold, lwd=2, lt=1, col="red")
  legend("topright", c("experimental", "random distribution"), fill = c(col2, col1),
         bty="n")
}

#' Histogram of crosslinked distance distributions.
#' distancePlot2 works with Touchstone datTab and differs as well as the Shiny interface.
#' @param datTab Parsed CLMS search results.
#' @param threshold the distance violation threshold.
#' @param binWidth binwidth for tabulating the histogram
#'
#' @returns A ggplot object
#' @export
#'
distancePlot2 <- function(datTab, threshold=35, binWidth=NA) {
  # col1 <- rgb(141, 90, 151, alpha = 150, maxColorValue = 255)
  # col2 <- rgb(184, 235, 208, alpha = 150, maxColorValue = 255)
  if (is.na(binWidth)) {
    maxDist <- max(datTab$distance, datTab$random.distance, na.rm=T)
    binWidth <- mmax(maxDist / 60, 2.5)
  }
  p <- datTab %>%
    pivot_longer(cols = c(.data$distance, .data$random.distance),
                 names_to = "distType",
                 values_to = "distance") %>%
    mutate(PDB = if_else(.data$distType == "random.distance", stringr::str_c(.data$PDB, ".random"), .data$PDB)) %>%
    ggplot(aes(.data$distance, fill=.data$PDB)) +
    geom_histogram(binwidth = binWidth, col="black", position="stack") +
    geom_vline(xintercept = threshold, col="red", size=1.25, linetype="dashed") +
    facet_grid(rows = "distType") +
    ggplot2::scale_fill_viridis_d(option="C") +
    # scale_fill_brewer(palette=4, type="seq", labels=abbreviate) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", size = 14)) +
    ggtitle("Crosslink Lengths") +
    xlab("C\u03B1-C\u03B1 Distance (\u00c5)")
  suppressWarnings(print(p))
}

#' Plot of number of crosslinked hits by FDR.
#' Internal plotting function used by the Shiny interface to Touchstone
#'
#' @param num.hits A data frame created by `generateErrorTable()`
#' @param threshold FDR threshold to plot
#'
#' @returns A ggplot object
#' @export
#' @seealso [generateErrorTable()]
numHitsPlot <- function(num.hits, threshold=0) {
  num.hits <- num.hits %>%
    filter(!is.infinite(.data$fdr), .data$total.hits >= 0.025 * max(.data$total.hits, na.rm=T)) %>%
    pivot_longer(cols=c(-.data$fdr,-.data$thresh),
                 names_to="crosslinkClass",
                 values_to="numHits") %>%
    filter(!is.na(.data$numHits))
  p <- num.hits %>% ggplot(aes(x=.data$fdr, y=.data$numHits)) +
    geom_path(aes(col=.data$crosslinkClass), na.rm=T, size=1.25) +
    theme_bw() +
    theme(legend.position="right") +
    scale_color_brewer(type="qual", palette="Paired") +
    xlab("FDR") +
    ylab("Num Hits") +
    ggtitle("Num Hits vs FDR") +
    theme(plot.title = element_text(face = "bold", size = 14),
    )
  p <- p + geom_vline(xintercept = threshold,
                      col="red", size=1.25)
  suppressWarnings(print(p))
}

#' Tabulates crosslink spectral counts (CSMs) by Protein Pair, typically for plotting.
#' Used by interactive, Shiny Touchstone interface.
#' @param datTab Parsed CLMS search results.
#'
#' @returns A data frame
#' @export
#'
summarizeProtData <- function(datTab) {
  datTab <- datTab %>% filter(.data$Decoy=="Target")
  datTab <- bestResPair(datTab)
  if (length(levels(datTab$Acc.1)) < 500) {
    datTab$Acc.1 <- forcats::fct_drop(datTab$Acc.1, only="decoy")
    datTab$Acc.2 <- forcats::fct_drop(datTab$Acc.2, only="decoy")
    matrixCounts <- datTab %>%
      group_by(.data$Acc.1, .data$Acc.2) %>%
      summarize(accCounts = sum(.data$numCSM), .groups = "drop") %>%
      complete(.data$Acc.1, .data$Acc.2) %>%
      unique() %>%
      pivot_wider(names_from = .data$Acc.2, values_from = .data$accCounts) %>%
      column_to_rownames("Acc.1") %>%
      as.matrix()
    matrixCounts[is.na(matrixCounts)] <- 0
    matrixCounts <- matrixCounts + t(matrixCounts)
    diag(matrixCounts) <- diag(matrixCounts) / 2
    accNames <- rownames(matrixCounts)
    matrixCounts <- clearAboveDiag(matrixCounts)
    matrixCounts <- as_tibble(matrixCounts)
    matrixCounts$Acc.1 <- accNames
    matrixCounts <- matrixCounts %>%
      pivot_longer(cols = -.data$Acc.1, names_to = "Acc.2", values_to = "counts") %>%
      mutate(counts = ifelse(.data$counts==0, NA, .data$counts))
    return(matrixCounts)
  } else {
    matrixCounts = NA
    return(matrixCounts)
  }
}

#' Tabulates crosslink spectral counts (CSMs) by Module Pair, typically for plotting.
#' Used by interactive, Shiny Touchstone interface.
#' @param datTab Parsed CLMS search results.
#' @param clearDiag Whether to de-duplicate Module Pairs, for plotting only half of symmetrical matrix.
#' @param modOrder Optional character vector showing order that modules should be arranged/plotted.
#'
#' @returns A data frame
#' @export
#'
summarizeModuleData <- function(datTab, clearDiag = F, modOrder = NULL) {
  if (!("Module.1" %in% names(datTab) & "Module.2" %in% names(datTab))) {
    stop("Modules must be assigned to create the Module Plot")
  }
  datTab <- datTab %>%
    filter(.data$Decoy=="Target") %>%
    mutate(Module.1 = forcats::fct_na_value_to_level(.data$Module.1, "other"),
           Module.2 = forcats::fct_na_value_to_level(.data$Module.2, "other"))
  if (is.null(modOrder)) {
    modOrder <- forcats::lvls_union(list(datTab$Module.1, datTab$Module.2))
    datTab[c("Module.1", "Module.2")] <- datTab %>%
      select(.data$Module.1, .data$Module.2) %>%
      as.list() %>%
      forcats::fct_unify(levels=modOrder) %>%
      bind_cols
  } else {
    if (!"other" %in% modOrder) {
      modOrder <- c(modOrder[1:length(modOrder)], "other")}
    datTab[c("Module.1", "Module.2")] <- datTab %>%
      select(.data$Module.1, .data$Module.2) %>%
      as.list() %>%
      forcats::fct_unify(levels = modOrder) %>%
      bind_cols
  }
  matrixCounts <- datTab %>%
    group_by(.data$Module.1, .data$Module.2) %>%
    summarize(modCounts = sum(.data$numCSM), .groups = "drop") %>%
    complete(.data$Module.1, .data$Module.2) %>%
    unique() %>%
    pivot_wider(names_from = .data$Module.2, values_from = .data$modCounts) %>%
    column_to_rownames("Module.1") %>%
    as.matrix()
  matrixCounts[is.na(matrixCounts)] <- 0
  matrixCounts <- matrixCounts + t(matrixCounts)
  diag(matrixCounts) <- diag(matrixCounts) / 2
  modNames <- rownames(matrixCounts)
  if (clearDiag) {matrixCounts <- clearAboveDiag(matrixCounts)}
  matrixCounts <- as_tibble(matrixCounts)
  matrixCounts$Module.1 <- modNames
  matrixCounts <- matrixCounts %>%
    pivot_longer(cols = -.data$Module.1, names_to = "Module.2", values_to = "counts") %>%
    mutate(counts = ifelse(.data$counts==0, NA, .data$counts),
           Module.1 = factor(.data$Module.1, levels = modOrder),
           Module.2 = factor(.data$Module.2, levels = modOrder))
  return(matrixCounts)
}

# developmental plotting function.  Doesn't really work yet.
.summarizeModuleDistance <- function(datTab, clearDiag = F, modOrder = NULL) {
  datTab <- datTab %>%
    filter(.data$Decoy=="Target") %>%
    mutate(Module.1 = forcats::fct_na_value_to_level(.data$Module.1, "other"),
           Module.2 = forcats::fct_na_value_to_level(.data$Module.2, "other"))
  if (is.null(modOrder)) {
    datTab[c("Module.1", "Module.2")] <- datTab %>%
      select(.data$Module.1, .data$Module.2) %>%
      as.list() %>%
      forcats::fct_unify() %>%
      bind_cols
  } else {
    datTab[c("Module.1", "Module.2")] <- datTab %>%
      select(.data$Module.1, .data$Module.2) %>%
      as.list() %>%
      forcats::fct_unify(levels = modOrder) %>%
      bind_cols
  }
  matrixCounts <- datTab %>%
    group_by(.data$Module.1, .data$Module.2) %>%
    nest() %>%
    mutate(vr = map_dbl(.data$data, function(x) violationRate(x, threshold=-10, classifier="SVM.score"))) %>%
    unnest(cols=c(.data$data)) %>%
    ungroup() %>%
    select(.data$Module.1, .data$Module.2, .data$vr) %>%
    complete(.data$Module.1, .data$Module.2) %>%
    unique() %>%
    pivot_wider(names_from = .data$Module.2, values_from = .data$vr) %>%
   column_to_rownames("Module.1") %>%
   as.matrix()
  # matrixCounts[is.na(vr)] <- 0
  # matrixCounts <- matrixCounts + t(matrixCounts)
  # diag(matrixCounts) <- diag(matrixCounts) / 2
   modNames <- rownames(matrixCounts)
  # if (clearDiag) {matrixCounts <- clearAboveDiag(matrixCounts)}
  matrixCounts <- as_tibble(matrixCounts)
  matrixCounts$Module.1 <- modNames
  matrixCounts <- matrixCounts %>%
    pivot_longer(cols = -.data$Module.1, names_to = "Module.2", values_to = "counts") %>%
    mutate(counts = ifelse(.data$counts==0, NA, .data$counts))
  return(matrixCounts)
}

#' Internal function used by some Touchstone plots to only display one half of the plot
#' in symmetrical square matrix representations of cross-linking data.
#'
#' @param sqMatrix A square matrix
#' @keywords internal
#' @returns A matrix
#' @export
#'
clearAboveDiag <- function(sqMatrix) {
  n <- dim(sqMatrix)[1]
  m <- dim(sqMatrix)[2]
  if (n != m) return("only works for square matrices")
  if (n == 1) return(sqMatrix)
  a <- 1 + (1:n-1)*n
  b <- (1:n-1)*(n+1)
  abovePositions <- unlist(purrr::map2(a, b, seq)[2:n])
  sqMatrix[abovePositions] <- NA
  return(sqMatrix)
}

#' Puts crosslink results into XiNet format for visualization.
#' see Combe et al, MCP 2015.
# 'https://doi.org/10.1074/mcp.O114.042259
#'
#' @param datTab Parsed CLMS search results.
#'
#' @returns A data frame
#' @export
#'
makeXiNetFile <- function(datTab) {
  datTab <- datTab %>%
    select(.data$SVM.score, .data$Acc.1, .data$Acc.2, .data$XLink.AA.1, .data$XLink.AA.2)
  names(datTab) <- c("Score", "Protein1", "Protein2", "LinkPos1", "LinkPos2")
  return(datTab)
}

#' Puts crosslink results into XVis format for visualization.
#' see Grimm et al, NAR 2017.
#' https://doi.org/10.1093/nar/gkv463
#'
#' @param datTab Parsed CLMS search results.
#'
#' @returns A data frame
#' @export
#'
XVisOutput <- function(datTab) {
  Protein1 <- pull(datTab, .data$Protein.1)
  Protein2 <- pull(datTab, .data$Protein.2)
  AbsPos1 <- pull(datTab, .data$XLink.AA.1)
  AbsPos2 <- pull(datTab, .data$XLink.AA.2)
  IDScore <- pull(datTab, .data$SVM.score)
  xVisXL <- data.frame("Protein1" = Protein1, "Protein2" = Protein2,
                       "AbsPos1" = AbsPos1, "AbsPos2"= AbsPos2, "Id-Score" = IDScore)
}

#' Internal convenience function for determining parameter ranges and bin sizes.
#'
#' @param x numeric parameter
#' @param base desired subdivision (integer)
#' @keywords internal
#' @returns numeric
#' @export
#'
mmax <- function(x, base=5) {
  base * ceiling(x/base)
}

#' Internal convenience function for determining parameter ranges and bin sizes.
#'
#' @param x numeric parameter
#' @param base desired subdivision (integer)
#' @keywords internal
#' @returns numeric
#' @export
#'
mmin <- function(x, base=5) {
  base * floor(x/base)
}

#' Reads the Protein Prospector parameter xml files.
#'
#' @param paramsFile xml file with the parameters submitted to Batch-Tag
#' @keywords internal
#' @returns An xml_document
#' @export
#'
readParamsFile <- function(paramsFile) {
  tryCatch(xml2::read_xml(paramsFile),
           error = function(e) {
             message(e)
             return(NA)
           }
  )
}

#' Plots the CSMs tabulated by Module-Pair in a Tile Plot
#'
#' @param datTab Parsed CLMS search results.
#' @param threshold Score threshold for classifying data.
#' @param title Plot title.
#' @param modBorders Numeric vector to specify optional borders between modules.
#' @param ... Parameters passed to `summarizeModuleData()`
#' @seealso [summarizeModuleData()]
#' @returns A ggplt object
#' @export
moduleTilePlot <- function(datTab, threshold=-100, title="Module Plot", modBorders = NULL, ...) {
  if (!("Module.1" %in% names(datTab) & "Module.2" %in% names(datTab))) {
    stop("Modules must be assigned to create the Module Plot")
    }
  tabulatedMods <- classifyDataset(datTab, threshold) %>%
    summarizeModuleData(...)
  gg <- tabulatedMods %>%
    tidyr::replace_na(list(counts=0)) %>%
    ggplot(aes(x=.data$Module.1, y=.data$Module.2, fill=.data$counts)) +
    ggplot2::geom_tile(color="black") +
    ggplot2::scale_fill_gradient(low="white", high="dodgerblue") +
    ggplot2::coord_fixed() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle(title)
  if (!is.null(modBorders)) {
    modBorders = modBorders + 0.5
    gg <- gg +
      geom_vline(xintercept = modBorders, size=1.1) +
      geom_hline(yintercept = modBorders, size=1.1)
  }
  gg <- gg +
    theme(axis.title.x=ggplot2::element_blank(),
          axis.title.y=ggplot2::element_blank())
  plot(gg)
}

#' Download STRING-DB scores for interProtein Crosslinks.
#'
#' Call the STRINGdb library and queries string for evidence supporting the
#' crosslinked PPIs. Attempts to automatically detect the species code from the
#' CLMS results, but this will only work for the most common research organisms.
#' Otherwise the user should provide the ncbiTaxonomy code:
#' `https://string-db.org/` under the `Oraganisms` link.
#' @param datTab Parsed CLMS search results.
#' @param ncbiTaxonomyCode NCBI format species code used by string-db
#' @returns A data frame
#' @export
getStringScores <- function(datTab, ncbiTaxonomyCode = NULL) {
  primarySpecies <- datTab %>%
    removeDecoys() %>%
    filter(.data$Score.Diff > 5) %>%
    count(.data$Species.1) %>%
    arrange(desc(.data$n)) %>%
    slice(1) %>%
    pull(.data$Species.1)
  tryCatch({
    if (is.null(ncbiTaxonomyCode)) {
      ncbiTaxonomyCode <- case_when(
        primarySpecies == "HUMAN" ~ 9606,
        primarySpecies == "MOUSE" ~ 10090,
        primarySpecies == "RAT" ~ 10116,
        primarySpecies == "ECOLI" ~ 511145,
        primarySpecies == "YEAST" ~ 4932,
        primarySpecies == "DROME" ~ 7227,
        primarySpecies == "ARATH" ~ 3702)
    }
    message(stringr::str_c("detected organism: ", primarySpecies, "\tncbi code:", ncbiTaxonomyCode))
  },
  error = function(cond) {
    message("Unknown species, please provide the ncbi taxonomy identifier")
    message("Original error message:")
    message(conditionMessage(cond))
    NA
  })
  string_db <- STRINGdb::STRINGdb$new(version = "12.0", network_type="full", link_data="combined_only",
                                      species = ncbiTaxonomyCode, score_threshold = 0, input_directory = "")
  datTab.inter <- datTab %>%
    filter(.data$xlinkClass == "interProtein")
  datTab.intra <- datTab %>%
    filter(.data$xlinkClass == "intraProtein")
  p.list.1 <- datTab.inter %>%
    removeDecoys() %>%
    pull(.data$Acc.1) %>%
    as.character()
  p.list.2 <- datTab.inter %>%
    removeDecoys() %>%
    pull(.data$Acc.2) %>%
    as.character()
  p.list <- unique(c(p.list.1, p.list.2))
  id.map <- string_db$map(data.frame(acc = p.list), "acc", removeUnmappedRows = F)
  datTab.intra$string.score <- NA
  datTab.inter <- datTab.inter %>%
    mutate(string.score = purrr::map2_dbl(.data$Acc.1, .data$Acc.2, function(x, y) {
      String.1 <- id.map %>% filter(.data$acc == x) %>% pull(.data$STRING_id)
      String.2 <- id.map %>% filter(.data$acc == y) %>% pull(.data$STRING_id)
      ppi <- string_db$get_interactions(c(String.1, String.2))
      if (length(ppi$combined_score)==0) {
        return(NA)
      } else {
        return(ppi$combined_score[[1]])
      }
    })
    )
  return(bind_rows(datTab.intra, datTab.inter))
}

#' Convenience function for working through the ribosome example dataset.
#'
#' @param path path to system file
#' @returns A character vector
#' @export
touchstone_example <- function(path = NULL) {
  if (is.null(path)) {
    dir(system.file("extdata", package = "touchstone"))
  } else {
    system.file("extdata", path, package = "touchstone", mustWork = TRUE)
  }
}
