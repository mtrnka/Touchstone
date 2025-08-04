#' Calculate FDR for a dataset above inter and intra-protein score thresholds
#'
#' `calculateFDR()` estimates the FDR for a dataset based on counting the number
#' of decoy hits. Treats inter- and intra-protein hits separately. If the
#' parameter `separateThresh` is provided it first thresholds the data.
#'
#' @param datTab Parsed CLMS search results.
#' @param threshold A list or a numeric of length 1. List must contain either `globalThresh` or both `interThresh` and `intraThresh` numeric elements.
#' @param classifier Column name in `datTab` used as the classifier to use to rank hits.
#' @param scalingFactor An integer k. The multiple by which decoy DB is larger than target DB
#' @return A numeric value
#' @export
#'
calculateFDR <- function(datTab,
                         threshold=list(intraThresh=-100,interThresh=-100),
                         classifier=.data$SVM.score,
                         scalingFactor=the$decoyScalingFactor) {
  classifier.quo <- enquo(classifier)
  datTab <- classifyDataset(datTab, threshold=threshold, classifier=!!classifier.quo)
  decoyFractions <- calculateDecoyFractions(datTab, scalingFactor)
  fdr <- (decoyFractions["ffTT"] + decoyFractions["ftTT"]) / decoyFractions["TT"]
  names(fdr) <- NULL
  return(fdr)
}

#' Calculate FDR for a dataset, above a single score threshold
#'
#' `calculateFDR.unseparated()` estimates the FDR for a dataset based on counting the number
#' of decoy hits. Does not separate inter- and intra-protein hits.
#' @param datTab Parsed CLMS search results.
#' @param threshold Score threshold.
#' @param classifier Column name in `datTab` used as the classifier to use to rank hits.
#' @param scalingFactor An integer k. The multiple by which decoy DB is larger than target DB
#' @return A numeric value
#' @keywords internal
#' @seealso [calculateFDR()]
#' @export
#'
calculateFDR.unseparated <- function(datTab, threshold=-100, classifier="SVM.score", scalingFactor=the$decoyScalingFactor) {
  datTab <- datTab[datTab[classifier] >= threshold,]
  decoyFractions <- calculateDecoyFractions(datTab, scalingFactor)
  fdr <- (decoyFractions["ffTT"] + decoyFractions["ftTT"]) / decoyFractions["TT"]
  names(fdr) <- NULL
  return(fdr)
}

#' Calculate the violation rate
#'
#' @param datTab parsed CLMS search results
#' @param threshold single score threshold
#' @param classifier which classifier to use to rank hits
#' @param maxDistance Calpha-Calpha distance which is considered a violation (Å)
#' @return A numeric value
#' @export
#'
violationRate <- function(datTab, threshold, classifier="Score.Diff", maxDistance = 35) {
  datTab <- datTab[datTab[[classifier]] >= threshold,]
  dists <- datTab$distance
  dists <- dists[!is.na(dists)]
  nAbove <- sum(dists > maxDistance)
  nBelow <- sum(dists <=maxDistance)
  return(nAbove / (nAbove + nBelow))
}

#' Calculate actual FDR when analyzing benchmarking datasets where ground truth is established.
#'
#' Requires the dataset to have a logical column named `groundTruth`. This is typically
#' for testing against benchmarking datasets. This particular function was written
#' to reproduce a specific analysis where each `Fraction` contained a replicate.
#' Will turn this into a more generally useful function at some point.
#'
#' @param datTab Parsed CLMS search results.
#' @return A data frame
#' @export
#'
calculateGroundTruthFDR <- function(datTab) {
  results <- datTab %>%
    filter(.data$Decoy=="Target") %>%
    group_by(.data$Fraction, .data$groundTruth) %>%
    count() %>%
    pivot_wider(names_from=.data$groundTruth, values_from=.data$n)
  names(results)[which(names(results)=="TRUE")] <- "Correct"
  names(results)[which(names(results)=="FALSE")] <- "Incorrect"
  results[is.na(results)] <- 0
  results <- results %>%
    mutate(FDR = round((100 * .data$Incorrect)/(.data$Incorrect + .data$Correct),2))
  means <- tibble(Fraction="means:",
                  Incorrect=mean(results$Incorrect, na.rm=T),
                  Correct=mean(results$Correct, na.rm=T),
                  FDR=mean(results$FDR, na.rm=T))
  sds <- tibble(Fraction="sds:",
                Incorrect=sd(results$Incorrect, na.rm=T),
                Correct=sd(results$Correct, na.rm=T),
                FDR=sd(results$FDR, na.rm=T))
  results <- bind_rows(results, means, sds)
  return(results)
}

#' Estimates fraction of TT hits with one or both peptides incorrectly identified.
#'
#' Internal function called by other functions in `fdr.R`
#'
#' @param datTab Parsed CLMS search results.
#' @param scalingFactor An integer k. The multiple by which decoy DB is larger than target DB
#' @return A numeric vector
#' @keywords internal
calculateDecoyFractions <- function(datTab, scalingFactor=the$decoyScalingFactor) {
  fdrTable <- table(datTab$Decoy)
  if (is.na(fdrTable["DoubleDecoy"])) {fdrTable["DoubleDecoy"] <- 0}
  if (is.na(fdrTable["Decoy"])) {fdrTable["Decoy"] <- 0}
  if (is.na(fdrTable["Target"])) {fdrTable["Target"] <- 0}
  ffTT <- fdrTable[["DoubleDecoy"]] / (scalingFactor ** 2)
  ftTT <- (fdrTable[["Decoy"]] / scalingFactor) -
    (2 * fdrTable[["DoubleDecoy"]] / (scalingFactor ** 2))
  TT <- fdrTable[["Target"]]
  if (ftTT < 0) {
    ffTT <- 0.5 * fdrTable[["Decoy"]] / scalingFactor
    ftTT <- 0
  }
  return(c("TT"=TT, "ftTT"=ftTT, "ffTT"=ffTT))
}

#' Model dependence of FDR on scoring function
#'
#' Internal function called by other functions in `fdr.R`
#'
#' @param datTab Parsed CLMS search results.
#' @param targetER desired error rate (typically FDR, but can be something else like violation rate or ground truth based error)
#' @param scalingFactor An integer k. The multiple by which decoy DB is larger than target DB.
#' @param classifier Column name in `datTab` used as the classifier to use to rank hits.
#' @param plot Show plot?
#' @return A data frame
#' @keywords internal
generateDecoyTable <- function(datTab,
                               targetER=0.01,
                               scalingFactor = the$decoyScalingFactor,
                               classifier=.data$SVM.score,
                               plot=T) {
  class = datTab %>% pull({{classifier}})
  class.min = floor(min(class))
  class.max = ceiling(max(class))
  range.spacing <- abs(class.max - class.min) / 500
  class.range <- seq(class.min, class.max-range.spacing, by=range.spacing)
  decTable <- class.range %>%
    purrr::map_dfr(function(x) {
      datTab <- datTab %>% filter({{classifier}} >= x, .data$xlinkClass=="interProtein")
      numDecs <- calculateDecoyFractions(datTab, scalingFactor = scalingFactor)
      tibble(thresh = x, TT = numDecs["TT"], ftTT = numDecs["ftTT"], ffTT = numDecs["ffTT"])
    })
  decTable <- decTable %>%
    mutate(numIncorrect = .data$ftTT + .data$ffTT,
           fdr.exp = if_else(.data$TT == 0, 0, .data$numIncorrect / .data$TT)
    )
  first0 <- which(decTable$fdr.exp < 0 & dplyr::lag(decTable$fdr.exp) >= 0)
  if (length(first0) == 0) {
    first0 = nrow(decTable)
  } else {
    first0 <- first0[1]
  }
  minFDR <- decTable %>% slice(first0) %>% pull(.data$fdr.exp)
  maxFDR <- max(decTable$fdr.exp[1:first0], na.rm=T)
  medFDR.post <- 0
  sdFDR.post <- 0.01 * abs(maxFDR - minFDR) / maxFDR
  decTable <- decTable %>%
    mutate(fdr.orig = .data$fdr.exp,
           delta.fdr = abs(.data$fdr.exp - dplyr::lag(.data$fdr.exp)),
           delta.roll = data.table::frollmean(.data$delta.fdr, 10)
    )
  medDroll <- decTable %>% pull(.data$delta.roll) %>% median(na.rm=T)
  decTable.problems <- which(decTable$delta.roll > medDroll)
  decTable.negs <- which(decTable$fdr.exp < 0)
  decTable.problems <- union(decTable.problems, decTable.negs) %>% sort
  decTable.problems <- decTable.problems[which(decTable.problems > first0)]
  decTable[decTable.problems, "fdr.exp"] <- stats::rnorm(length(decTable.problems), medFDR.post, sdFDR.post)
  firstGuess <- which.min(abs(decTable$fdr.exp[1:first0] - targetER))
  decTable <- decTable %>%
    mutate(fdr.weights = case_when(
      abs(.data$fdr.exp - targetER)/maxFDR < 0.1 ~ 10000,
      abs(.data$fdr.exp - targetER)/maxFDR < 0.25 ~ 1,
      abs(.data$fdr.exp - targetER)/maxFDR < 0.5 ~ .1,
      TRUE ~ 0.5))
  decTable.tailfit <- minpack.lm::nlsLM(fdr.exp ~ I(maxFDR/ (1 + exp(A*(decTable$thresh - C)))),
                                        data=decTable,
                                        start=list(A=1, C = 0),
                                        weights=decTable$fdr.weights,
                                        control=list(maxiter=100))
  decTable <- decTable %>%
    mutate(
      fdr = maxFDR /
        (1 + exp(stats::coef(decTable.tailfit)["A"] * (decTable$thresh - stats::coef(decTable.tailfit)["C"]))),
    )
  decTable.plot <- decTable %>%
    ggplot(aes(x=.data$thresh)) +
    geom_line(aes(y=.data$fdr.orig), color="green", linewidth=1.5) +
    geom_line(aes(y=.data$fdr.exp), color="violet", linewidth=1.5) +
    geom_line(aes(y=.data$fdr), color="purple", linewidth=1.5) +
    theme_bw() + xlim(c(-8,4)) +
    geom_vline(xintercept = decTable$thresh[firstGuess], color="red", linetype="dashed") +
    geom_hline(yintercept = targetER, color="red", linetype="dashed")
  if (plot) {suppressWarnings(plot(decTable.plot))}
  return(decTable)
}

#' Calculates number of interProtein and intraProtein hits above a single threshold.
#'
#' Internal function called by other functions in `fdr.R`
#'
#' @param datTab Parsed CLMS search results.
#' @param threshold Single score threshold
#' @param classifier Column name in `datTab` used as the classifier to use to rank hits.
#' @return A numeric vector
#' @keywords internal
calculateHits <- function(datTab, threshold=-100, classifier="SVM.score") {
  datTab <- datTab[datTab[[classifier]] >= threshold, ]
  intra <- sum(stringr::str_count(datTab$xlinkClass, "intraProtein"))
  inter <- sum(stringr::str_count(datTab$xlinkClass, "interProtein"))
  return(c("thresh"=threshold, "intra"=intra, "inter"=inter))
}


#' This is an older version of generateDecoyTable, which is still called a number times.
#'
#' @param datTab Parsed CLMS search results.
#' @param classifier Column name in `datTab` used as the classifier to use to rank hits.
#' @param errorFUN error function used to generate decoy table (eg, FDR or violation rate, or GT)
#' @param ... additional params passed by calling functions
#' @param class.min minimum value of classifier to consider
#' @param class.max maximum value of classifier to consider
#' @return A data frame
#' @keywords internal
generateErrorTable <- function(datTab,
                               classifier="SVM.score",
                               errorFUN=calculateFDR.unseparated,
                               ...,
                               class.min = NA,
                               class.max = NA) {
  if (is.na(class.min)) {class.min = floor(min(datTab[[classifier]]))}
  if (is.na(class.max)) {class.max = ceiling(max(datTab[[classifier]]))}
  range.spacing <- abs(class.max - class.min) / 500
  class.range <- seq(class.min, class.max-range.spacing, by=range.spacing)
  error.rates <- unlist(lapply(class.range, function(threshold) {
    errorFUN(datTab, threshold=threshold, classifier=classifier, ...)
  }))
  error.rates[is.na(error.rates)] <- 0
  #plot(error.rates ~ class.range, type="l")

  num.hits <- purrr::map_dfr(class.range, function(threshold) {
    calculateHits(removeDecoys(datTab), threshold=threshold, classifier=classifier)
  })
  num.hits$fdr <- error.rates
  num.hits$total.hits <- num.hits$inter + num.hits$intra
  num.hits <- num.hits %>% filter(.data$total.hits > 0)
  return(num.hits)
}

#' Determine SVM.score threshold to give desired FDR, using modeled FDR function
#'
#' `findThreshold()` uses calculated values at each threshold. Differs from
#' `findThresholdModelled()` which first models the error functions. Neither
#' of these functions separate the data into different classes. Can be used on
#' any summarization level.
#'
#' @param datTab Parsed CLMS search results.
#' @param targetER Desired FDR.
#' @param minThreshold Minimum acceptable SVM.score threshold.
#' @param scalingFactor An integer k. The multiple by which decoy DB is larger than target DB.
#' @param ... Additional params passed by calling function.
#' @return A list with the threshold and calculated FDR at this threshold
#' @seealso [findThreshold()], [findSeparateThresholds()]
#' @export
#'
findThresholdModelled <- function(datTab, targetER=0.01, minThreshold=-5, scalingFactor = the$decoyScalingFactor, ...) {
  num.hits <- generateDecoyTable(datTab, targetER=targetER, scalingFactor = scalingFactor, ...)
  nearestScore <- which.min(abs(num.hits$fdr - targetER))
  threshold = num.hits$thresh[nearestScore]
  modelledFDR = num.hits$fdr[nearestScore]
  return(list("globalThresh"=threshold, "correspondingFDR" = modelledFDR))
}

#' Determine SVM.score threshold to give desired FDR
#'
#' `findThreshold()` uses calculated values at each threshold. Differs from
#' `findThresholdModelled()` which first models the error functions. Neither
#' of these functions separate the data into different classes. Can be used on
#' any summarization level.
#'
#' @param datTab Parsed CLMS search results.
#' @param targetER Desired error rate (typically FDR).
#' @param minThreshold Minimum acceptable threshold.
#' @param classifier Column name in `datTab` used as the classifier to use to rank hits.
#' @param errorFUN Function used to estimate error.
#' @param ... Additional params passed by calling function. Typically scalingFactor.
#' @return A list with the threshold and calculated FDR at this threshold
#' @seealso [findThresholdModelled()], [findSeparateThresholds()], [findSeparateThresholds()],
#' @export
#'
findThreshold <- function(datTab, targetER=0.01, minThreshold=-5,
                          classifier="SVM.score", errorFUN=calculateFDR.unseparated, ...) {
  num.hits <- generateErrorTable(datTab, classifier, errorFUN, ...)
  error.rates <- num.hits$fdr
  class.range <- num.hits$thresh
  if (max(error.rates, na.rm=T) < targetER) {
    threshold = class.range[1]
  } else if (targetER < 0.15) {
    low.index <- which.min(abs(error.rates - 0.25))
    low.value <- class.range[low.index]
    high.index <- which.min(abs(error.rates[low.index:length(error.rates)]))
    high.value <- class.range[low.index + high.index]
    num.hits <- generateErrorTable(datTab, class.min = low.value, class.max = high.value,
                                   classifier = classifier, errorFUN = errorFUN, ...)
    error.rates <- num.hits$fdr
    class.range <- num.hits$thresh
  }
  threshold = class.range[which.min(abs(error.rates - targetER))]
  return(list("globalThresh"=threshold, "correspondingFDR" = errorFUN(datTab, threshold, classifier, ...)))
}

#' Determine separate SVM.score thresholds for inter and intra-protein hits that give the desired FDR
#'
#' `findSeparateThresholds()` uses calculated values at each threshold. Differs from
#' `findSeparateThresholdsModelled()` which first models the error functions. Both
#' of these functions separate inter-protein and intra-protein crosslinks into separate
#' classes and find two different thresholds that give the desired error rate.
#' Can be used on any summarization level.
#'
#' @param datTab Parsed CLMS search results.
#' @param targetER Desired error rate (typically FDR).
#' @param minThreshold Minimum acceptable score threshold.
#' @param classifier Column name in `datTab` used as the classifier to use to rank hits.
#' @param errorFUN Function used to estimate error.
#' @param ... Additional params passed by calling function. Typically scalingFactor.
#' @return A list of inter and intra-protein thresholds to give the desired error rate.
#' @seealso [findThreshold()], [findThresholdModelled()], [findSeparateThresholds()]
#' @export
#'
findSeparateThresholds <- function(datTab, targetER=0.01, minThreshold=-5,
                                   classifier="SVM.score", errorFUN=calculateFDR.unseparated, ...) {
  interTab <- datTab %>%
    filter(.data$xlinkClass=="interProtein")
  interThresh <- findThreshold(interTab, targetER, minThreshold, classifier, errorFUN, ...)[[1]]
  intraTab <- datTab %>%
    filter(.data$xlinkClass=="intraProtein")
  intraThresh <- findThreshold(intraTab, targetER, minThreshold, classifier, errorFUN, ...)[[1]]
  return(list("intraThresh"=intraThresh, "interThresh"=interThresh))
}

#' Determine separate SVM.score thresholds for inter and intra-protein hits that give the desired FDR using modeled FDR function.
#'
#' `findSeparateThresholds()` uses calculated values at each threshold. Differs from
#' `findSeparateThresholdsModelled()` which first models the error functions. Both
#' of these functions separate inter-protein and intra-protein crosslinks into separate
#' classes and find two different thresholds that give the desired error rate.
#' Can be used on any summarization level.
#'
#' @param datTab Parsed CLMS search results.
#' @param targetER Desired FDR.
#' @param minThreshold Minimum acceptable SVM.score threshold.
#' @param scalingFactor k, multiple by which which decoy DB is larger than target DB
#' @param plot Show decoy table plot?
#' @return A list of inter and intra-protein thresholds to give the desired error rate.
#' @seealso [findThreshold()], [findThresholdModelled()], [findSeparateThresholds()]
#' @export
#'
findSeparateThresholdsModelled <- function(datTab, targetER=0.01, minThreshold=-5, scalingFactor=the$decoyScalingFactor, plot=T) {
  interTab <- datTab %>%
    filter(.data$xlinkClass=="interProtein")
  intraTab <- datTab %>%
    filter(.data$xlinkClass=="intraProtein")
  interThresh <- findThresholdModelled(interTab, targetER, minThreshold, scalingFactor=scalingFactor, plot=plot)[[1]]
  intraThresh <- findThreshold(intraTab, targetER, minThreshold, classifier="SVM.score", scalingFactor=scalingFactor)[[1]]
  return(list("intraThresh"=intraThresh, "interThresh"=interThresh))
}

#' Classify CLMS datasets
#'
#' The `threshold` (a list) argument determines the behavior of `classifyDataset`.
#' If threshold contains a `globalThresh` element than `classifySingleThreshold()`
#' is called. If `threshold` contains `interThresh` and `intraThresh` elemenst than
#' inter-protein and intra-protein hits are treated separately and `classifySeparateThresholds()`
#' is called. `threshold` can also be a numeric value which call `classifySingleThreshold()`
#' @param datTab Parsed CLMS search results.
#' @param threshold A list or a numeric of length 1. List must contain either `globalThresh` or both `interThresh` and `intraThresh` numeric elements.
#' @param classifier Column name in `datTab` used as the classifier to use to rank hits.
#' @seealso [classifySeparateThresholds()], [classifySingleThreshold()]
#' @return A data frame
#' @export
#'
classifyDataset <- function(datTab, threshold=list(), classifier = "SVM.score") {
  classifier <- enquo(classifier)
  tryCatch(
    if (is.numeric(threshold) & length(threshold) == 1){
      datTab <- classifySingleThreshold(datTab, singleThresh = threshold, classifier = {{ classifier }})
    } else if (is.list(threshold) & !is.null(threshold$globalThresh)) {
      datTab <- classifySingleThreshold(datTab, singleThresh = threshold, classifier = {{ classifier }})
    } else if (is.list(threshold) & !is.null(threshold$intraThresh) & !is.null(threshold$interThresh)) {
      datTab <- classifySeparateThresholds(datTab, separateThresh = threshold, classifier = {{ classifier }})
    },
    error = function(e) {
      message(stringr::str_c("threshold must either be a single, numeric value or a list
                             with named interThresh and intraThresh values", as.character(e)))
      return(NA)
    }
  )
  return(datTab)
}

#' Classify CLMS dataset using a single threshold
#'
#' @param datTab Parsed CLMS search results.
#' @param singleThresh Either a list that contains `globalThresh` numeric elements or a single numeric value.
#' @param classifier Column name in `datTab` used as the classifier to use to rank hits.
#' @seealso [classifySeparateThresholds()], [classifyDataset()]
#' @return A data frame
#' @export
#'
classifySingleThreshold <- function(datTab, singleThresh, classifier = "SVM.score") {
  classifier <- enquo(classifier)
  if (is.numeric(singleThresh) & length(singleThresh==1)) {
    datTab <- datTab %>%
      filter({{classifier}} >= singleThresh)
  } else if (is.list(singleThresh) & !is.null(singleThresh$globalThresh)) {
    datTab <- datTab %>%
      filter({{classifier}} >= singleThresh$globalThresh)
  }
  #  datTab <- addXlinkCount(datTab)
  return(datTab)
}

#' Classify CLMS datasets using separate interProtein and intraProtein thresholds
#'
#' @param datTab Parsed CLMS search results.
#' @param separateThresh interProtein and intraProtein score thresholds as list(intraThresh = x, interThresh = y)
#' @param classifier Column name in `datTab` used as the classifier to use to rank hits.
#' @seealso [classifyDataset()], [classifySingleThreshold()]
#' @return A data frame
#' @export
#'
classifySeparateThresholds <- function(datTab, separateThresh, classifier = "SVM.score") {
  classifier = ensym(classifier)
  datTab <- datTab %>%
    filter(
      (.data$xlinkClass == "interProtein" & {{classifier}} >= separateThresh$interThresh) |
        (.data$xlinkClass == "intraProtein" & {{classifier}} >= separateThresh$intraThresh)
    )
  #  datTab <- addXlinkCount(datTab)
  return(datTab)
}

#' Remove decoy hits from tables
#'
#' Convenience function to remove decoys, typically before producing a report.
#' @param datTab Parsed CLMS search results.
#' @return A data frame
#' @export
#'
removeDecoys <- function(datTab) {
  datTabR <- datTab[datTab$Decoy=="Target",]
  return(datTabR)
}

#' Internal function called during hyperparameter optimization of SVM models.
#' @param datTab Parsed CLMS search results.
#' @param classifier Column name in `datTab` used as the classifier to use to rank hits.
#' @param errorFUN error function used to generate decoy table (eg, FDR or violation rate, or GT)
#' @param ... additional params passed by calling functions
#' @param class.min minimum value of classifier to consider
#' @param class.max maximum value of classifier to consider
#' @return A data frame
#' @keywords internal
#' @export
generateErrorTable.sep <- function(datTab,
                               classifier="SVM.score",
                               errorFUN=calculateFDR.unseparated,
                               ...,
                               class.min = NA,
                               class.max = NA) {
  if (is.na(class.min)) {class.min = floor(min(datTab[[classifier]]))}
  if (is.na(class.max)) {class.max = ceiling(max(datTab[[classifier]]))}
  range.spacing <- abs(class.max - class.min) / 500
  class.range <- seq(class.min, class.max-range.spacing, by=range.spacing)
  error.rates.inter <- unlist(lapply(class.range, function(threshold) {
    errorFUN(datTab[datTab$xlinkClass=="interProtein",], threshold=threshold, classifier=classifier, ...)
  }))
  error.rates.inter[is.na(error.rates.inter)] <- 0
  error.rates.intra <- unlist(lapply(class.range, function(threshold) {
    errorFUN(datTab[datTab$xlinkClass=="intraProtein",], threshold=threshold, classifier=classifier, ...)
  }))
  error.rates.intra[is.na(error.rates.intra)] <- 0
  num.hits <- purrr::map_dfr(class.range, function(threshold) {
    calculateHits(removeDecoys(datTab), threshold=threshold, classifier=classifier)
  })
  num.hits$fdr.intra <- error.rates.intra
  num.hits$fdr.inter <- error.rates.inter
  num.hits$total.hits <- num.hits$inter + num.hits$intra
  num.hits <- num.hits %>% filter(.data$total.hits > 0)
  return(num.hits)
}

#' Convenience function to display the confusion matrix.
#'
#' @param datTab Parsed CLMS search results.
#' @param threshold A list or a numeric of length 1. List must contain either `globalThresh` or both `interThresh` and `intraThresh` numeric elements.
#' @parm ... passed down to `classifyDataset()`
#' @returns A data frame
#' @export
countDecoys <- function(datTab, threshold=NULL, ...) {
  if (!is.null(threshold)) datTab <- classifyDataset(datTab, threshold, ...)
  datTab %>%
    deScaler(scalingFactor = the$decoyScalingFactor) %>%
    group_by(.data$xlinkClass, .data$Decoy) %>%
    count() %>%
    pivot_wider(names_from="Decoy", values_from="n")
}
