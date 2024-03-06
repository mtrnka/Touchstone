calculateDecoyFractions <- function(datTab, scalingFactor=10) {
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

generateDecoyTable <- function(datTab, targetER=0.01, scalingFactor = 10, classifier=SVM.score) {
  #Note, only really tested with SVM.score. Use generateErrorTable if having issues with other classifiers.
  class = datTab %>% pull({{ classifier }})
  class.min = floor(min(class))
  class.max = ceiling(max(class))
  range.spacing <- abs(class.max - class.min) / 500
  class.range <- seq(class.min, class.max-range.spacing, by=range.spacing)
  decTable <- class.range %>% map_dfr(function(x) {
    datTab <- datTab %>% filter({{ classifier }} >= x, xlinkClass=="interProtein")
    numDecs <- calculateDecoyFractions(datTab, scalingFactor = scalingFactor)
    tibble(thresh = x, TT = numDecs["TT"], ftTT = numDecs["ftTT"], ffTT = numDecs["ffTT"])
  })
  decTable <- decTable %>%
    mutate(numIncorrect = ftTT + ffTT,
           fdr.exp = if_else(TT == 0, 0, numIncorrect / TT)
    )
  first0 <- which(decTable$fdr.exp < 0 & lag(decTable$fdr.exp) >= 0)
  if (length(first0) == 0) {
    first0 = nrow(decTable)
  } else {
    first0 <- first0[1]
  }
  minFDR <- decTable %>% slice(first0) %>% pull(fdr.exp)
  maxFDR <- max(decTable$fdr.exp[1:first0], na.rm=T)
  medFDR.post <- 0
  sdFDR.post <- 0.01 * abs(maxFDR - minFDR) / maxFDR
  decTable <- decTable %>%
    mutate(fdr.orig = fdr.exp,
           delta.fdr = abs(fdr.exp - lag(fdr.exp)),
           delta.roll = data.table::frollmean(delta.fdr, 10)
    )
  medDroll <- decTable %>% pull(delta.roll) %>% median(na.rm=T)
  decTable.problems <- which(decTable$delta.roll > medDroll)
  decTable.negs <- which(decTable$fdr.exp < 0)
  decTable.problems <- union(decTable.problems, decTable.negs) %>% sort
  decTable.problems <- decTable.problems[which(decTable.problems > first0)]
  decTable[decTable.problems, "fdr.exp"] <- rnorm(length(decTable.problems), medFDR.post, sdFDR.post)
  firstGuess <- which.min(abs(decTable$fdr.exp[1:first0] - targetER))
  decTable <- decTable %>% mutate(fdr.weights = case_when(
    abs(fdr.exp - targetER)/maxFDR < 0.1 ~ 5000,
    abs(fdr.exp - targetER)/maxFDR < 0.25 ~ 500,
    abs(fdr.exp - targetER)/maxFDR < 0.5 ~ 5,
    TRUE ~ 0.5))
  decTable.tailfit <- nlsLM(fdr.exp ~ I(maxFDR * A / (A + C**thresh)),
                            data=decTable,
                            start=list(A=1, C = exp(1)),
                            lower=c("A"=0, "C"=1),
                            weights=decTable$fdr.weights,
                            control=list(maxiter=100))
  decTable <- decTable %>%
    mutate(
      fdr = maxFDR*coef(decTable.tailfit)["A"] /
        (coef(decTable.tailfit)["A"] + coef(decTable.tailfit)["C"]**thresh),
    )
  decTable.plot <- decTable %>%
    ggplot(aes(x=thresh)) +
    geom_line(aes(y=fdr.orig), color="green", linewidth=1.5) +
    geom_line(aes(y=fdr.exp), color="violet", linewidth=1.5) +
    geom_line(aes(y=fdr), color="purple", linewidth=1.5) +
    theme_bw() + xlim(c(-4,4)) +
    geom_vline(xintercept = decTable$thresh[firstGuess], color="red", linetype="dashed") +
    geom_hline(yintercept = 0.01, color="red", linetype="dashed")
  plot(decTable.plot)
  return(decTable)
}

calculateFDR <- function(datTab, threshold=list(intraThresh=-100,interThresh=-100), classifier="SVM.score", scalingFactor=10) {
  datTab <- classifySeparateThresholds(datTab, separateThresh=threshold, classifier=classifier)
  decoyFractions <- calculateDecoyFractions(datTab, scalingFactor)
  fdr <- (decoyFractions["ffTT"] + decoyFractions["ftTT"]) / decoyFractions["TT"]
  names(fdr) <- NULL
  return(fdr)
}

calculateFDR.old <- function(datTab, threshold=-100, classifier="SVM.score", scalingFactor=10) {
  datTab <- datTab[datTab[classifier] >= threshold,]
  decoyFractions <- calculateDecoyFractions(datTab, scalingFactor)
  fdr <- (decoyFractions["ffTT"] + decoyFractions["ftTT"]) / decoyFractions["TT"]
  names(fdr) <- NULL
  return(fdr)
}

calculateSeparateFDRs <- function(datTab, thresholdIntra = -100, thresholdInter = -100,
                                  classifier="SVM.score", scalingFactor=10) {
  intraHits <- datTab[str_detect(datTab$xlinkClass, "intraProtein"),]
  interHits <- datTab[str_detect(datTab$xlinkClass, "interProtein"),]
  intraHits <- intraHits[intraHits[[classifier]] >= thresholdIntra, ]
  interHits <- interHits[interHits[[classifier]] >= thresholdInter, ]
  intraFractions <- calculateDecoyFractions(intraHits, scalingFactor)
  interFractions <- calculateDecoyFractions(interHits, scalingFactor)
  decoyFractions <- intraFractions + interFractions
  fdr <- (decoyFractions["ffTT"] + decoyFractions["ftTT"]) / decoyFractions["TT"]
  names(fdr) <- NULL
  return(fdr)
}

calculateHits <- function(datTab, threshold=-100, classifier="SVM.score") {
  datTab <- datTab[datTab[[classifier]] >= threshold, ]
  intra <- sum(str_count(datTab$xlinkClass, "intraProtein"))
  inter <- sum(str_count(datTab$xlinkClass, "interProtein"))
  return(c("thresh"=threshold, "intra"=intra, "inter"=inter))
}

calculateGT <- function(datTab, threshold=-100, classifier="SVM.score") {
  datTab <- datTab[datTab[[classifier]] >= threshold &
                     datTab$Decoy=="Target", ]
  gtTable <- table(datTab$groundTruth)
  if (is.na(gtTable["FALSE"])) {gtTable["FALSE"] <- 0}
  if (is.na(gtTable["TRUE"])) {gtTable["TRUE"] <- 0}
  gt <- gtTable["FALSE"] / sum(gtTable)
  return(gt)
}

violationRate <- function(datTab, threshold, classifier="Score.Diff", maxDistance = 35) {
  datTab <- datTab[datTab[[classifier]] >= threshold,]
  dists <- datTab$distance
  dists <- dists[!is.na(dists)]
  nAbove <- sum(dists > maxDistance)
  nBelow <- sum(dists <=maxDistance)
  return(nAbove / (nAbove + nBelow))
}

generateErrorTable <- function(datTab,
                               classifier="SVM.score",
                               errorFUN=calculateFDR.old,
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

  num.hits <- map_dfr(class.range, function(threshold) {
    calculateHits(removeDecoys(datTab), threshold=threshold, classifier=classifier)
  })
  num.hits$fdr <- error.rates
  num.hits$total <- num.hits$inter + num.hits$intra
  num.hits <- num.hits %>% filter(total > 0)
  return(num.hits)
}

findThresholdModelled <- function(datTab, targetER=0.05, minThreshold=-5, scalingFactor = 10, ...) {
  num.hits <- generateDecoyTable(datTab, targetER=targetER, scalingFactor = scalingFactor, ...)
  nearestScore <- which.min(abs(num.hits$fdr - targetER))
  threshold = num.hits$thresh[nearestScore]
  modelledFDR = num.hits$fdr[nearestScore]
  return(list("thresh"=threshold, "correspondingFDR" = modelledFDR))
}

findThresholdOld <- function(datTab, targetER=0.05, minThreshold=-5,
                             classifier="SVM.score", errorFUN=calculateFDR.old, ...) {
  num.hits <- generateErrorTable(datTab, classifier, errorFUN, ...)
  error.rates <- num.hits$fdr
  class.range <- num.hits$thresh
  error.diff <- (error.rates - targetER)
  error.cross <- logical(length(class.range))
  for (i in 1:(length(class.range) - 1)) {
    error.cross[i] <- error.diff[i] > 0 & error.diff[i+1] <= 0
  }
  if (sum(error.cross)==0 & error.diff[1] < 0) {
    threshold <- minThreshold
  } else if (sum(error.cross)==0 & error.diff[1] > 0) {
    threshold <- max(class.range)
    print("no acceptable threshold found")
  } else {
    firstCross <- which(error.cross)[1] + 1
    threshold <- class.range[firstCross]
  }
  if (threshold < minThreshold) {
    threshold <- minThreshold
  }
  return(list("globalThresh"=threshold, "correspondingFDR" = errorFUN(datTab, threshold, classifier, ...)))
}

findThreshold <- function(datTab, targetER=0.05, minThreshold=-5,
                          classifier="SVM.score", errorFUN=calculateFDR.old, ...) {
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

findSeparateThresholds <- function(datTab, targetER=0.01, minThreshold=-5,
                                   classifier="SVM.score", errorFUN=calculateFDR.old, ...) {
  interTab <- datTab %>%
    filter(xlinkClass=="interProtein")
  interThresh <- findThreshold(interTab, targetER, minThreshold, classifier, errorFUN, ...)[[1]]
  intraTab <- datTab %>%
    filter(xlinkClass=="intraProtein")
  intraThresh <- findThreshold(intraTab, targetER, minThreshold, classifier, errorFUN, ...)[[1]]
  return(list("intraThresh"=intraThresh, "interThresh"=interThresh))
}

findSeparateThresholdsModelled <- function(datTab, targetER=0.01, minThreshold=-5, scalingFactor=10) {
  interTab <- datTab %>%
    filter(xlinkClass=="interProtein")
  intraTab <- datTab %>%
    filter(xlinkClass=="intraProtein")
  interThresh <- findThresholdModelled(interTab, targetER, minThreshold, scalingFactor=scalingFactor)[[1]]
  intraThresh <- findThreshold(intraTab, targetER, minThreshold, classifier="SVM.score", scalingFactor=scalingFactor)[[1]]
  return(list("intraThresh"=intraThresh, "interThresh"=interThresh))
}

calculateSeparateHits <- function(datTab, thresholds = c("intraThresh"=-100, "interThresh"=-100), classifier="SVM.score") {
  datTab.intra <- datTab[datTab[[classifier]] >= thresholds["intraThresh"], ]
  intra <- sum(str_count(datTab.intra$xlinkClass, "intraProtein"))
  datTab.inter <- datTab[datTab[[classifier]] >= thresholds["interThresh"], ]
  inter <- sum(str_count(datTab.inter$xlinkClass, "interProtein"))
  return(tibble("thresh"=thresholds["interThresh"], "intra"=intra, "inter"=inter))
}

classifySingleThresholds <- function(datTab, singleThresh, classifier = "SVM.score") {
  datTab <- datTab %>%
    filter(!!rlang::sym(classifier) >= singleThresh$globalThresh)
  #  datTab <- addXlinkCount(datTab)
  return(datTab)
}

classifySeparateThresholds <- function(datTab, separateThresh, classifier = "SVM.score") {
  datTab <- datTab %>%
    filter(
      (xlinkClass == "interProtein" & !!rlang::sym(classifier) >= separateThresh$interThresh) |
        (xlinkClass == "intraProtein" & !!rlang::sym(classifier) >= separateThresh$intraThresh)
    )
  #  datTab <- addXlinkCount(datTab)
  return(datTab)
}

removeDecoys <- function(datTab) {
  datTabR <- datTab[datTab$Decoy=="Target",]
  return(datTabR)
}

thresholdResults <- function(datTab, threshold, classifier="Score.Diff") {
  datTabR <- datTab[datTab[[classifier]] >= threshold,]
  return(datTabR)
}
