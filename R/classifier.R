#' Trains an SVM classifier for CLMS datasets. Performs basic feature selection with
#' different sets of features used for small database searches and larger one.
#' Performs some prefilitering of the data that optimizes performance based on database
#' search size. Performs hyperparamater optimization in a grid search using `tuneSVM()`
#' and chooses the model which produces the most inter-protein crosslinked residue pairs
#' in the range of 0.01 and 0.05 FDR.
#'
#' @param datTab Parsed CLMS search results
#' @param params Character vector specifying names of the features in `datTab` used to train model.
#' @param scoreName Name for the new scoring function.
#' @param scalingFactor An integer k. The multiple by which decoy DB is larger than target DB
#' @param targetER Desired FDR for classification of CSMs
#' @param sampleNo Size of the training dataset (integer).
#' @param cost_values Numeric vector of cost values used for hyperparameter tuning of the SVM model
#' @param gamma_values Numeric vector of gamma values used for hyperparameter tuning of the SVM model
#' @param sd_values Numeric vector of Score Diff values to use for prefilitering optimiziation.
#' @seealso [tuneSVM()], [tuneSVM.helper()], [buildSVM()]
#' @returns A list containing all of the SVM models at different cost and gamma values
#' as well as a summary table.
#' @export
trainCrosslinkScore <- function(datTab,
                                params = NULL,
                                scoreName="SVM.score",
                                scalingFactor = the$decoyScalingFactor,
                                targetER = 0.01,
                                sampleNo = 20000,
                                cost_values = c(1, 5, 10),
                                gamma_values = c(0.01, 0.05, 0.1),
                                sd_values = c(0,5,10,15,20)
) {
  datTab <- ungroup(datTab)
  # feature selection
  plausibleHits <- datTab %>%
    filter(.data$Decoy == "Target", .data$Score.Diff > 10, .data$numCSM > 1, .data$xlinkClass == "intraProtein") %>%
    group_by(.data$Acc.1) %>%
    count %>%
    arrange(desc(.data$n))
  if (is.null(params)) {
    params <- case_when(
      nrow(plausibleHits) <= 1 ~ list(params.best.nop),
      plausibleHits$n[1] > 100 * plausibleHits$n[2] ~ list(params.best.nop),
      !("Perc.Bond.Cleavage.1" %in% names(datTab)) ~ list(params.noPercBond),
      TRUE ~ list(params.best)
    ) %>% unlist()
  }
  # prefiltering
  if (length(unique(pull(removeDecoys(datTab), .data$Acc.1, .data$Acc.2))) > 1000) {
    preFiltered.dts <- sd_values %>%
      purrr::map(function(sd) {
        datTab.pre <- datTab %>%
          filter(.data$Score.Diff >= sd)
        tuneSVM.helper(datTab=datTab.pre,
                       targetER=targetER,
                       params=params,
                       scoreName=scoreName,
                       scalingFactor = scalingFactor,
                       sampleNo = sampleNo,
                       cost=10, gamma=0.1, kernel="radial")
      })
    preFilter.summary <- preFiltered.dts %>%
      purrr::imap_dfr(function(x, i) {
        data.frame("index" = i, "sd.thresh" = x$sd.thresh,
                    "interInt" = x$interInt, "interHits" = x$interHits)
        }) %>%
      arrange(desc(.data$interInt))
    bestPreFilter <- preFilter.summary %>%
      filter(dplyr::between(.data$interInt, 0.95 * max(.data$interInt), max(.data$interInt))) %>%
      pull(.data$sd.thresh) %>%
      min()
    #    return(list(preFiltered.dts, preFilter.summary, bestPreFilter))
    datTab <- datTab %>%
      filter(.data$Score.Diff >= bestPreFilter)
  }
  # hyperparamater optimziation
  tuned <- tuneSVM(datTab,
                   params=params,
                   scoreName=scoreName,
                   scalingFactor = scalingFactor,
                   targetER = targetER,
                   sampleNo = sampleNo,
                   cost_values = cost_values,
                   gamma_values = gamma_values)
  tuned.parse <- tuned %>%
    purrr::imap_dfr(function(x,i) {
      data.frame("index" = i, "cost" = x$cost, "gamma" = x$gamma,
                 "interInt" = x$interInt, "interHits" = x$interHits)
      }) %>%
    arrange(desc(.data$interInt))
  tuned.plot <- tuned %>%
    map_dfr(function(x) {
      df <- x$errorTable
      df <- df %>% mutate(
        kernel=x$kernel,
        cost=x$cost,
        gamma=x$gamma)
    }) %>%
    filter(.data$fdr.inter <= 0.05) %>%
    filter(.data$inter >= 0.05 * max(.data$inter)) %>%
    ggplot(aes(x=.data$fdr.inter, y=.data$inter, col=as.factor(.data$cost))) +
    geom_line(linewidth=1.2) +
    theme_bw() +
    xlim(0,0.05) +
    geom_vline(xintercept = targetER, color="red") +
    ggplot2::scale_color_viridis_d(option="C") +
    facet_grid(rows=ggplot2::vars(gamma), scales="free_y")
  bestModelIndex <- tuned.parse %>%
    filter(dplyr::between(.data$interInt, 0.975 * max(.data$interInt, na.rm=T), max(.data$interInt, na.rm=T))) %>%
    filter(.data$interHits == max(.data$interHits)) %>%
    pull(.data$index)
  print(tuned.parse)
  suppressWarnings(plot(tuned.plot))
  tuned[[length(tuned) + 1]] <- tuned.parse
  tuned[[length(tuned) + 1]] <- bestModelIndex
  bestModel <- tuned[[bestModelIndex]]
  return(list(
    "CSMs" = bestModel$CSMs,
    "URPs" = bestModel$URPs,
    "CSM.thresh" = findSeparateThresholdsModelled(bestModel$CSMs, targetER = targetER, scalingFactor = scalingFactor),
    "URP.thresh" = bestModel$thresh,
    "model.params" = list(
      "kernel" = bestModel$kernel,
      "cost" = bestModel$cost,
      "gamma" = bestModel$gamma,
      "sd.thresh" = bestModel$sd.thresh,
      "features" = bestModel$params)
  ))
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
#' @param gamma_values Numeric vector of gamma values used for hyperparameter tuning of the SVM model
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
                    gamma_values = c(0.01, 0.05, 0.1, 0.5)) {
  param_grid <- expand.grid(cost=cost_values, gamma=gamma_values)
  param_grid$kernel = "radial"
  param_grid = bind_rows(param_grid, data.frame(cost=cost_values, gamma=23, kernel="linear"))
  tuned <- purrr::pmap(param_grid, function(cost, gamma, kernel) {
    tuneSVM.helper(datTab=datTab,
                   targetER=targetER,
                   params=params,
                   scoreName=scoreName,
                   scalingFactor = scalingFactor,
                   sampleNo = sampleNo,
                   cost, gamma, kernel)
  })
  return(tuned)
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
                           cost, gamma, kernel) {
  datTab.csm <- buildSVM(datTab=datTab,
                         targetER=targetER,
                         params=params,
                         scoreName=scoreName,
                         sampleNo = sampleNo,
                         showTab = F,
                         cost=cost, gamma=gamma, kernel=kernel)
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

  list("CSMs" = datTab.csm,
       "URPs" = datTab.urp,
       "thresh" = datTab.urp.thresh,
       "intraHits" = intraHits,
       "interHits" = interHits,
       "errorTable" = errorTable,
       "interInt" = inter.integral,
       "cost" = cost,
       "gamma" = gamma,
       "kernel" = kernel,
       "sd.thresh" = min(datTab.csm$Score.Diff),
       "params" = params)
}

#' Basic function to build a new SVM classifier.  Doesn't do any feature selection or
#' hyperparamter tuning. Build two separate SVM models on different subsets of the data
#' and averages the results.
#'
#' @param datTab Parsed CLMS search results.
#' @param params Character vector specifying names of the features in `datTab` used to train model.
#' @param scoreName Name for the new scoring function.
#' @param sampleNo Size of the training dataset (integer).
#' @param showTab print classificaiton table?
#' @param ... paramters passed to `e1071:svm()` function
#' @seealso [trainCrosslinkScore()], [tuneSVM.helper()], [tuneSVM()]
#' @return A data frame, one column larger than the input containing the new score.
#' @export
buildSVM <- function(datTab,
                     params=params.best,
                     scoreName="SVM.score",
                     sampleNo = 20000,
                     showTab = F,
                     ...) {
  datTab$massError <- abs(datTab$ppm - mean(datTab$ppm))
  num.rows <- nrow(datTab)
  if ((num.rows) < 40000L) {
    sampleNo <- num.rows %/% 2
  } else {
    sampleNo <- 20000L
  }
  ind.1 <- sample(1:num.rows, sampleNo)
  ind.2 <- sample(c(1:num.rows)[-1*ind.1], sampleNo)
  train.1 <- datTab[ind.1,]
  train.2 <- datTab[ind.2,]
  test.1 <- datTab[-1 * ind.1,]
  test.2 <- datTab[-1 * ind.2,]
  wghts.1 <- numeric(0)
  wghts.2 <- numeric(0)
  wghts.1["Target"] <- table(test.1$Decoy2)["Decoy"] / sum(table(test.1$Decoy2),na.rm=T)
  wghts.1["Decoy"] <- table(test.1$Decoy2)["Target"] / sum(table(test.1$Decoy2),na.rm=T)
  wghts.2["Target"] <- table(test.2$Decoy2)["Decoy"] / sum(table(test.2$Decoy2),na.rm=T)
  wghts.2["Decoy"] <- table(test.2$Decoy2)["Target"] / sum(table(test.2$Decoy2),na.rm=T)
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
  if (stats::cor(datTab$Score.Diff, datTab$score.1) < 0) {datTab$score.1 <- -1 * datTab$score.1}
  if (stats::cor(datTab$Score.Diff, datTab$score.2) < 0) {datTab$score.2 <- -1 * datTab$score.2}
  datTab[ind.1, "score.1"] <- NA
  datTab[ind.2, "score.2"] <- NA
  datTab[[scoreName]] <- purrr::map2_dbl(datTab$score.1, datTab$score.2, function(x, y) mean(c(x, y), na.rm=T))
  if (showTab) {
    tab <- table(datTab$Decoy2, datTab[[scoreName]] > 0)
    print(tab)
    print(paste("specificity:", round(tab[1]/(tab[1]+tab[3]),2)))
  }
  return(datTab)
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
      nrow(plausibleHits) <= 1 ~ list(params.best.nop),
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
  requireNamespace(c("furrr","future"), quietly = TRUE)
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
    datTab <- classifySeparateThresholds(datTab, preFilter.thresh, classifier=.data$Score.Diff)
  }
  buildSVM(datTab, params, scoreName, sampleNo, ...)
}

