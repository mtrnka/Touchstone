setDecoyScalingFactor(10)
ribo.xl <- readProspectorXLOutput("inst/extdata/rRibo_DSSO_sthcd_scOut.txt")
ribo.tune <- trainCrosslinkScore(ribo.xl, params = params.best)
ribo.tune2 <- trainCrosslinkScore(ribo.xl, params = c("Score.Diff", "percMatched", "massError",
                                                    "z", "numURP", "numCSM", "xlinkClass",
                                                    "Perc.Bond.Cleavage.1", "Perc.Bond.Cleavage.2"))
ribo.lin <- buildSVM(ribo.xl, kernel = "linear", cost=1, gamma=1)

ribo.csm <- ribo.tune$CSMs
ribo.csm.1 <- ribo.tune$CSM.thresh

ribo2.csm <- ribo.tune2$CSMs
ribo2.csm.1 <- ribo.tune2$CSM.thresh

ribo.csm.lin <- ribo.lin
ribo.csm.lin.1 <- findSeparateThresholdsModelled(ribo.csm.lin, targetER = 0.01)

ribo.csm <- ribo.csm %>%
  processModuleFile("inst/extdata/rRibo_newMod_uniprot.txt")
ribo.urp <- ribo.csm %>%
  bestResPair()
ribo.urp.1 <- findSeparateThresholdsModelled(ribo.urp, targetER = 0.01)
ribo.ppi <- bestProtPair(ribo.csm)
ribo.ppi.1 <- findSeparateThresholdsModelled(ribo.ppi, targetER = 0.01)

ribo2.urp <- ribo2.csm %>%
  bestResPair()
ribo2.urp.1 <- findSeparateThresholdsModelled(ribo2.urp, targetER = 0.01)

ribo.urp.lin <- bestResPair(ribo.csm.lin)
ribo.urp.lin.1 <- findSeparateThresholdsModelled(ribo.urp.lin, targetER = 0.01)

fdrPlots(ribo.csm, ribo.csm.1)
calculateFDR(ribo.csm, ribo.csm.1)

ribo.csm %>%
  classifyDataset(ribo.csm.1) %>%
  countDecoys()

fdrPlots(ribo.urp, ribo.urp.1)
calculateFDR(ribo.urp, ribo.urp.1)

fdrPlots(ribo2.urp, ribo2.urp.1)
calculateFDR(ribo2.urp, ribo2.urp.1)

fdrPlots(ribo.urp.lin, ribo.urp.lin.1)
calculateFDR(ribo.urp.lin, ribo.urp.lin.1)

ribo.urp %>%
  classifyDataset(ribo.urp.1) %>%
  countDecoys()

ribo2.urp %>%
  classifyDataset(ribo2.urp.1) %>%
  countDecoys()

ribo.urp.lin %>%
  classifyDataset(ribo.urp.lin.1) %>%
  countDecoys()

ribo.urp %>%
  classifyDataset(ribo.urp.1) %>%
  distancePlot2(threshold = 35)

fdrPlots(ribo.ppi, ribo.ppi.1)
calculateFDR(ribo.ppi, ribo.ppi.1)

ribo.ppi %>%
  classifyDataset(ribo.ppi.1) %>%
  ggplot(aes(wtCSM)) +
  geom_histogram(color="white") +
  facet_grid(rows = vars(Decoy2), scales="free_y")

ribo.csm %>%
  classifyDataset(ribo.ppi.1) %>%
  calculatePairs() %>%
  bestProtPair() %>%
  ggplot(aes(wtCSM)) +
  geom_histogram(color="white") +
  facet_grid(rows = vars(Decoy2), scales="free_y")
