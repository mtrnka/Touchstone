library(tidyverse)
library(devtools)
load_all()

tuned.plot <- function(tuned) {
  tuned %>%
  map_dfr(function(x) {
    df <- x$errorTable
    df <- df %>% mutate(
      kernel=x$kernel,
      cost=x$cost,
      gamma=x$gamma)
  }) %>%
  filter(fdr.inter <= 0.15) %>%
  filter(inter >= 0.05 * max(inter)) %>%
  ggplot(aes(x=fdr.inter, y=inter, col=as.factor(cost))) +
  geom_line(linewidth=1.2) +
  theme_bw() +
#  xlim(0,0.15) +
  geom_vline(xintercept = 0.01, color="red") +
  scale_color_viridis_d(option="C") +
  facet_grid(rows=vars(gamma), scales="free_y")
}


setDecoyScalingFactor(10)


test.csm.0 <- readProspectorXLOutput("inst/extdata/rRibo_DSSO_sthcd_scOut.txt", minIons = 0)

test.tune.0 <- trainClassifer_new(test.csm.0,
                                  cost_values = c(0.1,1,10),
                                  gamma_values = c(0.001, 0.01, 0.1,1,10),
                                  preFilter = F)







bind_rows(
  mutate(train.test.sd10[[12]]$errorTable, exp="sd10"),
  mutate(train.test.sd5[[12]]$errorTable, exp="sd5"),
  mutate(train.test[[8]]$errorTable, exp="sd0")
  ) %>%
  filter(fdr.inter <= 0.1) %>%
  filter(inter >= 0.1 * max(inter)) %>%
  ggplot(aes(x=fdr.inter, y=inter, col=exp)) +
  geom_line(linewidth=1.2) +
  theme_bw() +
  xlim(0,0.1) +
  geom_vline(xintercept = 0.01, color="red") +
  scale_color_viridis_d(option="C")

t0 <- trainClassifer_new(filter(test.csm.0, Score.Diff >= 0), preFilter = F)
t5 <- trainClassifer_new(filter(test.csm.0, Score.Diff >= 5), preFilter = F)
t10 <- trainClassifer_new(filter(test.csm.0, Score.Diff >= 10), preFilter = F)

bind_rows(
  mutate(t0[[15]]$errorTable, exp="sd0"),
  mutate(t5[[11]]$errorTable, exp="sd5"),
  mutate(t10[[16]]$errorTable, exp="sd10")
) %>%
  filter(fdr.inter <= 0.1) %>%
  filter(inter >= 0.1 * max(inter)) %>%
  ggplot(aes(x=fdr.inter, y=inter, col=exp)) +
  geom_line(linewidth=1.2) +
  theme_bw() +
  xlim(0,0.1) +
  geom_vline(xintercept = 0.01, color="red") +
  scale_color_viridis_d(option="C")








projDir <- "~/Projects-Mine/HEKcompo/"
setDecoyScalingFactor(1)
h1 <- readProspectorXLOutput(file.path(projDir, "hek_hcd", "sc", "HEK_hcd_scout.txt"), minPepLen = 5) %>%
  mutate(acqMethod = "ms2.sthcd")
h1.tune <- trainClassifer_new(h1)

g1 <- readProspectorXLOutput(file.path(projDir, "hek_hcd2", "sc", "HEK_hcd2_scout.txt"), minPepLen = 5) %>%
  mutate(acqMethod = "ms2.sthcd")
g1.tune <- trainClassifer_new(g1)

# g1.tune.sd0 <- tuneSVM(g1)
# g1.tune.sd5 <- tuneSVM(filter(g1, Score.Diff >= 5))
# g1.tune.sd10 <- tuneSVM(filter(g1, Score.Diff >= 10))
# g1.tune.sd20 <- tuneSVM(filter(g1, Score.Diff >= 20))
# g1.tune.sd25 <- tuneSVM(filter(g1, Score.Diff >= 25))
#
# bind_rows(
#   mutate(g1.tune[[6]]$errorTable, sd=15),
#   mutate(g1.tune.sd5[[6]]$errorTable, sd=5),
#   mutate(g1.tune.sd10[[6]]$errorTable, sd=10),
#   mutate(g1.tune.sd20[[6]]$errorTable, sd=20),
#   mutate(g1.tune.sd25[[6]]$errorTable, sd=25)
# ) %>%
#   filter(fdr.inter <= 0.15) %>%
#   filter(inter >= 0.05 * max(inter)) %>%
#   ggplot(aes(x=fdr.inter, y=inter, col=as.factor(sd))) +
#   geom_line(linewidth=1.2) +
#   theme_bw() +
#   #  xlim(0,0.15) +
#   geom_vline(xintercept = 0.01, color="red") +
#   scale_color_viridis_d(option="C")
#
# g1.tune.sd20.thresh <- findSeparateThresholdsModelled(g1.tune.sd20[[11]]$URPs)

hcd1.comb <- bind_rows(classifyDataset(h1.tune[[9]]$URPs, h1.tune[[9]]$thresh),
                       classifyDataset(g1.tune[[9]]$URPs, g1.tune[[9]]$thresh)) %>%
  bestResPair()

classifyDataset(g1.tune[[9]]$URPs, g1.tune[[9]]$thresh) %>%
  count(Decoy, xlinkClass) %>% pivot_wider(names_from=xlinkClass, values_from=n)

classifyDataset(h1.tune[[9]]$URPs, h1.tune[[9]]$thresh) %>%
  count(Decoy, xlinkClass) %>% pivot_wider(names_from=xlinkClass, values_from=n)

fdrPlots(hcd1.comb, map2(h1.tune[[9]]$thresh, g1.tune[[9]]$thresh, function(x, y) mean(c(x,y))))
hcd1.comb %>% count(Decoy, xlinkClass) %>% pivot_wider(names_from=xlinkClass, values_from=n)

h1.prot <- bestProtPair(h1.tune[[9]]$CSMs)
g1.prot <- bestProtPair(g1.tune[[9]]$CSMs)

h1.prot.thresh <- findSeparateThresholdsModelled(h1.prot, targetER = 0.01)
g1.prot.thresh <- findSeparateThresholdsModelled(g1.prot, targetER = 0.01)

hcd1.comp.prot <- bind_rows(classifyDataset(h1.prot, h1.prot.thresh),
                            classifyDataset(g1.prot, g1.prot.thresh)) %>%
  bestProtPair()

fdrPlots(hcd1.comp.prot, map2(h1.prot.thresh, g1.prot.thresh, function(x, y) mean(c(x,y))))
hcd1.comp.prot %>% count(Decoy, xlinkClass) %>% pivot_wider(names_from=xlinkClass, values_from=n)
#hcd1.comp.prot %>% formatXLTable() %>% filter(xlinkClass=="interProtein") %>% View


library(STRINGdb)
string_db <- STRINGdb$new(version = "12.0", network_type="full", link_data="combined_only",
                          species = 9606, score_threshold = 0, input_directory = "")

p.list.1 <- hcd1.comp.prot %>%
  filter(xlinkClass=="interProtein") %>%
  removeDecoys() %>%
  pull(Acc.1) %>%
  as.character()
p.list.2 <- hcd1.comp.prot %>%
  filter(xlinkClass=="interProtein") %>%
  removeDecoys() %>%
  pull(Acc.2) %>%
  as.character()
p.list <- unique(c(p.list.1, p.list.2))

id.map <- string_db$map(data.frame(acc = p.list), "acc", removeUnmappedRows = F)

getStringDBscore <- function(Acc.1, Acc.2) {
  String.1 <- id.map %>% filter(acc == Acc.1) %>% pull(STRING_id)
  String.2 <- id.map %>% filter(acc == Acc.2) %>% pull(STRING_id)
  ppi <- string_db$get_interactions(c(String.1, String.2))
  if (length(ppi$combined_score)==0) {
    return(NA)
  } else {
    return(ppi$combined_score[[1]])
  }
}

hcd1.comp.prot <- hcd1.comp.prot %>%
  mutate(string.score = map2_dbl(Acc.1, Acc.2, getStringDBscore),
         string.ppi = case_when(is.na(string.score) ~ F,
                                string.score < 900 ~ F,
                                string.score >= 900 ~ T))

hcd1.comp.prot %>%
  filter(xlinkClass=="interProtein") %>%
  removeDecoys() %>%
  count(string.ppi)

ko_ppi.ppi %>%
  filter(xlinkClass=="interProtein") %>%
  removeDecoys() %>%
  count(string.ppi)


hcd1.comp.prot %>%
  filter(xlinkClass=="interProtein") %>%
  formatXLTable(extraCols = c("string.score","string.ppi")) %>%
  View

ko_ppi.ppi %>%
  filter(xlinkClass=="interProtein") %>%
  formatXLTable(extraCols = c("string.score","string.ppi")) %>%
  View


#save(hcd1.comp.prot, file=file.path(projDir, "hcd1.comp.prot.RData"))
#load(file.path(projDir, "hcd1.comp.prot.RData"))


############
train.test.too <- tuneSVM(test.csm.0)
train.test.three_proj <- tuneSVM(test.csm.0)
train.test.4_proj_revert_wghts <- tuneSVM(test.csm.0)
train.test.5_proj_fixed <- tuneSVM(test.csm.0)

train.test.6 <- trainClassifer_new(test.csm.0)

##################################################

reDir <- "~/Projects-Collaboration/Rappsilber/241101_re/"
scOut <- list.files(reDir, pattern = "ecoli[0-9]_6.5_scout\\.txt")
setDecoyScalingFactor(1)
ecAll <- map_dfr(scOut, function(x) readProspectorXLOutput(file.path(reDir, x), minPepLen = 5))
ecAll <- ecAll %>%
  mutate(entrapment = (stringr::str_detect(Species.1, "HUMAN") | stringr::str_detect(Species.2, "HUMAN")))

ecAll.inter <- ecAll %>%
  filter(xlinkClass == "interProtein") %>%
  calculatePairs()

params.best.noclass <- params.best[params.best != "xlinkClass"]

ecAll.tune <- trainClassifer_new(ecAll, targetER = 0.02, params = params.best)

ecAll.tuned <- ecAll.tune[[8]]
ecAll.tuned$URPs %>% pull(Score.Diff) %>% min
ecAll.tuned$thresh

ecAll.pp <- bestProtPair(ecAll.tuned$CSMs)
ecAll.pp.thresh <- findSeparateThresholdsModelled(ecAll.pp, targetER = 0.02)
fdrPlots(ecAll.pp, ecAll.pp.thresh)
ecAll.pp %>%
  classifyDataset(ecAll.pp.thresh) %>%
  group_by(xlinkClass, Decoy, entrapment) %>%
  count() %>%
  pivot_wider(names_from = Decoy, values_from = n)

p.list.1 <- ecAll.pp %>%
  filter(!entrapment, xlinkClass=="interProtein") %>%
  removeDecoys() %>%
  pull(Acc.1) %>%
  as.character()
p.list.2 <- ecAll.pp %>%
  filter(!entrapment, xlinkClass=="interProtein") %>%
  removeDecoys() %>%
  pull(Acc.2) %>%
  as.character()
p.list <- unique(c(p.list.1, p.list.2))

string_db <- STRINGdb$new(version = "12.0", network_type="full", link_data="combined_only",
                          species = 511145, score_threshold = 0, input_directory = "")

id.map <- string_db$map(data.frame(acc = p.list), "acc", removeUnmappedRows = F)

ppts.ppi <- ecAll.pp %>%
  classifyDataset(ecAll.pp.thresh) %>%
  mutate(MSMS.Info = str_extract(MSMS.Info, "[[0-9]]+$")) %>%
  filter(!entrapment, xlinkClass == "interProtein") %>%
  removeDecoys() %>%
  mutate(string.score = map2_dbl(Acc.1, Acc.2, getStringDBscore),
         string.ppi = case_when(is.na(string.score) ~ F,
                                string.score < 400 ~ F,
                                string.score >= 400 ~ T)
  )







getStringScores <- function(datTab, ncbiTaxonomyCode = NULL) {
  primarySpecies <- datTab %>%
    removeDecoys() %>%
    filter(Score.Diff > 5) %>%
    count(Species.1) %>%
    arrange(desc(n)) %>%
    slice(1) %>%
    pull(Species.1)
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
    print(str_c("detected organism: ", primarySpecies, "\tncbi code:", ncbiTaxonomyCode))
  },
  error = function(cond) {
    message("Unknown species, please provide the ncbi taxonomy identifier")
    message("Original error message:")
    message(conditionMessage(cond))
    NA
  })
  string_db <- STRINGdb$new(version = "12.0", network_type="full", link_data="combined_only",
                            species = ncbiTaxonomyCode, score_threshold = 0, input_directory = "")
  datTab.inter <- datTab %>%
    filter(xlinkClass == "interProtein")
  datTab.intra <- datTab %>%
    filter(xlinkClass == "intraProtein")
  p.list.1 <- datTab.inter %>%
    removeDecoys() %>%
    pull(Acc.1) %>%
    as.character()
  p.list.2 <- datTab.inter %>%
    removeDecoys() %>%
    pull(Acc.2) %>%
    as.character()
  p.list <- unique(c(p.list.1, p.list.2))
  id.map <- string_db$map(data.frame(acc = p.list), "acc", removeUnmappedRows = F)
  datTab.intra$string.score <- NA
  datTab.inter <- datTab.inter %>%
    mutate(string.score = map2_dbl(Acc.1, Acc.2, function(x, y) {
      String.1 <- id.map %>% filter(acc == x) %>% pull(STRING_id)
      String.2 <- id.map %>% filter(acc == y) %>% pull(STRING_id)
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

getStringDBscore <- function(Acc.1, Acc.2) {
  String.1 <- id.map %>% filter(acc == Acc.1) %>% pull(STRING_id)
  String.2 <- id.map %>% filter(acc == Acc.2) %>% pull(STRING_id)
  ppi <- string_db$get_interactions(c(String.1, String.2))
  if (length(ppi$combined_score)==0) {
    return(NA)
  } else {
    return(ppi$combined_score[[1]])
  }
}

names(test)
summary(test$string.score)

