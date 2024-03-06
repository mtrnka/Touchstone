#' Read Search Compare Output from Protein Prospector
#'
#' @param inputFile A tab delimited search compare file generated using tstoneMS2.json param file
#' @param minPepLen minimum length of each peptide in amino acids.
#' @param minPepScore minimum prospector score of each peptide.
#' @param minScoreDiff minimum Score Difference for the CSM.
#' @param minIons mimimum number of product ions matched for each peptide.
#'
#' @return A data frame
#' @export
#'
#' @examples
readProspectorXLOutput <- function(inputFile, minPepLen = 3, minPepScore = 0, minScoreDiff = 0, minIons = 3){

  datTab <- read_tsv(inputFile, guess_max = 10000)
  header <- names(datTab) %>%
    str_replace("_[[0-9]]$", "") %>%
    str_replace_all("[[:space:]]",".") %>%
    str_replace_all("#", "Num") %>%
    str_replace_all("%", "Perc")
  acc_pos <- str_which(header, "Acc")
  header[acc_pos] <- c("Acc.1", "Acc.2")
  spec_pos <- str_which(header, "Species")
  header[spec_pos] <- c("Species.1", "Species.2")
  prot_pos <- str_which(header, "Protein.Name")
  header[prot_pos] <- c("Protein.1", "Protein.2")
  if (sum(str_detect(header, "Protein.MW")) == 2) {
    mw_pos <- str_which(header, "Protein.MW")
    header[mw_pos] <- c("MW.1", "MW.2")
  }
  if (sum(str_detect(header, "Protein.Length")) == 2) {
    len_pos <- str_which(header, "Protein.Length")
    header[len_pos] <- c("Protein.len.1", "Protein.len.2")
  }
  if (sum(str_detect(header, "MSMS.Ions")) == 2) {
    len_pos <- str_which(header, "MSMS.Ions")
    header[len_pos] <- c("MSMS.Ions.1", "MSMS.Ions.2")
  }
  if (sum(str_detect(header, "MSMS.M/Zs")) == 2) {
    len_pos <- str_which(header, "MSMS.M/Zs")
    header[len_pos] <- c("MSMS.MZs.1", "MSMS.MZs.2")
  }
  if (sum(str_detect(header, "MSMS.Intensities")) == 2) {
    len_pos <- str_which(header, "MSMS.Intensities")
    header[len_pos] <- c("MSMS.Intensities.1", "MSMS.Intensities.2")
  }
  if (sum(str_detect(header, "MSMS.Errors")) == 2) {
    len_pos <- str_which(header, "MSMS.Errors")
    header[len_pos] <- c("MSMS.Errors.1", "MSMS.Errors.2")
  }
  if (sum(str_detect(header, "Perc.Bond.Cleavage")) == 2) {
    len_pos <- str_which(header, "Perc.Bond.Cleavage")
    header[len_pos] <- c("Perc.Bond.Cleavage.1", "Perc.Bond.Cleavage.2")
  }
  names(datTab) <- header
  if (!"Spectrum" %in% names(datTab)) {
    datTab$Spectrum <- 1
  }
  if (!"distance" %in% names(datTab)) {
    datTab$distance <- NA_real_
  }
  datTab <- datTab %>% mutate(
    Acc.1 = as.character(Acc.1),
    Acc.2 = as.character(Acc.2)
  )
  datTab <- datTab %>% mutate(Protein.1 = if_else(is.na(Protein.1), Acc.1, Protein.1),
                              Protein.2 = if_else(is.na(Protein.2), Acc.2, Protein.2))
  datTab <- calculateDecoys(datTab)
  datTab <- assignXLinkClass(datTab)
  datTab <- calculatePercentMatched(datTab)
  datTab <- calculatePeptideLengths(datTab)
  datTab <- lengthFilter(datTab, minLen = minPepLen, maxLen = 40)
  datTab <- scoreFilter(datTab, minScore = minPepScore)
  datTab <- datTab %>% filter(Score.Diff >= minScoreDiff)
  datTab <- calculatePairs(datTab)
  if ("numProdIons.1" %in% names(datTab) & "numProdIons.1" %in% names(datTab)) {
    datTab <- productIonFilter(datTab, minProducts.1 = minIons, minProducts.2 = minIons) }
  datTab <- datTab %>%
    filter(xlinkedResPair != "0.decoy::0.decoy")
  return(datTab)
}

#' determine whether a CSM is a decoy or target
#'
#' @param datTab parsed CLMS search results
#'
#' @return
#' @export
#'
#' @examples
calculateDecoys <- function(datTab) {
  datTab$Decoy <- "Target"
  decoyReg <- "(^r[[:digit:]]_|^dec|^DECOY)"
  datTab[grepl(decoyReg, datTab$Acc.1) |
           grepl(decoyReg, datTab$Acc.2),
         "Decoy"] <- "Decoy"
  datTab[grepl(decoyReg, datTab$Acc.1) &
           grepl(decoyReg, datTab$Acc.2),
         "Decoy"] <- "DoubleDecoy"
  datTab$Decoy2 <- "Target"
  datTab[grepl("Decoy", datTab$Decoy), "Decoy2"] <- "Decoy"
  datTab$Decoy2 <- factor(datTab$Decoy2, levels=c("Decoy","Target"))
  datTab$Decoy <- factor(datTab$Decoy, levels=c("DoubleDecoy","Decoy","Target"))
  return(datTab)
}

#' determine residue, peptide,and protein pairs
#'
#' @param datTab
#'
#' @return
#' @export
#'
#' @examples
calculatePairs <- function(datTab){
  datTab <- datTab %>%
    mutate(Acc.1 = as.character(Acc.1),
           Acc.2 = as.character(Acc.2))
  datTab$xlinkedProtPair <- ifelse(datTab$Decoy2 == "Decoy",
                                   ifelse(datTab$Protein.1 <= datTab$Protein.2,
                                          paste("decoy", datTab$Protein.1, datTab$Protein.2, sep="::"),
                                          paste("decoy", datTab$Protein.2, datTab$Protein.1, sep="::")),
                                   ifelse(datTab$Acc.1 <= datTab$Acc.2,
                                          paste(datTab$Acc.1, datTab$Acc.2, sep="::"),
                                          paste(datTab$Acc.2, datTab$Acc.1, sep="::")))
  accs <- unique(c(datTab$Acc.1, datTab$Acc.2)) %>% str_sort()
  datTab <- datTab %>% mutate(Acc.1 = factor(Acc.1, levels = accs),
                              Acc.2 = factor(Acc.2, levels = accs))
  datTab$Res.1 <- paste(datTab$XLink.AA.1, datTab$Acc.1, sep=".")
  datTab$Res.2 <- paste(datTab$XLink.AA.2, datTab$Acc.2, sep=".")
  datTab$xlinkedResPair <- ifelse(datTab$Res.1 <= datTab$Res.2,
                                  paste(datTab$Res.1, datTab$Res.2, sep="::"),
                                  paste(datTab$Res.2, datTab$Res.1, sep="::"))
  datTab$xlinkedPepPair <- ifelse(datTab$DB.Peptide.1 <= datTab$DB.Peptide.2,
                                  paste(datTab$DB.Peptide.1, datTab$DB.Peptide.2, sep="::"),
                                  paste(datTab$DB.Peptide.2, datTab$DB.Peptide.1, sep="::"))
  datTab$xlinkedResPair <- as.factor(datTab$xlinkedResPair)
  datTab$xlinkedProtPair <- as.factor(datTab$xlinkedProtPair)
  datTab$xlinkedPepPair <- as.factor(datTab$xlinkedPepPair)
  datTab <- datTab %>%
    add_count(xlinkedResPair, name="numCSM")
  uniqueProtCount <- datTab %>%
    select(xlinkedProtPair, xlinkedResPair) %>%
    group_by(xlinkedProtPair, xlinkedResPair) %>%
    summarize(.groups="drop") %>%
    group_by(xlinkedProtPair) %>%
    add_count(name = "numURP") %>%
    ungroup()
  datTab <- left_join(select(datTab, -any_of("numURP")), uniqueProtCount, by=c("xlinkedProtPair","xlinkedResPair"))
  if ("Module.1" %in% names(datTab) & "Module.2" %in% names(datTab)) {
    datTab <- datTab %>%
      mutate(Module.1 = as.character(Module.1),
             Module.2 = as.character(Module.2))
    mods <- unique(c(datTab$Module.1, datTab$Module.2)) %>% str_sort()
    datTab <- datTab %>%
      mutate(xlinkedModulPair = ifelse(Module.1 <= Module.2,
                                       str_c(Module.1, Module.2, sep="::"),
                                       str_c(Module.2, Module.1, sep="::"))
      ) %>%
      mutate(xlinkedModulPair = as_factor(xlinkedModulPair),
             Module.1 = factor(Module.1, levels = mods),
             Module.2 = factor(Module.2, levels = mods)) #%>%
    #         add_count(xlinkedModulPair, name="numMPSM")
  }
  if ("MSMS.Ions.1" %in% names(datTab) & "MSMS.Ions.2" %in% names(datTab)) {
    datTab <- calculateProductIons(datTab)
  }
  return(datTab)
}

#' count number of CSMs per URP and number of URPs per PP
#'
#' @param datTab
#'
#' @return
#' @export
#'
#' @examples
addXlinkCount <- function(datTab) {
  datTab <- datTab %>%
    add_count(xlinkedResPair, name="numCSM")
  uniqueProtCount <- datTab %>%
    select(xlinkedProtPair, xlinkedResPair) %>%
    group_by(xlinkedProtPair, xlinkedResPair) %>%
    summarize(.groups="drop") %>%
    group_by(xlinkedProtPair) %>%
    count(name = "numURP") %>%
    ungroup()
  datTab <- left_join(select(datTab, -any_of("numURP")), uniqueProtCount, by="xlinkedProtPair")
}

#' is the crosslink interprotein (hetromeric) or intraprotein (homomeric)
#'
#' @param datTab
#'
#' @return
#' @export
#'
#' @examples
assignXLinkClass <- function(datTab) {
  intraProteinLinks <- datTab$Protein.1 == datTab$Protein.2
  interProteinLinks <- datTab$Protein.1 != datTab$Protein.2
  if (sum(c("Start.1", "Start.2", "End.1", "End.2") %in% names(datTab)==4)) {
    s1 <- datTab$Start.1
    s2 <- datTab$Start.2
    e1 <- datTab$End.1
    e2 <- datTab$End.2
    condition1 <- (s1 >= s2) & (s1 <= e2)
    condition2 <- (e1 >= s2) & (e1 <= e2)
    interHomomerLinks <- condition1 | condition2
    datTab$xlinkClass <- ifelse(intraProteinLinks,
                                ifelse(interHomomerLinks,
                                       "interProtein, homomeric",
                                       "intraProtein"),
                                "interProtein, heteromeric")
  } else {
    datTab$xlinkClass <- ifelse(intraProteinLinks, "intraProtein", "interProtein")
  }
  return(datTab)
}

#' If MS-Product reporting of %TIC matched is present use that.  Otherwise approximate with number of matched ions.
#'
#' @param datTab
#'
#' @return
#' @export
#'
#' @examples
calculatePercentMatched <- function(datTab) {
  if ("Match.Int" %in% names(datTab)) {
    datTab$percMatched <- datTab$Match.Int
  } else if ("Num.Pks" %in% names(datTab) & "Num.Unmat" %in% names(datTab)) {
    num.Matched <- datTab$Num.Pks - datTab$Num.Unmat
    datTab$percMatch <- num.Matched / datTab$Num.Pks
  }
  return(datTab)
}

#' calculate length of each peptide in amino acid residues.
#'
#' @param datTab
#'
#' @return
#' @export
#'
#' @examples
calculatePeptideLengths <- function(datTab) {
  datTab$Len.Pep.1 <- nchar(datTab$DB.Peptide.1)
  datTab$Len.Pep.2 <- nchar(datTab$DB.Peptide.2)
  return(datTab)
}

#' are diagnostic MH*/# ions from cleavable crosslinking reagents detected?
#'
#' @param msms.ions columns containing list of detected ions, must be present in Search Compare output.
#' @param msms.int  columns with corresponding intensity vaues fro detected ions.
#'
#' @return
#' @export
#'
#' @examples
getDiagnosticMatches <- function(msms.ions, msms.int) {
  ions <- unlist(str_split(msms.ions, ";"))
  ints <- as.numeric(unlist(str_split(msms.int, ";")))
  total.int <- sum(ints, na.rm = T)
  #   max.int <- max(ints, na.rm=T)
  mh_a <- str_which(ions, "MH\\#")
  ion_a <- ions[mh_a]
  int_a <- sum(ints[mh_a]/total.int)
  mh_s <- str_which(ions, "MH\\*")
  ion_s <- ions[mh_s]
  int_s <- sum(ints[mh_s]/total.int)
  if (length(mh_a) == 0) ion_a <- NA
  if (length(mh_a) == 0) int_a <- 0
  if (length(mh_s) == 0) ion_s <- NA
  if (length(mh_a) == 0) int_a <- 0
  diagnostics <- list(ion_a, int_a, ion_s, int_s)
}

#' determined whether diagnostic ion pair from cleavable crosslinker was deteceted
#'
#' @param datTab
#' @param int.thresh minimum intensity for diagnostic ion.
#'
#' @return
#' @export
#'
#' @examples
calculateDiagnosticPairs <- function(datTab, int.thresh = 0.05) {
  datTab <- datTab %>%
    mutate(diagnostics.1 = map2(MSMS.Ions.1, MSMS.Intensities.1, getDiagnosticMatches),
           diagnostics.2 = map2(MSMS.Ions.2, MSMS.Intensities.2, getDiagnosticMatches),
           aaPair = map2_lgl(diagnostics.1, diagnostics.2, function(x, y) {
             if_else(x[[2]] > int.thresh & y[[2]] > int.thresh, T, F)}),
           ttPair = map2_lgl(diagnostics.1, diagnostics.2, function(x, y) {
             if_else(x[[4]] > int.thresh & y[[4]] > int.thresh, T, F)}),
           atPair = map2_lgl(diagnostics.1, diagnostics.2, function(x, y) {
             if_else(x[[2]] > int.thresh & y[[4]] > int.thresh, T, F)}),
           taPair = map2_lgl(diagnostics.1, diagnostics.2, function(x, y) {
             if_else(x[[4]] > int.thresh & y[[2]] > int.thresh, T, F)})
    )
  datTab <- datTab %>%
    mutate(diagnosticPair = pmap_lgl(list(datTab$aaPair,datTab$ttPair,datTab$taPair,datTab$atPair),
                                     function(a,b,c,d) sum(a,b,c,d, na.rm=T) >= 1)
    )
}

#' function for sequential scan data, e.g. sequential MS2-CID and MS2-ETD determined triggered on the same precursor ion.
#'
#' @param datTab.etd parsed search table that is not expected to contain diagnostic ions (does not necessarily need to be ETD)
#' @param datTab.cid parsed search table expected to contain diagnostic ions (typical MS2-CID)
#' @param masterScanFile master scan file that charts the dependencies between different scans.
#'
#' @return
#' @export
#'
#' @examples
calculateCrossDiagnosticPair <- function(datTab.etd, datTab.cid, masterScanFile) {
  datTab.etd <- datTab.etd %>%
    mutate(Fraction = str_replace_all(Fraction, "_FTMSms2[[a-z]]+$", ""))
  datTab.cid <- datTab.cid %>%
    mutate(Fraction = str_replace_all(Fraction, "_FTMSms2[[a-z]]+$", ""))
  masterScanFile <- masterScanFile %>%
    mutate(fraction = str_replace_all(fraction, "_FTMSms2[[a-z]]+$", ""))
  datTab.etd <- left_join(datTab.etd, select(masterScanFile, fraction, ms2cidScanNo, ms2etdScanNo),
                          by=c("Fraction"="fraction", "MSMS.Info"="ms2etdScanNo"))
  datTab.etd <- left_join(datTab.etd, select(datTab.cid, Fraction, MSMS.Info, diagnosticPair),
                          by=c("Fraction"="Fraction", "ms2cidScanNo"="MSMS.Info"),
                          suffix=c("",".cross"))
  datTab.etd <- datTab.etd %>%
    mutate(diagnosticPair.cross = if_else(is.na(diagnosticPair.cross), F, diagnosticPair.cross))
  return(datTab.etd)
}

#' Analog of getDiagnostic Matches for non-cleavable crosslinkers. Looks for P, PL, PLK type ions.
#'
#' @param msms.ions
#' @param msms.int
#'
#' @return
#' @export
#'
#' @examples
getDiagnosticMatchesNonCleavable <- function(msms.ions, msms.int) {
  ions <- unlist(str_split(msms.ions, ";"))
  ints <- as.numeric(unlist(str_split(msms.int, ";")))
  total.int <- sum(ints, na.rm = T)
  mh_p <- str_which(ions, "P[^\\+\\-]")
  ion_p <- ions[mh_p]
  int_p <- sum(ints[mh_p]/total.int)
  if (length(mh_p) == 0) ion_p <- NA
  if (length(mh_p) == 0) int_p <- 0
  diagnostics <- list(ion_p, int_p)
}

#' analog of calcualte diagnostic pairs for non-cleavable crosslinkers
#'
#' @param datTab
#'
#' @return
#' @export
#'
#' @examples
calculateDiagnosticPairsNonCleavable <- function(datTab) {
  datTab <- datTab %>%
    mutate(diagnostics.1 = map2(MSMS.Ions.1, MSMS.Intensities.1, getDiagnosticMatchesNonCleavable),
           diagnostics.2 = map2(MSMS.Ions.2, MSMS.Intensities.2, getDiagnosticMatchesNonCleavable),
           ppPair = map2_lgl(diagnostics.1, diagnostics.2, function(x, y) {
             if_else(x[[2]] > 0.01 & y[[2]] > 0.01, T, F)}),
           diagnosticPair = ppPair
    )
}

#' parses the list matched product ions and translate to peptide bond position.
#'
#' @param msms.ions list of product ion matches found in Search Compare output
#' @param pep.len length of peptide
#'
#' @return
#' @export
#'
#' @examples
getProductIonMatches <- function(msms.ions, pep.len) {
  # will break if pep.len > 99
  ions <- unlist(str_split(msms.ions, ";"))
  n_indicies <- str_extract_all(ions, "(?<=^(b|c).?)[[0-9]]{1,2}(?=($|\\())") %>%
    unlist %>%
    unique %>%
    as.numeric
  c_indicies <- str_extract_all(ions, "(?<=^(y|z).?)[[0-9]]{1,2}(?=($|\\())") %>%
    unlist %>%
    unique %>%
    as.numeric
  c_indicies <- pep.len - c_indicies %>%
    sort
  ion_indicies <- union(n_indicies, c_indicies) %>% sort
  return(ion_indicies)
}

#' calculates number of distinct bonds cleavages observed as product ions.
#'
#' @param datTab
#'
#' @return
#' @export
#'
#' @examples
calculateProductIons <- function(datTab) {
  datTab <- datTab %>%
    mutate(numProdIons.1 = map2_int(MSMS.Ions.1, Len.Pep.1, function(x,y) length(getProductIonMatches(x,y))),
           numProdIons.2 = map2_int(MSMS.Ions.2, Len.Pep.2, function(x,y) length(getProductIonMatches(x,y))))
}

#' exclude CSMs below a given number of bond cleavages matched per peptide.
#'
#' @param datTab
#' @param minProducts.1 number of bonds cleaved in peptide 1
#' @param minProducts.2 number of bonds cleaved in peptide 2
#'
#' @return
#' @export
#'
#' @examples
productIonFilter <- function(datTab, minProducts.1 = 3, minProducts.2 = 3) {
  datTab <- datTab %>%
    filter(numProdIons.1 >= minProducts.1, numProdIons.2 >= minProducts.2)
}

#' prefilter CSMs by length of each peptide.
#'
#' @param datTab
#' @param minLen
#' @param maxLen
#'
#' @return
#' @export
#'
#' @examples
lengthFilter <- function(datTab, minLen, maxLen) {
  datTab <- datTab[datTab$Len.Pep.2 >= minLen & datTab$Len.Pep.1 >= minLen,]
  datTab <- datTab[datTab$Len.Pep.2 <= maxLen & datTab$Len.Pep.1 <= maxLen,]
  return(datTab)
}

#' prefilter CSMs by Prospector Score of each peptide
#'
#' @param datTab
#' @param minScore
#'
#' @return
#' @export
#'
#' @examples
scoreFilter <- function(datTab, minScore=0) {
  datTab <- datTab[datTab$Sc.1 >= minScore & datTab$Sc.2 >= minScore,]
  return(datTab)
}
