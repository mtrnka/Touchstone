#' @import rlang

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))

#' Read Search Compare Output from Protein Prospector
#'
#' Function for importing Crosslinked Search Compare output into Touchstone.
#' Search Compare output should be saved in tab delimited text format ideally
#' created with tstone specific param files, such as `tstoneMS2.json`.
#'
#' @param inputFile A tab delimited search compare file generated using tstoneMS2.json param file
#' @param minPepLen Minimum length of each peptide in amino acids.
#' @param minPepScore Minimum Prospector score of each peptide.
#' @param minScoreDiff Minimum Score Difference for the CSM.
#' @param minIons Mimimum number of product ions matched for each peptide.
#' @return A data frame
#' @export
#'
readProspectorXLOutput <- function(inputFile, minPepLen = 3, minPepScore = 0, minScoreDiff = 0, minIons = 0){
  datTab <- readr::read_tsv(inputFile, guess_max = 10000)
  header <- names(datTab) %>%
    stringr::str_replace("_[[0-9]]$", "") %>%
    stringr::str_replace_all("[[:space:]]",".") %>%
    stringr::str_replace_all("#", "Num") %>%
    stringr::str_replace_all("%", "Perc")
  acc_pos <- stringr::str_which(header, "Acc")
  header[acc_pos] <- c("Acc.1", "Acc.2")
  spec_pos <- stringr::str_which(header, "Species")
  header[spec_pos] <- c("Species.1", "Species.2")
  prot_pos <- stringr::str_which(header, "Protein.Name")
  header[prot_pos] <- c("Protein.1", "Protein.2")
  if (sum(stringr::str_detect(header, "Protein.MW")) == 2) {
    mw_pos <- stringr::str_which(header, "Protein.MW")
    header[mw_pos] <- c("MW.1", "MW.2")
  }
  if (sum(stringr::str_detect(header, "Protein.Length")) == 2) {
    len_pos <- stringr::str_which(header, "Protein.Length")
    header[len_pos] <- c("Protein.len.1", "Protein.len.2")
  }
  if (sum(stringr::str_detect(header, "MSMS.Ions")) == 2) {
    len_pos <- stringr::str_which(header, "MSMS.Ions")
    header[len_pos] <- c("MSMS.Ions.1", "MSMS.Ions.2")
  }
  if (sum(stringr::str_detect(header, "MSMS.M/Zs")) == 2) {
    len_pos <- stringr::str_which(header, "MSMS.M/Zs")
    header[len_pos] <- c("MSMS.MZs.1", "MSMS.MZs.2")
  }
  if (sum(stringr::str_detect(header, "MSMS.Intensities")) == 2) {
    len_pos <- stringr::str_which(header, "MSMS.Intensities")
    header[len_pos] <- c("MSMS.Intensities.1", "MSMS.Intensities.2")
  }
  if (sum(stringr::str_detect(header, "MSMS.Errors")) == 2) {
    len_pos <- stringr::str_which(header, "MSMS.Errors")
    header[len_pos] <- c("MSMS.Errors.1", "MSMS.Errors.2")
  }
  if (sum(stringr::str_detect(header, "Perc.Bond.Cleavage")) == 2) {
    len_pos <- stringr::str_which(header, "Perc.Bond.Cleavage")
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
    Acc.1 = as.character(.data$Acc.1),
    Acc.2 = as.character(.data$Acc.2)
  )
  datTab <- datTab %>%
    mutate(Protein.1 = if_else(is.na(.data$Protein.1), .data$Acc.1, .data$Protein.1),
           Protein.2 = if_else(is.na(.data$Protein.2), .data$Acc.2, .data$Protein.2))
  datTab <- calculateDecoys(datTab)
  datTab <- assignXLinkClass(datTab)
  datTab <- calculatePercentMatched(datTab)
  datTab <- calculatePeptideLengths(datTab)
  datTab <- lengthFilter(datTab, minLen = minPepLen, maxLen = 40)
  datTab <- scoreFilter(datTab, minScore = minPepScore)
  datTab <- datTab %>%
    filter(.data$Score.Diff >= minScoreDiff)
  datTab <- calculatePairs(datTab)
  if ("numProdIons.1" %in% names(datTab) & "numProdIons.1" %in% names(datTab)) {
    datTab <- productIonFilter(datTab, minProducts.1 = minIons, minProducts.2 = minIons) }
  datTab <- datTab %>%
    filter(.data$xlinkedResPair != "0.decoy::0.decoy")
  return(datTab)
}

#' Assign decoy/target status for each line of a crosslinked dataframe.
#'
#' @param datTab Parsed CLMS search results.
#' @return A data frame
#' @export
#'
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

#' Determine residue, peptide,and protein pairs
#'
#' @param datTab Parsed CLMS search results.
#' @return A data frame
#' @export
#'
calculatePairs <- function(datTab){
  datTab <- datTab %>%
    mutate(Acc.1 = as.character(.data$Acc.1),
           Acc.2 = as.character(.data$Acc.2))
  datTab$xlinkedProtPair <- ifelse(
    datTab$Decoy2 == "Decoy",
    ifelse(datTab$Protein.1 <= datTab$Protein.2,
           paste("decoy", datTab$Protein.1, datTab$Protein.2, sep="::"),
           paste("decoy", datTab$Protein.2, datTab$Protein.1, sep="::")),
    ifelse(datTab$Acc.1 <= datTab$Acc.2,
           paste(datTab$Acc.1, datTab$Acc.2, sep="::"),
           paste(datTab$Acc.2, datTab$Acc.1, sep="::"))
  )
  accs <- unique(c(datTab$Acc.1, datTab$Acc.2)) %>%
    stringr::str_sort()
  datTab <- datTab %>%
    mutate(Acc.1 = factor(.data$Acc.1, levels = accs),
           Acc.2 = factor(.data$Acc.2, levels = accs))
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
    add_count(.data$xlinkedResPair, name="numCSM") %>%
    add_count(.data$xlinkedResPair, name="wtCSM", wt = .data$Score.Diff)
  uniqueProtCount <- datTab %>%
    select(.data$xlinkedProtPair, .data$xlinkedResPair, .data$Score.Diff) %>%
    group_by(.data$xlinkedProtPair, .data$xlinkedResPair) %>%
    summarize(max.sd = max(.data$Score.Diff), .groups="drop") %>%
    group_by(.data$xlinkedProtPair) %>%
    add_count(name = "numURP") %>%
    add_count(name = "wtURP", wt=.data$max.sd) %>%
    ungroup()
  datTab <- left_join(select(datTab, -any_of("numURP")), uniqueProtCount, by=c("xlinkedProtPair","xlinkedResPair"))
  if ("Module.1" %in% names(datTab) & "Module.2" %in% names(datTab)) {
    datTab <- datTab %>%
      mutate(Module.1 = as.character(.data$Module.1),
             Module.2 = as.character(.data$Module.2))
    mods <- unique(c(datTab$Module.1, datTab$Module.2)) %>% stringr::str_sort()
    datTab <- datTab %>%
      mutate(xlinkedModulPair = ifelse(.data$Module.1 <= .data$Module.2,
                                       stringr::str_c(.data$Module.1, .data$Module.2, sep="::"),
                                       stringr::str_c(.data$Module.2, .data$Module.1, sep="::"))
      ) %>%
      mutate(xlinkedModulPair = forcats::as_factor(.data$xlinkedModulPair),
             Module.1 = factor(.data$Module.1, levels = mods),
             Module.2 = factor(.data$Module.2, levels = mods)) #%>%
    #         add_count(xlinkedModulPair, name="numMPSM")
  }
  if ("MSMS.Ions.1" %in% names(datTab) & "MSMS.Ions.2" %in% names(datTab)) {
    datTab <- calculateProductIons(datTab)
  }
  return(datTab)
}

#' Count number of CSMs per URP and the number of URPs per PP
#'
#' @param datTab Parsed CLMS search results
#' @return A data frame
#' @export
#'
addXlinkCount <- function(datTab) {
  datTab <- datTab %>%
    add_count(.data$xlinkedResPair, name="numCSM")
  uniqueProtCount <- datTab %>%
    select(.data$xlinkedProtPair, .data$xlinkedResPair) %>%
    group_by(.data$xlinkedProtPair, .data$xlinkedResPair) %>%
    summarize(.groups="drop") %>%
    group_by(.data$xlinkedProtPair) %>%
    count(name = "numURP") %>%
    ungroup()
  datTab <- left_join(select(datTab, -any_of("numURP")), uniqueProtCount, by="xlinkedProtPair")
}

#' Determines whether crosslink is inter-protein (hetromeric) or intra-protein (homomeric)
#'
#' @param datTab Parsed CLMS search results.
#' @return A data frame
#' @export
#'
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

#' Determines the percent of the total ion current in a product ion spectra that
#' is matched by the crosslinked peptide annotation.
#'
#' If MS-Product is reporting `%TIC matched` that value is used. Otherwise
#' the value is calculated based on the number of ions matches.
#'
#' @param datTab Parsed CLMS search results.
#' @return A data frame
#' @export
#'
calculatePercentMatched <- function(datTab) {
  if ("Match.Int" %in% names(datTab)) {
    datTab$percMatched <- datTab$Match.Int
  } else if ("Num.Pks" %in% names(datTab) & "Num.Unmat" %in% names(datTab)) {
    num.Matched <- datTab$Num.Pks - datTab$Num.Unmat
    datTab$percMatch <- num.Matched / datTab$Num.Pks
  }
  return(datTab)
}

#' Calculate length of each peptide in amino acid residues.
#'
#' @param datTab Parsed CLMS search results.
#' @return A data frame
#' @export
#'
calculatePeptideLengths <- function(datTab) {
  datTab$Len.Pep.1 <- nchar(datTab$DB.Peptide.1)
  datTab$Len.Pep.2 <- nchar(datTab$DB.Peptide.2)
  return(datTab)
}

#' Are diagnostic/reporter ions (MH*/#) from cleavable crosslinking reagents detected?
#'
#' @param msms.ions Columns containing list of detected ions, must be present in Search Compare output.
#' @param msms.int  Columns with corresponding intensity values for detected ions.
#' @return A list of diagnostic ions and their intensities.
#' @seealso [calculateDiagnosticPairs()]
getDiagnosticMatches <- function(msms.ions, msms.int) {
  ions <- unlist(stringr::str_split(msms.ions, ";"))
  ints <- as.numeric(unlist(stringr::str_split(msms.int, ";")))
  total.int <- sum(ints, na.rm = T)
  #   max.int <- max(ints, na.rm=T)
  mh_a <- stringr::str_which(ions, "MH\\#")
  ion_a <- ions[mh_a]
  int_a <- sum(ints[mh_a]/total.int)
  mh_s <- stringr::str_which(ions, "MH\\*")
  ion_s <- ions[mh_s]
  int_s <- sum(ints[mh_s]/total.int)
  if (length(mh_a) == 0) ion_a <- NA
  if (length(mh_a) == 0) int_a <- 0
  if (length(mh_s) == 0) ion_s <- NA
  if (length(mh_a) == 0) int_a <- 0
  diagnostics <- list(ion_a, int_a, ion_s, int_s)
}

#' Determines which crosslinks contain matches to diagnostic/reporter ions (MH*/#)
#' from cleavable crosslinking reagents.
#'
#' @param datTab Parsed CLMS search results.
#' @param int.thresh Minimum intensity for diagnostic ion.
#' @return A data frame
#' @seealso [getDiagnosticMatches()]
#' @export
#'
calculateDiagnosticPairs <- function(datTab, int.thresh = 0.05) {
  datTab <- datTab %>%
    mutate(diagnostics.1 = purrr::map2(.data$MSMS.Ions.1, .data$MSMS.Intensities.1, getDiagnosticMatches),
           diagnostics.2 = purrr::map2(.data$MSMS.Ions.2, .data$MSMS.Intensities.2, getDiagnosticMatches),
           aaPair = purrr::map2_lgl(.data$diagnostics.1, .data$diagnostics.2, function(x, y) {
             if_else(x[[2]] > int.thresh & y[[2]] > int.thresh, T, F)}),
           ttPair = purrr::map2_lgl(.data$diagnostics.1, .data$diagnostics.2, function(x, y) {
             if_else(x[[4]] > int.thresh & y[[4]] > int.thresh, T, F)}),
           atPair = purrr::map2_lgl(.data$diagnostics.1, .data$diagnostics.2, function(x, y) {
             if_else(x[[2]] > int.thresh & y[[4]] > int.thresh, T, F)}),
           taPair = purrr::map2_lgl(.data$diagnostics.1, .data$diagnostics.2, function(x, y) {
             if_else(x[[4]] > int.thresh & y[[2]] > int.thresh, T, F)})
    )
  datTab <- datTab %>%
    mutate(diagnosticPair = purrr::pmap_lgl(list(datTab$aaPair,datTab$ttPair,datTab$taPair,datTab$atPair),
                                            function(a,b,c,d) sum(a,b,c,d, na.rm=T) >= 1)
    )
}

#' Determines if diagnostic/reporter ions from cleavable crosslinkers are present
#' across linked MS2 scans.
#'
#' For instance, when sequential MS2-CID and MS2-ETD are triggered on the same precursor ion.
#' `calculateCrossDiagnosticPairs()` looks for the `MH*/#-ion` in the CID scan, but
#' annotated it as having been found in both scans. Requires a `masterScanFile` typically
#' generated from PAVA processed mgf files.
#'
#' This function is mostly intended for internal testing and development.
#'
#' @param datTab.etd Parsed CLMS table that is not expected to contain diagnostic ions (does not need to be ETD)
#' @param datTab.cid Parsed CLMS table that is expected to contain diagnostic ions (typical MS2-CID)
#' @param masterScanFile Master scan file that charts the dependencies between different scans.
#' @return A data frame
#' @seealso [calculateDiagnosticPairs()]
#' @export
#'
calculateCrossDiagnosticPairs <- function(datTab.etd, datTab.cid, masterScanFile) {
  datTab.etd <- datTab.etd %>%
    mutate(Fraction = stringr::str_replace_all(.data$Fraction, "_FTMSms2[[a-z]]+$", ""))
  datTab.cid <- datTab.cid %>%
    mutate(Fraction = stringr::str_replace_all(.data$Fraction, "_FTMSms2[[a-z]]+$", ""))
  masterScanFile <- masterScanFile %>%
    mutate(fraction = stringr::str_replace_all(.data$fraction, "_FTMSms2[[a-z]]+$", ""))
  datTab.etd <- left_join(datTab.etd, select(masterScanFile, .data$fraction, .data$ms2cidScanNo, .data$ms2etdScanNo),
                          by=c("Fraction"="fraction", "MSMS.Info"="ms2etdScanNo"))
  datTab.etd <- left_join(datTab.etd, select(datTab.cid, .data$Fraction, .data$MSMS.Info, .data$diagnosticPair),
                          by=c("Fraction"="Fraction", "ms2cidScanNo"="MSMS.Info"),
                          suffix=c("",".cross"))
  datTab.etd <- datTab.etd %>%
    mutate(diagnosticPair.cross = if_else(is.na(.data$diagnosticPair.cross), F, .data$diagnosticPair.cross))
  return(datTab.etd)
}

#' Determines which crosslinks contain matches to diagnostic/reporter ions (P,PL,PLK)
#' from non-cleavable crosslinking reagents such as DSS/BS3.
#'
#' Analog of `getDiagnosticMatches()` for non-cleavable crosslinkers. Looks for P, PL, PLK type ions.
#'
#' @param msms.ions Columns containing list of detected ions, must be present in Search Compare output.
#' @param msms.int  Columns with corresponding intensity vaues fro detected ions.
#' @return A list of p-type ions and their intensities
#' @seealso [getDiagnosticMatches()], [calculateDiagnosticPairsNonCleavable]
getDiagnosticMatchesNonCleavable <- function(msms.ions, msms.int) {
  ions <- unlist(stringr::str_split(msms.ions, ";"))
  ints <- as.numeric(unlist(stringr::str_split(msms.int, ";")))
  total.int <- sum(ints, na.rm = T)
  mh_p <- stringr::str_which(ions, "P[^\\+\\-]")
  ion_p <- ions[mh_p]
  int_p <- sum(ints[mh_p]/total.int)
  if (length(mh_p) == 0) ion_p <- NA
  if (length(mh_p) == 0) int_p <- 0
  diagnostics <- list(ion_p, int_p)
}

#' Determines which crosslinks contain matches immonium related ions (P, PL, PLK)
#' resulting from dissociation of non-cleavable crosslinking reagents such as BS3.
#'
#' Analog of `calculateDiagnosticPairs()` for non-cleavable crosslinkers
#'
#' @param datTab parsed CLMS search results
#' @return A data frame
#' @seealso [getDiagnosticMatchesNonCleavable()], [calculateDiagnosticPairs()]
#' @export
#'
calculateDiagnosticPairsNonCleavable <- function(datTab) {
  datTab <- datTab %>%
    mutate(diagnostics.1 = purrr::map2(.data$MSMS.Ions.1, .data$MSMS.Intensities.1, getDiagnosticMatchesNonCleavable),
           diagnostics.2 = purrr::map2(.data$MSMS.Ions.2, .data$MSMS.Intensities.2, getDiagnosticMatchesNonCleavable),
           ppPair = purrr::map2_lgl(.data$diagnostics.1, .data$diagnostics.2, function(x, y) {
             if_else(x[[2]] > 0.01 & y[[2]] > 0.01, T, F)}),
           diagnosticPair = .data$ppPair
    )
}

#' Parses the list matched product ions and translate to peptide bond position.
#'
#' @param msms.ions List of product ion matches found in Search Compare output
#' @param pep.len Length of peptide
#' @return A numeric vector of bond cleavage indicies
#' @seealso [calculateProductIons()]
getProductIonMatches <- function(msms.ions, pep.len) {
  # will break if pep.len > 99
  ions <- unlist(stringr::str_split(msms.ions, ";"))
  n_indicies <- stringr::str_extract_all(ions, "(?<=^[bc][\\*\\#]?)([[0-9]]{1,2})") %>%
    unlist %>%
    unique %>%
    as.numeric
  c_indicies <- stringr::str_extract_all(ions, "(?<=^[yz][\\*\\#]?)([[0-9]]{1,2})") %>%
    unlist %>%
    unique %>%
    as.numeric
  c_indicies <- pep.len - c_indicies %>%
    sort
  ion_indicies <- union(n_indicies, c_indicies) %>% sort
  return(ion_indicies)
}

#' Calculates the number of distinct backbone bond-cleavages observed as in the product ion spectrum.
#'
#' N-terminal ions (b-, c-) and C-terminal ions (y-, z-) are translated into a bond index.
#' It is the number of unique bond indicies that is calculated here.
#'
#' @param datTab Parsed CLMS search results.
#' @return A data frame
#' @export
#'
calculateProductIons <- function(datTab) {
  datTab <- datTab %>%
    mutate(numProdIons.1 = purrr::map2_int(.data$MSMS.Ions.1, .data$Len.Pep.1, function(x,y) length(getProductIonMatches(x,y))),
           numProdIons.2 = purrr::map2_int(.data$MSMS.Ions.2, .data$Len.Pep.2, function(x,y) length(getProductIonMatches(x,y))))
}

#' Exclude crosslinks or CSMs below a threshold number of bond cleavages matched per peptide.
#'
#' @param datTab Parsed CLMS search results.
#' @param minProducts.1 Number of bonds cleaved in peptide 1.
#' @param minProducts.2 Number of bonds cleaved in peptide 2.
#' @return A data frame
#' @seealso [calculateProductIons()], [getProductIonMatches()]
#' @export
#'
productIonFilter <- function(datTab, minProducts.1 = 3, minProducts.2 = 3) {
  datTab <- datTab %>%
    filter(.data$numProdIons.1 >= minProducts.1, .data$numProdIons.2 >= minProducts.2)
}

#' Prefilter crosslinks or CSMs by length of each peptide.
#'
#' @param datTab Parsed CLMS search results.
#' @param minLen Mimum peptide length, as number of amino acid residues.
#' @param maxLen Maximum peptide length, as number of amino acid residues.
#' @return A data frame
#' @export
#'
lengthFilter <- function(datTab, minLen, maxLen) {
  datTab <- datTab[datTab$Len.Pep.2 >= minLen & datTab$Len.Pep.1 >= minLen,]
  datTab <- datTab[datTab$Len.Pep.2 <= maxLen & datTab$Len.Pep.1 <= maxLen,]
  return(datTab)
}

#' Prefilter crosslinks or CSMs by Prospector Score of each peptide.
#'
#' @param datTab Parsed CLMS search results.
#' @param minScore Minimum Protein Prospector Score for each individual peptide
#' @return A data frame
#' @export
#'
scoreFilter <- function(datTab, minScore=0) {
  datTab <- datTab[datTab$Sc.1 >= minScore & datTab$Sc.2 >= minScore,]
  return(datTab)
}

#' Parent function for summarization. CLMS data is grouped by summarization variable
#' and the best scoring CSM from the group is chosen as the representative
#'
#' @param datTab Parsed CLMS search results.
#' @param summarizationVar variable to summarize on
#' @param classifier Column name in `datTab` used as the classifier to use to rank hits.
#' @param retainGroups Should pre-existing grouping of datTab be retained?
#' @keywords internal
#' @returns A data frame
#' @export
#'
bestPair <- function(datTab,
                     summarizationVar=.data$xlinkedResPair,
                     classifier=.data$SVM.score,
                     retainGroups=T){
  datTab <- datTab %>%
    group_by({{summarizationVar}}, .add=retainGroups) %>%
    filter(n() > 0) %>%
    filter({{classifier}}==max({{classifier}})) %>%
    slice(1) %>%
    ungroup()
  return(datTab)
}

#' Summarize on Unique Residue Pairs (URPs)
#'
#' The default behaviour is to summarize residue-pairs within existing groups.
#' This would typically be different experimental conditions. To overwrite and
#' calculate the residue-pairs on the entire data table seq the `retainGroup` option
#' to `False`.
#' @param datTab Parsed CLMS search results.
#' @param classifier Column name in `datTab` used as the classifier to use to rank hits.
#' @param retainGroups Will retain pre-existing group in the data frame when apply the summarization.
#' @return A data frame
#' @seealso [bestPepPair()], [bestProtPair()], [bestModPair()]
#' @export
#'
bestResPair <- function(datTab, classifier=.data$SVM.score, retainGroups=T){
  classifier.quo <- enquo(classifier)
  bestPair(datTab, summarizationVar=.data$xlinkedResPair, !!classifier.quo, retainGroups)
}

# bestResPair <- function(datTab, classifier=.data$SVM.score, retainGroups=T){
#   quoClass <- enquo(classifier)
#   datTab <- datTab %>%
#     group_by(.data$xlinkedResPair, .add=retainGroups) %>%
#     filter(n() > 0) %>%
#     filter(!! quoClass==max(!! quoClass)) %>%
#     slice(1) %>%
#     ungroup()
#   return(datTab)
# }

#' Summarize on Peptide Pairs
#'
#' The default behaviour is to summarize peptide-pairs within existing groups.
#' This would typically be different experimental conditions. To overwrite and
#' calculate the peptide-pairs on the entire data table seq the `retainGroup` option
#' to `False`.
#' @param datTab Parsed CLMS search results.
#' @param classifier Column name in `datTab` used as the classifier to use to rank hits.
#' @param retainGroups Will retain pre-existing group in the data frame when apply the summarization.
#' @return A data frame
#' @seealso [bestResPair()], [bestProtPair()], [bestProtPair()]
#' @export
#'
bestPepPair <- function(datTab, classifier=.data$SVM.score, retainGroups=T){
  classifier.quo <- enquo(classifier)
  bestPair(datTab, summarizationVar=.data$xlinkedPepPair, !!classifier.quo, retainGroups)
}

#' Summarize on Protein Pairs (PPs)
#'
#' The default behaviour is to summarize protein-pairs within existing groups.
#' This would typically be different experimental conditions. To overwrite and
#' calculate the protein-pairs on the entire data table seq the `retainGroup` option
#' to `False`.
#' @param datTab Parsed CLMS search results.
#' @param classifier Column name in `datTab` used as the classifier to use to rank hits.
#' @param retainGroups Will retain pre-existing group in the data frame when apply the summarization.
#' @return A data frame
#' @seealso [bestResPair()], [bestPepPair()], [bestModPair()]
#' @export
#'
bestProtPair <- function(datTab, classifier=.data$SVM.score, retainGroups=T){
  classifier.quo <- enquo(classifier)
  bestPair(datTab, summarizationVar=.data$xlinkedProtPair, !!classifier.quo, retainGroups)
}

#' Summarize on Module Pairs (MPs)
#'
#' The default behaviour is to summarize moudle-pairs within existing groups.
#' This would typically be different experimental conditions. To overwrite and
#' calculate the module-pairs on the entire data table seq the `retainGroup` option
#' to `False`.
#' @param datTab Parsed CLMS search results.
#' @param classifier Column name in `datTab` used as the classifier to use to rank hits.
#' @param retainGroups Will retain pre-existing group in the data frame when apply the summarization.
#' @return A data frame
#' @seealso [bestResPair()], [bestPepPair()], [bestProtPair()]
#' @export
#'
bestModPair <- function(datTab, classifier=.data$SVM.score, retainGroups=T){
  classifier.quo <- enquo(classifier)
  bestPair(datTab, summarizationVar=.data$xlinkedModulPair, !!classifier.quo, retainGroups)
}

