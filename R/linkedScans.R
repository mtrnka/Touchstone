#' Calculate precursor mass from m/z and charge
#'
#' @param pepMass m/z from mgf header.
#' @param charge  z from mgf header.
#' @return A numeric value
#' @export
#'
calculatePrecMass <- function(pepMass, charge) {
  precMass <- pepMass * charge - atomic.weight.da[["proton"]] * charge
  return(precMass)
}

#' Calculate the mass error between theoretical and measured values in ppm.
#'
#' @param targetMass theoretical mass of peptide.
#' @param testMass experimental mass of peptide.
#' @return A numeric value
#' @export
#'
getMassError <- function(targetMass, testMass) {
  # reports mass error in ppm
  massDiff <- 1e6 * (targetMass - testMass) / targetMass
  return(massDiff)
}

#' Gathers information from the headers of an MGF formatted peaklist file.
#'
#' This function assumes PAVA formatted mgf file headers. There is some variability
#' in the header format depending on the program used to generate the mgf and I
#' haven't tested other programs output. `getIndicesforFile()` is looking for
#' the Scan No in the `TITLE` linke as well as a `Master_Scan_Number()`.
#'
#' @param file Path to the mgf peaklist.
#' @return A data frame
#' @seealso [createMasterScanFile()], [createMasterScanFile2()]
#' @export
#'
getIndicesForFile <- function(file) {
  scanNoLine <- system2("grep", c("'TITLE=Scan'", file), stdout = T)
  scanNo <- stringr::str_extract(scanNoLine, "(?<=TITLE\\=Scan\\s)[[0-9]]+")
  scanNo <- as.integer(scanNo)
  retTime <- stringr::str_extract(scanNoLine, "(?<=\\(rt\\=)[[0-9]]+\\.[[0-9]]+")
  retTime <- as.numeric(retTime)
  masterScanNo <- system2("grep", c("'Master_Scan_Number'", file), stdout = T)
  masterScanNo <- stringr::str_extract(masterScanNo, "(?<=Master_Scan_Number\\=)([[0-9]]+)")
  masterScanNo <- as.integer(masterScanNo)
  pepMass <- system2("grep", c("'PEPMASS'", file), stdout = T)
  pepMass <- stringr::str_extract(pepMass, "(?<=PEPMASS\\=)[[0-9]]+\\.[[0-9]]+")
  pepMass <- as.numeric(pepMass)
  charge <- system2("grep", c("'CHARGE'", file), stdout = T)
  charge <- stringr::str_extract(charge, "(?<=CHARGE\\=)[[0-9]](?=\\+)")
  charge <- as.integer(charge)
  extract <- data.frame(scanNo=scanNo, masterScanNo=masterScanNo,
                        pepmass=pepMass, charge=charge, retTime = retTime)
  return(extract)
}

#' Create a master scan file which maps scan dependencies in mgf peakists.
#'
#' `createMasterScanFile()` is for scan sequences that contain an MS3 scan. May
#' contain an additional MS2 scan in the acquisition sequence (assumed to be ETD
#' in this function, but can be any). Requires using the PAVA formating conventions
#' for the peaklists filenames. For instance, PAVA creates separate peaklist
#' files for each scan type and appends a suffic such as: `_FTMSms2ethcd.txt`
#' for high-res MS2-EThcD scans, or `_ITMSms3cid.txt` for low-res MS3-CID scans.
#' The base filename is the same for all the linked peaklist files, corresponding
#' to the name of the .raw file.
#'
#' @param ms3PAVApeaklist path to MS3 peaklist.
#' @param ms2PAVApeaklistCID path to MS2 CID peaklist (parent of MS3 scans).
#' @param ms2PAVApeaklistETD path to MS2 ETD peaklist (optional).
#' @return A data frame
#' @seealso [getIndicesForFile()],[createMasterScanFile2()]
#' @export
#'
createMasterScanFile <- function(ms3PAVApeaklist,
                                 ms2PAVApeaklistCID,
                                 ms2PAVApeaklistETD=NA) {
  #Peaklists need to be pava generated. Can take vector input.
  numFiles <- length(ms3PAVApeaklist)
  if (numFiles != length(ms2PAVApeaklistCID)) {
    stop("Different numbers of MS2 and MS3 files provided")
  }
  masterScans <- list()
  for (i in 1:numFiles) {
    ms3masterscansCID <- getIndicesForFile(ms3PAVApeaklist[i])
    names(ms3masterscansCID) <- c("ms3ScanNo", "ms2cidScanNo",
                                  "pepmass.ms3", "charge.ms3", "retTime.ms3")
    ms2masterscansCID <- getIndicesForFile(ms2PAVApeaklistCID[i])
    if (!is.na(ms2PAVApeaklistETD[i])) {
      ms2masterscansETD <- getIndicesForFile(ms2PAVApeaklistETD[i])
      ms2masterscans <- full_join(ms2masterscansCID, ms2masterscansETD,
                                  by=c("masterScanNo", "pepmass","charge"))
      names(ms2masterscans) <- c("ms2cidScanNo", "ms1MasterScanNo",
                                 "pepmass.ms2", "charge.ms2", "ms2etdScanNo", "retTime.ms2")
    } else {
      ms2masterscans <- ms2masterscansCID
      names(ms2masterscans) <- c("ms2cidScanNo", "ms1MasterScanNo",
                                 "pepmass.ms2", "charge.ms2", "retTime.ms2")
    }
    ms2masterscans$precMass.ms2 <- purrr::map2_dbl(ms2masterscans$pepmass.ms2,
                                                   ms2masterscans$charge.ms2,
                                                   calculatePrecMass)
    ms3masterscansCID$precMass.ms3 <- purrr::map2_dbl(ms3masterscansCID$pepmass.ms3,
                                                      ms3masterscansCID$charge.ms3,
                                                      calculatePrecMass)
    ms2masterscans <- full_join(ms2masterscans, ms3masterscansCID,
                                by="ms2cidScanNo")
    if (!is.na(ms2PAVApeaklistETD[i])) {
      ms2masterscans <- ms2masterscans[,c(10000,2,1,5:6,3:4,9,7:8)]
    } else {
      ms2masterscans <- ms2masterscans[,c(2,1,7,3:6,8:11)]
    }
    fract <- stringr::str_replace(basename(ms3PAVApeaklist[i]), "\\.(txt|mgf)$", "")
    ms2masterscans$fraction <- fract
    masterScans[[fract]] <- ms2masterscans
  }
  masterScans <- do.call(rbind, masterScans)
  return(masterScans)
}

#' Create a master scan file which maps scan dependencies in mgf peakists.  For scan sequences with two MS2 files.
#'
#' `createMasterScanFile2()` is for scan sequences that contain two dependent
#' MS2 scans and no MS3 scans. They are annotted in the help file as MS2-CID and
#' MS2-ETD peaklists but can be any type. Requires using the PAVA formating conventions
#' for the peaklists filenames. For instance, PAVA creates separate peaklist
#' files for each scan type and appends a suffic such as: `_FTMSms2ethcd.txt`
#' for high-res MS2-EThcD scans, or `_FTMSms2cid.txt` for high-res MS2-CID scans.
#' The base filename is the same for all the linked peaklist files, corresponding
#' to the name of the .raw file.
#' @param ms2PAVApeaklistCID path to MS2 CID peaklist (parent of MS3 scans).
#' @param ms2PAVApeaklistETD path to MS2 ETD peaklist (optional).
#' @return A data frame
#' @seealso [getIndicesForFile()],[createMasterScanFile2()]
#' @export
#'
createMasterScanFile2 <- function(ms2PAVApeaklistCID=NA,
                                  ms2PAVApeaklistETD=NA) {
  #Peaklists need to be pava generated. Can take vector input.
  numFiles <- length(ms2PAVApeaklistCID)
  if (numFiles != length(ms2PAVApeaklistETD)) {
    stop("Different numbers of MS2 files provided")
  }
  masterScans <- list()
  for (i in 1:numFiles) {
    ms2masterscansCID <- getIndicesForFile(ms2PAVApeaklistCID[i]) %>%
      rename(ms2cidScanNo = .data$scanNo)
    if (!is.na(ms2PAVApeaklistETD[i])) {
      ms2masterscansETD <- getIndicesForFile(ms2PAVApeaklistETD[i]) %>%
        rename(ms2etdScanNo = .data$scanNo)
      ms2masterscans <- full_join(ms2masterscansCID, ms2masterscansETD,
                                  by=c("masterScanNo", "pepmass","charge"),
                                  suffix=c(".cid", ".etd")) %>%
        select(.data$masterScanNo, .data$ms2cidScanNo, .data$ms2etdScanNo, .data$pepmass, .data$charge, everything()) %>%
        rename("pepmass.ms2" = .data$pepmass, "charge.ms2" = .data$charge)
    }
    fract <- stringr::str_replace(basename(ms2PAVApeaklistCID[i]), "\\.(txt|mgf)$", "")
    ms2masterscans$fraction <- fract
    masterScans[[fract]] <- ms2masterscans
  }
  masterScans <- do.call(rbind, masterScans)
  return(masterScans)
}

#' Reads Prospector formatted search compare output from MS3 searches of linear peptides.
#'
#' This is a developmental function for generating results from MS3-based search
#' strategies using cleavable linkers.
#'
#' @param ms3SearchTable Search Compare Output from a linear Protein Prospector linear peptide search.
#' @return A data frame.
#' @export
#'
readMS3results <- function(ms3SearchTable){
  ms3results <- read_tsv(ms3SearchTable, skip=2)
  if (!"Spectrum" %in% names(ms3results)) ms3results <- ms3results %>% add_column("Spectrum" = 1, .after="RT")
  ms3results <- ms3results %>%
    select(c("DB Peptide", "Peptide", "Protein Mods",
             "Fraction", "RT", "Spectrum", "MSMS Info", "Score",
             "Expect", "Protein Name", "Acc #", "Species"))
  names(ms3results) <- c("DB.Peptide.ms3", "Peptide.ms3", "Protein.Mods.ms3",
                         "Fraction.ms3", "RT.ms3", "Spectrum.ms3", "ms3ScanNo",
                         "Score.ms3", "Expect.ms3", "Protein.ms3", "Acc.ms3",
                         "Species.ms3")
  ms3results$Decoy.ms3 <- ifelse(stringr::str_detect(ms3results$Acc.ms3, "decoy"),
                                 "DECOY", "Target")
  ms3results <- ms3results %>%
    #        filter(stringr::str_detect(Protein.Mods.ms3, "Xlink:DSSO_[[as]]_fragment"))
    filter(stringr::str_detect(.data$Protein.Mods.ms3, "XL:A-[[a-zA-Z\\(\\)]]+"))
  return(ms3results)
}

#' Add information from the master scan file to an MS3 search result
#'
#' Adds information from scan-dependencies generated with `createMasterScanFile()`.
#' @param ms3results ms3 level results table as parsed by `readMS3results()`
#' @param masterScanFile master scan file as created by `createMasterScanFile()`
#' @return A data frame
#' @seealso [createMasterScanFile()],[readMS3results()]
#'
addMasterScanInfo <- function(ms3results, masterScanFile) {
  ms3results <- left_join(ms3results, masterScanFile,
                          by = c("ms3ScanNo"="ms3ScanNo", "Fraction.ms3"="fraction"))
  ms3results$precMass.ms3 <- purrr::map2_dbl(ms3results$pepmass.ms3,
                                             ms3results$charge.ms3, calculatePrecMass)
  return(ms3results)
}

#' Recombines linked MS3-search results into crosslinked Peptides.
#'
#' `processMS3xlinkResults()` takes Protein Prospector search results for linear
#' peptide from MS3-scans of cleavable crosslinkers and searches links scans for
#' possible crosslinked spectral matches (CSMs). This is developmental and currently
#' only works with DSSO. Requires more testing. Requires the MS3-search results
#' to have `Master Scan File` information appended.
#'
#' @param ms3results Parsed MS3-search results.
#' @param tol precursor Mass tolerance in ppm.
#' @return A data frame.
#' @seealso [ms3ionFragFind()], [addMasterScanInfo()]
#' @export
#'
processMS3xlinkResults <- function(ms3results, tol=25) {
  n_ms3 <- ms3results %>%
    group_by(.data$ms2cidScanNo, .data$Fraction.ms3) %>%
    nest()
  ms3crosslinksFlat <- purrr::pmap_dfr(list(n_ms3$Fraction.ms3,
                                            n_ms3$ms2cidScanNo,
                                            n_ms3$data), ms3ionFragFind, tol)
  ms3crosslinksFlatMis <- purrr::pmap_dfr(list(n_ms3$Fraction.ms3,
                                               n_ms3$ms2cidScanNo,
                                               n_ms3$data), ms3ionFragFind, tol, mis=atomic.weight.da[["proton"]])
  ms3crosslinksFlatMis2 <- purrr::pmap_dfr(list(n_ms3$Fraction.ms3,
                                                n_ms3$ms2cidScanNo,
                                                n_ms3$data), ms3ionFragFind, tol, mis=2*atomic.weight.da[["proton"]])
  ms3crosslinksFlatMisM1 <- purrr::pmap_dfr(list(n_ms3$Fraction.ms3,
                                                 n_ms3$ms2cidScanNo,
                                                 n_ms3$data), ms3ionFragFind, tol, mis=-1*atomic.weight.da[["proton"]])
  ms3crosslinksFlat <- rbind(ms3crosslinksFlatMisM1, ms3crosslinksFlatMis2, ms3crosslinksFlatMis, ms3crosslinksFlat)
  ms3crosslinksFlat$Score.Diff <- ifelse(ms3crosslinksFlat$Sc.2 < ms3crosslinksFlat$Sc.1,
                                         ms3crosslinksFlat$Sc.2, ms3crosslinksFlat$Sc.1)
  ms3crosslinksFlat <- assignXLinkClass(ms3crosslinksFlat)
  ms3crosslinksFlat <- lengthFilter(ms3crosslinksFlat, minLen = 3, maxLen = 40)
  ms3crosslinksFlat <- scoreFilter(ms3crosslinksFlat, minScore = 0)
  ms3crosslinksFlat <- calculateDecoys(ms3crosslinksFlat)
  ms3crosslinksFlat <- calculatePairs(ms3crosslinksFlat)
  if (!"distance" %in% names(ms3crosslinksFlat)) {
    ms3crosslinksFlat$distance <- NA_real_
  }
  return(ms3crosslinksFlat)
}

#' Takes linked MS3-level PSMs and tests whether there are complementary scans
#' that sum up to the precursor mass to detect putative crosslinks. This is
#' developmental and currently only works with DSSO. Requires more testing.
#' Requires the MS3-search results to have `Master Scan File` information appended.
#'
#' @param fraction Base name of RAW file.
#' @param masterScanNo Scan number of the triggering scan (typically MS2-CID).
#' @param scanGroupTibble Nested ms3 level linear search results.
#' @param tol precursor Mass tolerance in ppm.
#' @param mis misannotated precursor ion index.
#' @seealso [processMS3xlinkResults()], [addMasterScanInfo()]
#' @return A data frame
ms3ionFragFind <- function(fraction, masterScanNo, scanGroupTibble, tol = 25, mis = 0) {
  ms2Mass <- scanGroupTibble$precMass.ms2[1] - mis
  scanGroupSize <- nrow(scanGroupTibble)
  if (scanGroupSize < 2) {
    return(data.frame())
  }
  scanGroupTibble$alkene <- stringr::str_count(scanGroupTibble$Peptide.ms3, stringr::fixed("XL:A-Alkene"))
  scanGroupTibble$thiol <- stringr::str_count(scanGroupTibble$Peptide.ms3, stringr::fixed("XL:A-Thiol(Unsaturated)"))
  scanGroupTibble$sulfenic <- stringr::str_count(scanGroupTibble$Peptide.ms3, stringr::fixed("XL:A-Sulfenic"))

  #    scanGroupTibble$alkene <- stringr::stringr::str_detect(scanGroupTibble$Protein.Mods.ms3, "Alkene@[[0-9\\|]]+")
  #    scanGroupTibble$thiol <- stringr::stringr::str_detect(scanGroupTibble$Protein.Mods.ms3, "Thiol\\(Unsaturated\\)@[[0-9\\|]]+")
  #    scanGroupTibble$sulfenic <- stringr::stringr::str_detect(scanGroupTibble$Protein.Mods.ms3, "Sulfenic@[[0-9\\|]]+")
  summedPairs <- scanGroupTibble$precMass.ms3 %o% rep(1, scanGroupSize) +
    rep(1, scanGroupSize) %o% scanGroupTibble$precMass.ms3

  massErrors.AT <- getMassError(ms2Mass, summedPairs + H2O)
  matchesMS2 <- abs(massErrors.AT) <= tol
  allowed <- scanGroupTibble$alkene %o% scanGroupTibble$thiol
  ATpairs <- which(matchesMS2 * allowed == 1)

  massErrors.AS <- getMassError(ms2Mass, summedPairs)
  matchesMS2 <- abs(massErrors.AS) <= tol
  allowed <- scanGroupTibble$alkene %o% scanGroupTibble$sulfenic
  ASpairs <- which(matchesMS2 * allowed == 1)

  massErrors.AA <- getMassError(ms2Mass, summedPairs + atomic.weight.da[["sulfur"]] + H2O)
  matchesMS2 <- abs(massErrors.AA) <= tol
  allowed <- scanGroupTibble$alkene %o% scanGroupTibble$alkene
  AApairs <- which(matchesMS2 * allowed == 1)

  massErrors.TT <- getMassError(ms2Mass, summedPairs - atomic.weight.da[["sulfur"]] + H2O)
  matchesMS2 <- abs(massErrors.TT) <= tol
  allowed <- scanGroupTibble$thiol %o% scanGroupTibble$thiol
  TTpairs <- which(matchesMS2 * allowed == 1)

  numMatches <- length(c(ATpairs, ASpairs, AApairs, TTpairs))
  dimensions <- dim(allowed)
  if (numMatches==0) return(data.frame())
  placeHolder <- vector("list", numMatches)
  index = 1
  for (CSM in ATpairs) {
    massError <- massErrors.AT[CSM]
    peptidePair <- arrayInd(CSM, dimensions)
    flatMS3pair <- flattenMS3pair(scanGroupTibble[peptidePair[,1],],
                                  scanGroupTibble[peptidePair[,2],],
                                  fraction, masterScanNo,
                                  frag1="Alkene", frag2="Thiol\\(Unsaturated\\)")
    flatMS3pair$ppm <- massError
    missing <- is.na(flatMS3pair$XLink.AA.1) | is.na(flatMS3pair$XLink.AA.2)
    flatMS3pair <- flatMS3pair[!missing, ]
    if (nrow(flatMS3pair)==0) {
      placeHolder <- placeHolder[-index]
    } else {
      placeHolder[[index]] <- flatMS3pair
    }
    index = index + 1
  }
  for (CSM in ASpairs) {
    massError <- massErrors.AS[CSM]
    peptidePair <- arrayInd(CSM, dimensions)
    flatMS3pair <- flattenMS3pair(scanGroupTibble[peptidePair[,1],],
                                  scanGroupTibble[peptidePair[,2],],
                                  fraction, masterScanNo,
                                  frag1="Alkene", frag2="Sulfenic")
    flatMS3pair$ppm <- massError
    missing <- is.na(flatMS3pair$XLink.AA.1) | is.na(flatMS3pair$XLink.AA.2)
    flatMS3pair <- flatMS3pair[!missing, ]
    if (nrow(flatMS3pair)==0) {
      placeHolder <- placeHolder[-index]
    } else {
      placeHolder[[index]] <- flatMS3pair
    }
    index = index + 1
  }
  for (CSM in AApairs) {
    massError <- massErrors.AA[CSM]
    peptidePair <- arrayInd(CSM, dimensions)
    flatMS3pair <- flattenMS3pair(scanGroupTibble[peptidePair[,1],],
                                  scanGroupTibble[peptidePair[,2],],
                                  fraction, masterScanNo,
                                  frag1="Alkene", frag2="Alkene")
    flatMS3pair$ppm <- massError
    missing <- is.na(flatMS3pair$XLink.AA.1) | is.na(flatMS3pair$XLink.AA.2)
    flatMS3pair <- flatMS3pair[!missing, ]
    if (nrow(flatMS3pair)==0) {
      placeHolder <- placeHolder[-index]
    } else {
      placeHolder[[index]] <- flatMS3pair
    }
    index = index + 1
  }
  for (CSM in TTpairs) {
    massError <- massErrors.TT[CSM]
    peptidePair <- arrayInd(CSM, dimensions)
    flatMS3pair <- flattenMS3pair(scanGroupTibble[peptidePair[,1],],
                                  scanGroupTibble[peptidePair[,2],],
                                  fraction, masterScanNo,
                                  frag1="Thiol\\(Unsaturated\\)", frag2="Thiol\\(Unsaturated\\)")
    flatMS3pair$ppm <- massError
    missing <- is.na(flatMS3pair$XLink.AA.1) | is.na(flatMS3pair$XLink.AA.2)
    flatMS3pair <- flatMS3pair[!missing, ]
    if (nrow(flatMS3pair)==0) {
      placeHolder <- placeHolder[-index]
    } else {
      placeHolder[[index]] <- flatMS3pair
    }
    index = index + 1
  }
  if (length(placeHolder)==0) {
    return(data.frame())
  } else {
    return(do.call(rbind, placeHolder))
  }
}

#' Helper function called by `ms3ionFragFind`.  Combined two linked linear MS3
#' PSMs into a Prospector-Touchstone MS2 crosslinked data frame.
#' @param ms3CSM1 The first MS3-PSM.
#' @param ms3CSM2 The second MS3-PSM.
#' @param fraction Base name of RAW file.
#' @param masterScanNo Scan number of the triggering scan (typically MS2-CID).
#' @param frag1 Regex to detect the variable modification corresponding to the cleaved crosslinked fragment.
#' @param frag2 Regex to detect the variable modification corresponding to the other cleaved crosslinked fragment.
#' @seealso [ms3ionFragFind()], [processMS3xlinkResults()], [addMasterScanInfo()]
#' @returns A data frame
#' @export
flattenMS3pair <- function(ms3CSM1, ms3CSM2, fraction, masterScanNo=NA,
                           frag1="Alkene", frag2="Thiol\\(Unsaturated\\)") {
  if (ms3CSM1$Expect.ms3 > ms3CSM2$Expect.ms3) {
    tempCSM <- ms3CSM1
    ms3CSM1 <- ms3CSM2
    ms3CSM2 <- tempCSM
    tempFrag <- frag1
    frag1 <- frag2
    frag2 <- tempFrag
  }
  flatCSM <- tibble(MSMS.Info=masterScanNo,
                    MSMS.Info.1=ms3CSM1$`ms3ScanNo`,
                    MSMS.Info.2=ms3CSM2$`ms3ScanNo`,
                    DB.Peptide.1=ms3CSM1$`DB.Peptide.ms3`,
                    DB.Peptide.2=ms3CSM2$`DB.Peptide.ms3`,
                    Len.Pep.1=nchar(ms3CSM1$`DB.Peptide.ms3`),
                    Len.Pep.2=nchar(ms3CSM2$`DB.Peptide.ms3`),
                    Peptide.1=ms3CSM1$Peptide.ms3,
                    Peptide.2=ms3CSM2$Peptide.ms3,
                    Score=ms3CSM1$Score.ms3 + ms3CSM2$Score.ms3,
                    Score.Diff=NA,
                    Fraction=fraction,
                    RT.1=ms3CSM1$RT.ms3,
                    RT.2=ms3CSM2$RT.ms3,
                    Spectrum.1=ms3CSM1$Spectrum.ms3,
                    Spectrum.2=ms3CSM2$Spectrum.ms3,
                    ppm=NA,
                    Elemental.Composition=NA,
                    Sc.1=ms3CSM1$Score.ms3,
                    Sc.2=ms3CSM2$Score.ms3,
                    Expect.1=ms3CSM1$Expect.ms3,
                    Expect.2=ms3CSM2$Expect.ms3,
                    Protein.1=ms3CSM1$`Protein.ms3`,
                    Protein.2=ms3CSM2$`Protein.ms3`,
                    Acc.1=ms3CSM1$`Acc.ms3`,
                    Acc.2=ms3CSM2$`Acc.ms3`,
                    Decoy.1=ms3CSM1$Decoy.ms3,
                    Decoy.2=ms3CSM2$Decoy.ms3,
                    z=ms3CSM1$charge.ms2,
                    z.1=ms3CSM1$charge.ms3,
                    z.2=ms3CSM2$charge.ms3,
                    AA.1=stringr::str_match_all(ms3CSM1$`Protein.Mods.ms3`,
                                                stringr::str_c("XL:A-", frag1, "@([[0-9|]]+)"))[[1]][,2],
                    AA.2=stringr::str_match_all(ms3CSM2$`Protein.Mods.ms3`,
                                                stringr::str_c("XL:A-", frag2, "@([[0-9|]]+)"))[[1]][,2],
                    XLink.AA.1=stringr::str_c(.data$AA.1, sep="|"),
                    XLink.AA.2=stringr::str_c(.data$AA.2, sep="|")
  )
  return(flatCSM)
}
