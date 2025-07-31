#' Prepare PDB file for crosslinked distance measurements
#'
#' Searches for pdb file in the pdbFileDir and if it fails to find it there,
#' downloads from the pdb repository. Works with either pdb or mmCIF formatted
#' files. Retains only Calpha coordinates. Wrapper for functions in the `Bio3d` package.
#'
#' @param pdbID Code for pdb/cif file
#' @param pdbFileDir Path to local directory containin pdb/cif file
#' @return A data frame.
#' @export
#'
parsePDB<- function(pdbID=NA, pdbFileDir=getwd()) {
  if (missing(pdbID)) {
    stop("read.pdb: please specify a PDB 'file' for reading")
  }
  pdbSplit <- unlist(stringr::str_split(pdbID, stringr::fixed(".")))
  if (length(pdbSplit) > 1) {
    pdbExtension <- pdbSplit[length(pdbSplit)]
    if (!pdbExtension %in% c("cif", "pdb")) pdbExtension <- "pdb"
    pdbFront <- pdbSplit[-length(pdbSplit)]
  } else {
    pdbFront <- pdbID
    pdbExtension <- "pdb"
  }
  pdbFile <- stringr::str_c(pdbFront, pdbExtension, sep=".")
  if (file.exists(file.path(pdbFileDir, pdbFile))) {
    parsedPDB <- tryCatch(
      switch(pdbExtension,
             "cif" = bio3d::read.cif(file.path(pdbFileDir, pdbFile), verbose=F),
             "pdb" = bio3d::read.pdb(file.path(pdbFileDir, pdbFile), verbose=F)
      ),
      error = function(e) {
        message(stringr::str_c("Trouble reading pdb id: ", pdbCode, "\n", e))
        return(NA)
      }
    )
  } else {
    pdbCode = stringr::str_extract(pdbFront, "^[[0-9A-Za-z]]{4}")
    if (is.na(pdbCode)) {
      message(stringr::str_c(pdbID, " is not a valid pdb id"))
      return(NA)
    }
    parsedPDB <- tryCatch(
      switch(pdbExtension,
             "cif" = bio3d::read.cif(pdbCode, verbose=F),
             "pdb" = bio3d::read.pdb(pdbCode, verbose=F)
      ),
      error = function(e) {
        message(stringr::str_c("Trouble reading pdb id: ", pdbCode, "\n", e))
        return(NA)
      }
    )
  }
  if (sum(class(parsedPDB) %in% c("pdb", "cif")) == 0) return(NA)
  ca.ind <- bio3d::atom.select(parsedPDB, "calpha")
  parsedPDB <- parsedPDB$atom[ca.ind$atom,]
  return(parsedPDB)
}

#' Downloads and prepares all pdb files specified in the module file.
#'
#' @param modFile A parsed Touchstone module file.
#' @param pdbFileDir Path to local directory containing pdb/cif file
#' @return A list of Calpha coordinates for each pdb id.
#' @seealso [parsePDB()]
#' @export
#'
gatherPDBs <- function(modFile, pdbFileDir=getwd()) {
  pdbs <- modFile$PDB.code %>% unique
  pdbs <- pdbs[!is.na(pdbs)]
  structFiles <- purrr::map(pdbs, parsePDB, pdbFileDir)
  names(structFiles) <- pdbs
  return(structFiles)
}

#' Measure euclidean distance between Calpha atoms of a crosslinked residue pair
#'
#' @param XLink.AA.1 Sequence index of 1st residue.
#' @param PDB.chain.1 PDB identifier for 1st residue.
#' @param XLink.AA.2 Sequence index of 2nd residue.
#' @param PDB.chain.2 PDB identifier for 2nd residue.
#' @param parsedPDB A parsed Touchstone pdb file.
#' @seealso [parsePDB()], [euclideanDistanceHomologous()], [euclideanHelper()]
#' @return A numeric.
#' @export
#'
euclideanDistance <- function(XLink.AA.1, PDB.chain.1, XLink.AA.2, PDB.chain.2, parsedPDB) {
  if (is.na(PDB.chain.1) | is.na(PDB.chain.2)) {
    return(NA_real_)
  } else {
    coord1 <- parsedPDB[parsedPDB$chain==PDB.chain.1 & parsedPDB$resno==XLink.AA.1,
                        c("x","y","z")]
    if (nrow(coord1)==0) coord1 <- NA
    coord2 <- parsedPDB[parsedPDB$chain==PDB.chain.2 & parsedPDB$resno==XLink.AA.2,
                        c("x","y","z")]
    if (nrow(coord2)==0) coord2 <- NA
    distance <- sqrt(sum((coord1 - coord2)**2))
    return(round(distance,2))
  }
}

#' Measure crosslinked euclidean distances over the entire results file.
#'
#' @param datTab Parsed CLMS search results.
#' @param parsedPDB A parsed Touchstone pdb file.
#' @seealso [parsePDB()], [euclideanDistanceHomologous()]
#' @return A data frame.
#' @export
#'
euclideanHelper <- function(datTab, parsedPDB) {
  if (!is.null(parsedPDB)) {
    datTabDist <- datTab %>%
      #      mutate(distance = purrr::pmap_dbl(., euclideanDistanceHomologous, parsedPDB))
      mutate(distance = purrr::pmap_dbl(
        list(.data$XLink.AA.1, .data$PDB.chain.1, .data$XLink.AA.2, .data$PDB.chain.2),
        euclideanDistanceHomologous,
        parsedPDB)
      )
  } else {
    datTabDist <- datTab %>%
      mutate(distance = NA_real_)
  }
  return(datTabDist)
}

#' Measure euclidean distance between Calpha atoms of a crosslinked residue pair when there are homologous results.
#'
#' `euclideanDistanceHomologous()` expects that all possible chain ids are supplied as
#' a comma separated character string. For instance is the first crosslinked residue can
#' be either position 100 on either chain A, D, G, `XLink.AA.1` value should be `100` and
#' the `PDB.chain.1` should be given as `"A, D, G"`. All possible euclidean distances
#' are evaluated, and the shortest distance is reported. This function is meant to
#' handle cases in which a protein complex has multiple copies of the same subunit,
#' rather than situations where the crosslinked peptide appears in multiple different
#' homologous proteins.
#' @param XLink.AA.1 Sequence index of 1st residue.
#' @param PDB.chain.1 PDB identifier for 1st residue.
#' @param XLink.AA.2 Sequence index of 2nd residue.
#' @param PDB.chain.2 PDB identifier for 2nd residue.
#' @param parsedPDB A parsed Touchstone pdb file.
#' @seealso [parsePDB()], [euclideanDistance()], [euclideanHelper()]
#' @return A numeric.
#' @export
#'
euclideanDistanceHomologous <- function(XLink.AA.1, PDB.chain.1, XLink.AA.2, PDB.chain.2, parsedPDB) {
  PDB.chain.1 <- stringr::str_remove_all(PDB.chain.1, "\\s+")
  PDB.chain.2 <- stringr::str_remove_all(PDB.chain.2, "\\s+")
  if (is.na(PDB.chain.1) | is.na(PDB.chain.2)) {
    return(NA_real_)
  } else {
    distanceBox <- numeric()
    index <- 1
    for (chain1 in stringr::str_split(PDB.chain.1, ",")[[1]]) {
      for (chain2 in stringr::str_split(PDB.chain.2, ",")[[1]]) {
        coord1 <- parsedPDB[parsedPDB$chain==chain1 & parsedPDB$resno==XLink.AA.1,
                            c("x","y","z")]
        if (nrow(coord1)==0) coord1 <- NA
        coord2 <- parsedPDB[parsedPDB$chain==chain2 & parsedPDB$resno==XLink.AA.2,
                            c("x","y","z")]
        if (nrow(coord2)==0) coord2 <- NA
        distance <- sqrt(sum((coord1 - coord2)**2))
        distanceBox[index] <- distance
        index <- index + 1
      }
    }
    if (sum(is.na(distanceBox)) == length(distanceBox)) {
      return(NA)
    } else {
      return(round(min(distanceBox, na.rm=T), 2))
    }
  }
}

#' Sample random C-alpha distances on a PDB structure.
#'
#' Generates a distribution of random C-alpha, C-alpha distances on a PDB structure.
#' The size of the distibution is given my `numCrosslinks`. By default, only Lys-Lys
#' pairs are sampled, but this can be changed with the `resType` parameters.
#'
#' @param parsedPDB A parsed Touchstone pdb file.
#' @param numCrosslinks Number of random crosslinks generated.
#' @param resType1 A character vector containing three letter amino acid codes, from which the random crosslinks are selected.
#' @param resType2 A character vector containing three letter amino acid codes, from which the random crosslinks are selected.
#' @param chains If not `null`, restrict the random distribution to certain chains, given in a character vector.
#' @seealso [randomHelper()]
#' @return A numeric vector.
#' @export
#'
randomCrosslinks <- function(parsedPDB, numCrosslinks,
                             resType1="LYS", resType2="LYS",
                             chains=NULL) {
  if (!is.null(chains)) {
    parsedPDB <- parsedPDB %>% filter(.data$chain %in% chains)
  }
  pdb1 <- parsedPDB %>%
    filter(.data$resid %in% resType1) %>%
    dplyr::slice_sample(n=numCrosslinks, replace=T) %>%
    select(XLink.AA.1=.data$resno, PDB.chain.1=.data$chain)
  pdb2 <- parsedPDB %>%
    filter(.data$resid %in% resType2) %>%
    dplyr::slice_sample(n=numCrosslinks, replace=T) %>%
    select(XLink.AA.2=.data$resno, PDB.chain.2=.data$chain)
  pdb <- cbind(pdb1, pdb2)
  randomDists <- pdb %>%
    purrr::pmap_dbl(euclideanDistance, parsedPDB)
  return(randomDists)
}

#' Generate random distributions of crosslinked euclidean distances across a dataset.
#'
#' @param datTab Parsed CLMS search results.
#' @param parsedPDB A parsed Touchstone pdb file.
#' @param resType1 A character vector containing three letter amino acid codes, from which the random crosslinks are selected.
#' @param resType2 A character vector containing three letter amino acid codes, from which the random crosslinks are selected.
#' @param chains If not `null`, restrict the random distribution to certain chains, given in a character vector.
#' @seealso [randomCrosslinks()]
#' @return A data frame.
#' @export
#'
randomHelper <- function(datTab, parsedPDB, resType1="LYS", resType2="LYS", chains=NULL) {
  if (!is.null(parsedPDB)) {
    numCrosslinks <- nrow(datTab)
    datTab$random.distance <- randomCrosslinks(parsedPDB, numCrosslinks,
                                               resType1, resType2, chains=NULL)
  } else {
    datTab$random.distance <- NA_real_
  }
  return(datTab)
}

#' Measure crosslinked and randomized distributions of euclidean distances for a dataset.
#'
#' For each pdb id assigned in the module file, euclidean crosslinked and randomized
#' Calpha-Calpha distance distributions are measured and assigned to the dataset.
#'
#' @param datTabMod Parsed CLMS search results, with modules and chain ids assigned.
#' @param pdbList A list of pdb ids.
#' @seealso [processModuleFile()], [euclideanHelper()], [randomHelper()]
#' @return a data frame
#' @export
#'
measureCrosslinkDistance <- function(datTabMod, pdbList) {
  nestedDatTabMod <- datTabMod %>%
    group_by(.data$PDB) %>%
    nest(.key="data.nested")
  if (sum(is.na(nestedDatTabMod$PDB)) > 0) {
    nestedDatTabNA <- nestedDatTabMod %>%
      filter(is.na(.data$PDB)) %>%
      unnest(cols=c(.data$data.nested)) %>%
      ungroup() %>%
      mutate(random.distance = NA_real_)
  } else {
    nestedDatTabNA <- data.frame()
  }
  nestedDatTabMod <- nestedDatTabMod %>%
    filter(!is.na(.data$PDB)) %>%
    mutate(parsedPDB = pdbList[.data$PDB]) %>%
    mutate(data.nested=purrr::map2(.data$data.nested, .data$parsedPDB, euclideanHelper)) %>%
    mutate(data.nested=purrr::map2(.data$data.nested, .data$parsedPDB, randomHelper)) %>%
    select(-.data$parsedPDB) %>%
    unnest(cols=c(.data$data.nested)) %>%
    ungroup()
  bind_rows(nestedDatTabMod, nestedDatTabNA) %>%
    relocate(.data$PDB, .after=.data$PDB.chain.2)
}

#' Parse a Touchstone formatted module file.
#'
#' @param pathToModFile Location of module file.
#' @return A data frame.
#' @export
#'
readModuleFile <- function(pathToModFile) {
  modFile <- read_tsv(pathToModFile,
                      col_types = cols(
                        Acc = col_character(),
                        Protein.name = col_character(),
                        Module = col_character(),
                        PDB.code = col_character(),
                        PDB.chain = col_character(),
                        first.AA = col_integer(),
                        last.AA = col_integer()
                      )
  )
  if (sum(is.na(modFile$last.AA)) > 0) {

    # If the first/last AA columns are missing from mod file, create them.
    # get protein lengths from Uniprot... it might be more efficient, to force
    # Prospector to output sequence length in search compare and use instead.
    if (!"first.AA" %in% names(modFile)) {
      modFile$first.AA <- 0L
    }
    if (!"last.AA" %in% names(modFile)) {
      modFile$last.AA <- purrr::map_int(modFile$Acc, getUniprotSeqLength(modFile$Acc))
    }

    # If some entries don't have first/last NA specified, eg... there is some
    # domain of particular interest, but rest of the protein should default to some
    # other module, then first.AA needs to be NA.  Fill in last.AA to protein length
    # for purpose of displaying XiNet visualizations correctly.
    modFile <- modFile %>%
      mutate(last.AA = purrr::map2_int(.data$last.AA, .data$Acc, function(last.AA, Acc) {
        if(is.na(last.AA)) {
          return(getUniprotSeqLength(Acc))
        } else {
          return(as.integer(last.AA))
        }
      }))
    tryCatch(readr::write_tsv(modFile, pathToModFile),
             error = function(e) {
               message(stringr::str_c("Could not write module file: ", e))
               return(modFile)
             }
    )
  }
  modFile <- modFile %>%
    select(.data$Acc, .data$Protein.name, .data$Module, .data$first.AA, .data$last.AA,
           any_of(c("PDB.code", "PDB.chain")))
  return(modFile)
}

#' Returns the length of a protein sequence in number of AA residues.
#'
#' @param uniprotAcc A valid Uniprot Accession Number.
#' @returns an integer
#' @export
#'
getUniprotSeqLength <- function(uniprotAcc) {
  seqLen <- tryCatch(
    {seq <- bio3d::uniprot(uniprotAcc)$sequence
    proteinLen <- nchar(seq)
    },
    error = function(e) {
      message(stringr::str_c("trouble accessing Uniprot ID: ", uniprotAcc, "\n", e))
      return(50000L)
    }
  )
}

#' Assigns crosslinks to modules defined in a Touchstone formatted module file
#'
#' @param datTab Parsed CLMS search results.
#' @param moduleFile Parsed, Touchstone format module file
#' @returns a data frame
#' @export
#'
assignModules <- function(datTab, moduleFile) {
  datTab <- datTab %>%
    select(-any_of(c(starts_with("xlinkedModulPair"), starts_with("Module."), starts_with("PDB.code."), starts_with("PDB.chain."))))
  mods <- unique(moduleFile$Module) %>% stringr::str_sort()
  accs <- sort(unique(as.character(c(datTab$Acc.1, datTab$Acc.2))))
  datTab <- datTab %>% group_by(.data$Decoy2)
  datTabList <- datTab %>% group_split()
  names(datTabList) <- datTab %>% group_keys() %>% pull(.data$Decoy2)
  datTab.target <- datTabList$Target %>%
    left_join(moduleFile, by=c("Acc.1"="Acc"), relationship = "many-to-many") %>%
    mutate("keepEntry" = suppressWarnings(purrr::pmap_lgl(select(., .data$XLink.AA.1, .data$first.AA, .data$last.AA),
                                                         function(XLink.AA.1, first.AA, last.AA)
                                                           dplyr::between(XLink.AA.1, first.AA, last.AA))
    )) %>%
    filter(.data$keepEntry | is.na(.data$keepEntry)) %>%
    group_by(.data$Fraction, .data$MSMS.Info, .data$xlinkedResPair) %>%
    add_count(name="redundance") %>%
    filter(.data$redundance == 1 | (.data$redundance > 1 & .data$keepEntry)) %>%
    ungroup() %>%
    select(-.data$first.AA, -.data$last.AA, -.data$keepEntry, -.data$Protein.name, -.data$redundance) %>%
    left_join(moduleFile, by=c("Acc.2"="Acc"), suffix=c(".1", ".2"), relationship = "many-to-many") %>%
    mutate("keepEntry" = suppressWarnings(purrr::pmap_lgl(select(., .data$XLink.AA.2, .data$first.AA, .data$last.AA),
                                                         function(XLink.AA.2, first.AA, last.AA)
                                                           dplyr::between(XLink.AA.2, first.AA, last.AA))
    )) %>%
    filter(.data$keepEntry | is.na(.data$keepEntry)) %>%
    group_by(.data$Fraction, .data$MSMS.Info, .data$xlinkedResPair) %>%
    add_count(name="redundance") %>%
    filter(.data$redundance == 1 | (.data$redundance > 1 & .data$keepEntry)) %>%
    ungroup() %>%
    select(-.data$first.AA, -.data$last.AA, -.data$keepEntry, -.data$Protein.name, -.data$redundance) %>%
    mutate(xlinkedModulPair = ifelse(.data$Module.1 <= .data$Module.2,
                                     stringr::str_c(.data$Module.1, .data$Module.2, sep="::"),
                                     stringr::str_c(.data$Module.2, .data$Module.1, sep="::"))
    ) %>%
    mutate(xlinkedModulPair = forcats::as_factor(.data$xlinkedModulPair),
           Module.1 = factor(.data$Module.1, levels = mods),
           Module.2 = factor(.data$Module.2, levels = mods),
           Acc.1 = factor(.data$Acc.1, levels = accs),
           Acc.2 = factor(.data$Acc.2, levels = accs)) #%>%
  #         add_count(xlinkedModulPair, name="numMPSM")
  if (sum(stringr::str_detect(names(datTab.target), "PDB.code")) > 0) {
    datTab.target <- datTab.target %>%
      mutate(PDB = ifelse(!is.na(.data$PDB.code.1) & !is.na(.data$PDB.code.2) & .data$PDB.code.1 == .data$PDB.code.2,
                          .data$PDB.code.1, NA)
      )
  }
  if (length(datTabList) == 2) {
    datTab.decoy <- datTabList$Decoy
    decoyLength <- nrow(datTab.decoy)
    colsToAdd <- setdiff(names(datTab.target), names(datTab.decoy))
    emptyCols <- matrix(nrow=decoyLength, ncol=length(colsToAdd))
    colnames(emptyCols) <- colsToAdd
    emptyCols <- as_tibble(emptyCols)
    datTab.decoy <- bind_cols(datTab.decoy, emptyCols)
    datTab <- bind_rows(datTab.target, datTab.decoy)
  } else {
    datTab <- datTab.target
  }
  return(datTab)
}

#' Wrapper function that assigns modules and distances defined in a Touchstone Module File
#'
#' @param datTab Parsed CLMS search results.
#' @param pathToModFile Location of module file
#' @param pdbFileDir Path to local directory containing pdb/cif file
#' @returns a data frame
#' @export
#'
processModuleFile <- function(datTab, pathToModFile, pdbFileDir = getwd()) {
  modFile <- readModuleFile(pathToModFile)
  datTab <- assignModules(datTab, modFile)
  if ("PDB.code" %in% names(modFile) & sum(!is.na(datTab$PDB) > 0)) {
    pdbFiles <- gatherPDBs(modFile, pdbFileDir = pdbFileDir)
    datTab <- measureCrosslinkDistance(datTab, pdbFiles)
  }
  return(datTab)
}
