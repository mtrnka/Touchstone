parsePDB <- function(pdbID=NA, pdbFileDir=getwd()) {
  if (is.na(pdbID)) {
    message("pdb id missing")
    return(NA)
  }
  pdbSplit <- unlist(str_split(pdbID, fixed(".")))
  if (length(pdbSplit) > 1) {
    pdbExtension <- pdbSplit[length(pdbSplit)]
    if (!pdbExtension %in% c("cif", "pdb")) pdbExtension <- "pdb"
    pdbFront <- pdbSplit[-length(pdbSplit)]
  } else {
    pdbFront <- pdbID
    pdbExtension <- "pdb"
  }
  pdbFile <- str_c(pdbFront, pdbExtension, sep=".")
  if (file.exists(file.path(pdbFileDir, pdbFile))) {
    pdbCode <- file.path(pdbFileDir, pdbFile)
  } else {
    pdbCode = str_extract(pdbFront, "^[[0-9A-Za-z]]{4}")
  }
  if (is.na(pdbCode)) {
    message(str_c(pdbID, " is not a valid pdb id"))
    return(NA)
  }
  parsedPDB <- tryCatch(
    switch(pdbExtension,
           "cif" = read_cif(pdbCode, verbose=F, cacheDir=pdbFileDir),
           "pdb" = read_pdb(pdbCode, verbose=F, cacheDir=pdbFileDir)
    ),
    error = function(e) {
      message(str_c("Trouble reading pdb id: ", pdbCode, "\n", e))
      return(NA)
    }
  )
  if (sum(class(parsedPDB) %in% c("pdb", "cif")) == 0) return(NA)
  # if (!file.exists(file.path(pdbFileDir, pdbFile)) & ("pdb" %in% class(parsedPDB))) {
  #     write.pdb(parsedPDB, file.path(pdbFileDir, str_c(pdbCode, ".pdb")))
  # }
  ca.ind <- bio3d::atom.select(parsedPDB, "calpha")
  parsedPDB <- parsedPDB$atom[ca.ind$atom,]
  return(parsedPDB)
}

gatherPDBs <- function(modFile, pdbFileDir=getwd()) {
  pdbs <- modFile %>% pluck("PDB.code") %>% unique
  pdbs <- pdbs[!is.na(pdbs)]
  structFiles <- map(pdbs, parsePDB, pdbFileDir)
  names(structFiles) <- pdbs
  return(structFiles)
}

euclideanDistance <- function(XLink.AA.1, PDB.chain.1, XLink.AA.2, PDB.chain.2, parsedPDB, ...) {
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

euclideanHelper <- function(datTab, parsedPDB) {
  if (!is.null(parsedPDB)) {
    datTabDist <- datTab %>%
      mutate(distance = pmap_dbl(., euclideanDistanceMultimer, parsedPDB))
  } else {
    datTabDist <- datTab %>%
      mutate(distance = NA_real_)
  }
  return(datTabDist)
}

euclideanDistanceMultimer <- function(XLink.AA.1, PDB.chain.1, XLink.AA.2, PDB.chain.2, parsedPDB, ...) {
  PDB.chain.1 <- str_remove_all(PDB.chain.1, "\\s+")
  PDB.chain.2 <- str_remove_all(PDB.chain.2, "\\s+")
  if (is.na(PDB.chain.1) | is.na(PDB.chain.2)) {
    return(NA_real_)
  } else {
    distanceBox <- numeric()
    index <- 1
    for (chain1 in str_split(PDB.chain.1, ",")[[1]]) {
      for (chain2 in str_split(PDB.chain.2, ",")[[1]]) {
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

randomCrosslinks <- function(parsedPDB, numCrosslinks,
                             resType1="LYS", resType2="LYS",
                             chains=NULL) {
  if (!is.null(chains)) {
    parsedPDB <- parsedPDB %>% filter(chain %in% chains)
  }
  pdb1 <- parsedPDB %>%
    filter(resid %in% resType1) %>%
    sample_n(numCrosslinks, replace=T) %>%
    select(XLink.AA.1=resno, PDB.chain.1=chain)
  pdb2 <- parsedPDB %>%
    filter(resid %in% resType2) %>%
    sample_n(numCrosslinks, replace=T) %>%
    select(XLink.AA.2=resno, PDB.chain.2=chain)
  pdb <- cbind(pdb1, pdb2)
  randomDists <- pdb %>%
    pmap_dbl(euclideanDistance, parsedPDB)
  return(randomDists)
}

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

measureCrosslinkDistance <- function(datTabMod, pdbList) {
  nestedDatTabMod <- datTabMod %>%
    group_by(PDB) %>%
    nest()
  if (sum(is.na(nestedDatTabMod$PDB)) > 0) {
    nestedDatTabNA <- nestedDatTabMod %>%
      filter(is.na(PDB)) %>%
      unnest(col=c(data, PDB)) %>%
      ungroup() %>%
      mutate(random.distance = NA_real_)
  } else {
    nestedDatTabNA <- data.frame()
  }
  nestedDatTabMod <- nestedDatTabMod %>%
    filter(!is.na(PDB)) %>%
    mutate(parsedPDB = pdbList[PDB]) %>%
    mutate(data=map2(data, parsedPDB, euclideanHelper)) %>%
    mutate(data=map2(data, parsedPDB, randomHelper)) %>%
    select(-parsedPDB) %>%
    unnest(col=c(data, PDB)) %>%
    ungroup()
  bind_rows(nestedDatTabMod, nestedDatTabNA) %>%
    relocate(PDB, .after=PDB.chain.2)
}

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
      modFile$last.AA <- map_int(modFile$Acc, getUniprotSeqLength(Acc))
    }

    # If some entries don't have first/last NA specified, eg... there is some
    # domain of particular interest, but rest of the protein should default to some
    # other module, then first.AA needs to be NA.  Fill in last.AA to protein length
    # for purpose of displaying XiNet visualizations correctly.
    modFile <- modFile %>%
      mutate(last.AA = map2_int(last.AA, Acc, function(last.AA, Acc) {
        if(is.na(last.AA)) {
          return(getUniprotSeqLength(Acc))
        } else {
          return(as.integer(last.AA))
        }
      }))
    tryCatch(write_tsv(modFile, pathToModFile),
             error = function(e) {
               message(str_c("Could not write module file: ", e))
               return(modFile)
             }
    )
  }
  modFile <- modFile %>%
    select(Acc, Protein.name, Module, first.AA, last.AA,
           any_of(c("PDB.code", "PDB.chain")))
  return(modFile)
}

getUniprotSeqLength <- function(uniprotAcc) {
  seqLen <- tryCatch(
    {seq <- bio3d::uniprot(uniprotAcc)$sequence
    proteinLen <- nchar(seq)
    },
    error = function(e) {
      message(str_c("trouble accessing Uniprot ID: ", uniprotAcc, "\n", e))
      return(50000L)
    }
  )
}

assignModules <- function(datTab, moduleFile) {
  mods <- unique(moduleFile$Module) %>% str_sort
  accs <- levels(datTab$Acc.1)
  datTab <- datTab %>% group_by(Decoy2)
  datTabList <- datTab %>% group_split()
  names(datTabList) <- datTab %>% group_keys() %>% pull(Decoy2)
  datTab.target <- datTabList$Target %>%
    left_join(moduleFile, by=c("Acc.1"="Acc"), relationship = "many-to-many") %>%
    mutate("keepEntry" = suppressWarnings(pmap_lgl(select(., XLink.AA.1, first.AA, last.AA),
                                                   function(XLink.AA.1, first.AA, last.AA)
                                                     dplyr::between(XLink.AA.1, first.AA, last.AA))
    )) %>%
    filter(keepEntry | is.na(keepEntry)) %>%
    group_by(Fraction, MSMS.Info, xlinkedResPair) %>%
    add_count(name="redundance") %>%
    filter(redundance == 1 | (redundance > 1 & keepEntry)) %>%
    ungroup() %>%
    select(-first.AA, -last.AA, -keepEntry, -Protein.name, -redundance) %>%
    left_join(moduleFile, by=c("Acc.2"="Acc"), suffix=c(".1", ".2"), relationship = "many-to-many") %>%
    mutate("keepEntry" = suppressWarnings(pmap_lgl(select(., XLink.AA.2, first.AA, last.AA),
                                                   function(XLink.AA.2, first.AA, last.AA)
                                                     dplyr::between(XLink.AA.2, first.AA, last.AA))
    )) %>%
    filter(keepEntry | is.na(keepEntry)) %>%
    group_by(Fraction, MSMS.Info, xlinkedResPair) %>%
    add_count(name="redundance") %>%
    filter(redundance == 1 | (redundance > 1 & keepEntry)) %>%
    ungroup() %>%
    select(-first.AA, -last.AA, -keepEntry, -Protein.name, -redundance) %>%
    mutate(xlinkedModulPair = ifelse(Module.1 <= Module.2,
                                     str_c(Module.1, Module.2, sep="::"),
                                     str_c(Module.2, Module.1, sep="::"))
    ) %>%
    mutate(xlinkedModulPair = as_factor(xlinkedModulPair),
           Module.1 = factor(Module.1, levels = mods),
           Module.2 = factor(Module.2, levels = mods),
           Acc.1 = factor(Acc.1, levels = accs),
           Acc.2 = factor(Acc.2, levels = accs)) #%>%
  #         add_count(xlinkedModulPair, name="numMPSM")
  if (sum(str_detect(names(datTab.target), "PDB.code")) > 0) {
    datTab.target <- datTab.target %>%
      mutate(PDB = ifelse(!is.na(PDB.code.1) & !is.na(PDB.code.2) & PDB.code.1 == PDB.code.2,
                          PDB.code.1, NA)
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

processModuleFile <- function(datTab, pathToModFile, pdbFileDir = getwd()) {
  modFile <- readModuleFile(pathToModFile)
  datTab <- assignModules(datTab, modFile)
  if ("PDB.code" %in% names(modFile) & sum(!is.na(datTab$PDB) > 0)) {
    pdbFiles <- gatherPDBs(modFile, pdbFileDir = pdbFileDir)
    datTab <- measureCrosslinkDistance(datTab, pdbFiles)
  }
  return(datTab)
}
