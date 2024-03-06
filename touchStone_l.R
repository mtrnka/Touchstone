setClass(Class="PPsearchCompareXL",
         slots=c(
             parentDir="character",
             dataTable="data.frame",
             peaklistDir="character",
             fastaDB="character",
             modulFile="list",
             chainMap="list",
             caTracedPDB="data.frame",
             reductionState="factor",
             expSource="character")
)

setMethod("initialize", "PPsearchCompareXL",
          function(.Object, directory, dataFile, modFile, 
                   pdbDir, pdbFile, chainMapFile, 
                   redState=factor("spectralMatchPairs"), 
                   expSource=NA_character_){
              .Object@expSource <- expSource
              levels(.Object@reductionState) <- c("spectralMatchPairs",
                                                  "xlinkedResPairs",
                                                  "xlinkedPepPairs")
              .Object@reductionState = redState
              cat("*** Reading Prospector XL Search Table *** \n")
              .Object@parentDir <- directory
              setwd(directory)
              .Object@dataTable <- readProspectorXLOutput(dataFile)
              cat("*** Reading Modules *** \n")
              .Object@modulFile <- readModuleFile(modFile)
              .Object@dataTable <- populateModules(.Object@dataTable,
                                                   .Object@modulFile)
              cat("*** Calculating Decoys *** \n")
              .Object@dataTable <- calculateDecoys(.Object@dataTable)
              cat("*** Calculating Xlinked Pairs *** \n")
              .Object@dataTable <- calculatePairs(.Object@dataTable)
              cat("*** Calculating Distances *** \n")
              chainFile <- paste(pdbDir,"/",chainMapFile,sep="")
              pdbFile <- paste(pdbDir,"/",pdbFile,sep="")
              .Object@chainMap <- readChainMap(chainFile)
              .Object@caTracedPDB <- parsePDB(pdbFile)
              .Object@dataTable <- measureDistances(.Object@dataTable,
                                                    .Object@caTracedPDB,
                                                    .Object@chainMap)
              #Add validity check on .Object
              #validObject(.Object)
              return(.Object)
          }
)

setMethod("show", "PPsearchCompareXL",
          function(object){
              cat("*** Class: PPsearchCompareXL *** \n")
              cat("*** Truncated Data table *** \n")
              print(head(object@dataTable, 5))
          }
)

setGeneric("getSearchTable", function(object) {
    standardGeneric("getSearchTable")
    }
)

setMethod("getSearchTable", "PPsearchCompareXL",
          function(object){
              return(object@dataTable)
          }
)

setGeneric("bestResPair", function(object, ...) {
    standardGeneric("bestResPair")
    }
)

setMethod("bestResPair", "PPsearchCompareXL",
          function(object, classifier=c("Score.Diff","Score",
                                        "Peptide.2","Peptide.1",
                                        "MSMS.Info")) {
              datTab <- object@dataTable
              datTabR <- object@dataTable
              xlinks <- unique(datTab$xlinkedResPair)
              for (param in classifier) {
                  datTabR <- reduceTo(datTabR, "xlinkedResPair", param)
              }
              cat(c(nrow(datTabR), length(xlinks),
                    sum(datTabR$num), nrow(datTab), "\n"), sep="\t")
              object@reductionState = factor("xlinkedResPairs")
              object@dataTable <- datTabR
              return(object)
          }
)

setGeneric("bestPepPair", function(object, ...) {
    standardGeneric("bestPepPair")
    }
)

setMethod("bestPepPair", "PPsearchCompareXL",
          function(object, classifier=c("Score.Diff","Score",
                                        "Peptide.2","Peptide.1",
                                        "MSMS.Info")) {
              datTab <- object@dataTable
              datTabR <- object@dataTable
              xlinks <- unique(datTab$xlinkedPepPair)
              for (param in classifier) {
                  datTabR <- reduceTo(datTabR, "xlinkedPepPair", param)
              }
              cat(c(nrow(datTabR), length(xlinks),
                    sum(datTabR$num), nrow(datTab), "\n"), sep="\t")
              return(datTabR)
          }
)

reduceTo <- function(datTab, category, classifier) {
    if (class(datTab[[classifier]]) == "factor") {
        datTab[[classifier]] = as.character(datTab[[classifier]])
    }
    bestScores <- aggregate(datTab[[classifier]] ~ datTab[[category]], data=datTab, max)
    names(bestScores) <- c(category, classifier)
    outTab <- merge(bestScores, datTab, by=c(category, classifier))
    if (!("num" %in% colnames(outTab))) {
        xlinks <- table(datTab[[category]])
        outTab$num <- sapply(outTab[[category]], function(x) xlinks[x])
    }
    return(outTab)
}

readProspectorXLOutput <- function(inputFile){
    inFile <- file(inputFile,open="r")
    header <- readLines(inFile,n=1)
    dataTable <- readLines(inFile,n=-1)
    close(inFile)
    #Clean up header:
    header <- lineSplit(header)
    header <- gsub("[[:space:]]",".",header)
    acc_pos <- grep("Acc",header)
    header[acc_pos[1]] <- "Acc.1"
    header[acc_pos[2]] <- "Acc.2"
    spec_pos <- grep("Species",header)
    header[spec_pos[1]] <- "Species.1"
    header[spec_pos[2]] <- "Species.2"
    prot_pos <- grep("Protein",header)
    header[prot_pos[1]] <- "Protein.1"
    header[prot_pos[2]] <- "Protein.2"
    #Clean up data:
    dataTable <- rbind(sapply(as.list(dataTable),lineSplit))
    dataTable <- as.data.frame(t(dataTable),stringsAsFactors=F)
    names(dataTable) <- header
    return(cleanTypes(dataTable))
}

readModuleFile <- function(modFile) {
    inFile <- file(modFile,open="r")
    header <- readLines(inFile,n=1)
    header <- lineSplit(header)
    dataTable <- character(4)
    while (T) {
        line <- readLines(inFile, n=1)
        line <- lineSplit(line)
        if (length(line)==0) {break}
        else if (length(line) < 4) {length(line)=4}
        dataTable <- rbind(dataTable,line)}
    close(inFile)
    row.names(dataTable) <- NULL
    dataTable <- as.data.frame(dataTable[-1,],stringsAsFactors=F)
    names(dataTable) <- header
    dataTable <- cleanTypes(dataTable)
    modules <- list()
    for (subunit in unique(dataTable$Subunit)) {
        conv <- dataTable[dataTable$Subunit==subunit,
                          c("Range_low","Range_high","Module")]
        modules[[subunit]] <- conv[order(conv$Range_low),]
    }
    return(modules)
}

readChainMap <- function(chainFile) {
    chainMap <- read.table(chainFile,header=F,sep="\t",stringsAsFactors=F)
    names(chainMap) <- c("Subunit","Chain")
    chains <- list()
    for (su in unique(chainMap$Subunit)) {
        chains[[su]]=unlist(chainMap[chainMap$Subunit==su,"Chain"])
    }
    return(chains)
}

parsePDB <- function(pdbAcc) {
    require(bio3d)
    struct <- read.pdb(pdbAcc)
    ca.ind <- atom.select(struct,"calpha")
    struct <- struct$atom[ca.ind$atom,]
    return(struct)
}

euclideanDistance <- function(parsedPDB, residue1, chain1, residue2, chain2) {
    coord1 <- parsedPDB[parsedPDB$chain==chain1 & parsedPDB$resno==residue1,
                        c("x","y","z")]
    if (nrow(coord1)==0) coord1 <- NA
    coord2 <- parsedPDB[parsedPDB$chain==chain2 & parsedPDB$resno==residue2,
                        c("x","y","z")]
    if (nrow(coord2)==0) coord2 <- NA
    distance <- sqrt(sum((coord1 - coord2)**2))
    return(round(distance,2))
}

multiEuclideanDistance <- function(parsedPDB, residue1, chain1, residue2, chain2) {
    minDist <- NA
    minChainPair <- NA
    for (chains1 in chain1) { for (chains2 in chain2) {
        dist <- euclideanDistance(parsedPDB, residue1, chains1, residue2, chains2)
        if (is.na(minDist) | dist < minDist) {
            minDist <- dist
            if (chains1 <= chains2) {
                minChainPair <- paste(chains1, chains2, sep="")
            } else {
                minChainPair <- paste(chains2, chains1, sep="")
            }
        }
    }}
    return(minDist)
    #    return(list(minDist, minChainPair))
}

measureDistances <- function(searchTable, parsedPDB, chainMap) {
    chainLookup <- function(proteinName) {chainMap[[proteinName]]}
    pdbList <- list(parsedPDB)
    chains1 <- sapply(searchTable$Protein.1, chainLookup)
    chains2 <- sapply(searchTable$Protein.2, chainLookup)
    distance <- mapply(
        multiEuclideanDistance,
        pdbList,
        searchTable$XLink.AA.1,
        chains1,
        searchTable$XLink.AA.2,
        chains2)
    searchTable$distance <- distance
    return(searchTable)
}

populateModules <- function(searchTable, moduleDefinitions) {
    modules <- list(moduleDefinitions)
    pep1.modul <- mapply(
        assignModule,
        searchTable$XLink.AA.1,
        searchTable$Protein.1,
        modules)
    pep2.modul <- mapply(
        assignModule,
        searchTable$XLink.AA.2,
        searchTable$Protein.2,
        modules)
    searchTable$Modul.1 <- as.factor(unlist(pep1.modul))
    searchTable$Modul.2 <- as.factor(unlist(pep2.modul))
    return(searchTable)
}

assignModule <- function(seqPosition,protein,modList) {
    if (protein %in% names(modList)) {
        n <- findInterval(as.integer(seqPosition),
                          modList[[as.character(protein)]]$Range_low)
        if (n != 0) {
            return(modList[[protein]]$Module[n])
        } else return(NA)
    }
    else return(NA)
}

calculateDecoys <- function(searchTable) {
    searchTable$Decoy <- "Target"
    decoyReg <- "(^r[[:digit:]]_|^dec)"
    searchTable[grepl(decoyReg,searchTable$Protein.1) |
                    grepl(decoyReg,searchTable$Protein.2),
                "Decoy"] <- "Decoy"
    searchTable[grepl(decoyReg,searchTable$Protein.1) &
                    grepl(decoyReg,searchTable$Protein.2),
                "Decoy"] <- "DoubleDecoy"
    searchTable$Decoy <- as.factor(searchTable$Decoy)
    return(searchTable)
}

calculatePairs <- function(searchTable){
    searchTable$Res.1 <- paste(searchTable$XLink.AA.1, searchTable$Acc.1, sep=".")
    searchTable$Res.2 <- paste(searchTable$XLink.AA.2, searchTable$Acc.2, sep=".")
    searchTable$xlinkedResPair <- ifelse(searchTable$Res.1 <= searchTable$Res.2, 
                                         paste(searchTable$Res.1, searchTable$Res.2, sep="::"),
                                         paste(searchTable$Res.2, searchTable$Res.1, sep="::"))
    searchTable$xlinkedResPair <- as.factor(searchTable$xlinkedResPair)
    searchTable$xlinkedPepPair <- ifelse(searchTable$DB.Peptide.1 <= searchTable$DB.Peptide.2,
                                         paste(searchTable$DB.Peptide.1, searchTable$DB.Peptide.2, sep="::"),
                                         paste(searchTable$DB.Peptide.2, searchTable$DB.Peptide.1, sep="::"))
    searchTable$xlinkedPepPair <- as.factor(searchTable$xlinkedPepPair)
    return(searchTable)
}

lineSplit <- function(line) {
    return(unlist(strsplit(line,"\t")))
}

cleanTypes <- function(dataTable) {
    return(as.data.frame(
        lapply(dataTable, function(x) {
            type.convert(x, as.is=TRUE)}),
        stringsAsFactors=FALSE))
}

calculateFDR <- function(datTab, scalingFactor=10) {
    fdrTable <- table(datTab$Decoy)
    fdr <- (fdrTable[["Decoy"]] - fdrTable[["DoubleDecoy"]]) /
        (scalingFactor * fdrTable[["Target"]])
    #print(fdrTable)
    #return(sprintf("FDR: %0.3f", fdr))
    return(fdr)
}

setGeneric("prepPairs", function(object, threshold) {
    standardGeneric("prepPairs")
    }
)

setMethod("prepPairs", "PPsearchCompareXL",
          function(object, threshold){
              object@dataTable <- thresholdResults(object@dataTable, threshold)
              object@dataTable <- removeDecoys(object@dataTable)
              return(object)
          }
)

removeDecoys <- function(datTab) {
    datTabR <- datTab[datTab$Decoy=="Target",]
    return(datTabR)
}

thresholdResults <- function(datTab, threshold, classifier="Score.Diff") {
    datTabR <- datTab[datTab[[classifier]] >= threshold,]
    return(datTabR)
}

setGeneric("splitSources", function(object, sourceTab) {
    standardGeneric("splitSources")
    }
)

setMethod("splitSources", "PPsearchCompareXL",
#sourceTab is a list that maps Fractions onto expSources
#st <- list("source1" = c("Fraction1, Fraction2, ...), "source2" = c(Fract))
          function(object, sourceTab){
              datTab <- object@dataTable 
              splitResult <- lapply(sourceTab, function (x) {
                  datTab2 <- datTab[datTab$Fraction %in% x, ]
                  obj <- object
                  obj@dataTable <- datTab2
                  obj }
                  )
              for (result in names(splitResult)) {
                  splitResult[[result]]@expSource <- result
              }
              return(splitResult)
          }
)

setGeneric("renameProtein", function(object, oldName, newName) {
    standardGeneric("renameProtein")
    }
)

setMethod("renameProtein", "PPsearchCompareXL",
          function(object, oldName, newName){
              object@dataTable <- 
                  .renameProtein(object@dataTable, oldName, newName)
              object@modulFile <- 
                  .renameModTable(object@modulFile, oldName, newName)
              return(object)
          }
)

.renameProtein <- function(datTab, oldName, newName) {
    datTab$Acc.1 <- gsub(oldName, newName, datTab$Acc.1)
    datTab$Acc.2 <- gsub(oldName, newName, datTab$Acc.2)
    datTab <- calculatePairs(datTab)
    return(datTab)
}

.renameModTable <- function(modTab, oldName, newName) {
    names(modTab) <- gsub(oldName, newName, names(modTab))
    return(modTab)
}

setGeneric("compareTwo", function(data1, data2) {
    standardGeneric("compareTwo")
    }
)

setMethod("compareTwo", "PPsearchCompareXL",
          function(data1, data2){
              if ((!data1@reductionState %in% c("xlinkedResPairs")) & 
                  (!data2@reductionState %in% c("xlinkedResPairs")) &
                  (!data2@reductionState == data1@reductionState)) {
               
             #Check that reduction state of two dataTables is the same...
              #Maybe several different methods for doing the comparison:
              dt2 <- data2@dataTable
              dt1 <- data1@dataTable
              }
          }
)


setGeneric("pairPlot", function(object, color="lightseagreen", 
                                removeMods=NA_character_, displayEmpty=T) {
    standardGeneric("pairPlot")
    }
)

setMethod("pairPlot", "PPsearchCompareXL",
          function(object, color, removeMods, displayEmpty){
              .pairPlot(object@dataTable, object@modulFile, color, 
                        removeMods, displayEmpty)
          }
)

.rejiggerMods <- function(modTab) {
    #modTab <- lapply(modTab, function(x) x[-nrow(x),])
    dummy <- function(x) {
        first <- x[1,]
        first[,2] <- first[,1]
        return(rbind(first, x))
    }
    modTab <- lapply(modTab, dummy)
    modTab <- do.call(rbind, modTab)
    modTab$Acc.1 <- gsub("\\..+$","", row.names(modTab))
    row.names(modTab) <- NULL
    names(modTab) <- c("Range_low","XLink.AA.1","Module","Acc.1")
    modTab <- modTab[,c(2,4)]
    modIndex <- expand.grid(1:nrow(modTab),1:nrow(modTab))
    modTab <- cbind(modTab[modIndex$Var1,], modTab[modIndex$Var2,])
    names(modTab) <- c("XLink.AA.1","Acc.1","XLink.AA.2","Acc.2")
    return(modTab)
}

.pairPlot <- function(datTab, modTab, col="lightseagreen", 
                      removeMods=NA_character_, displayEmpty=T) {
    #Function takes a dataTable which is reduced to Residue Pairs and plot
    #Dot for each unique crosslink with size equal to num of instances.  The
    #color can reflect the source of the crosslink to facilitate comparisons.
    #Basic skeleton of the function that will get built up later.  Ultimately
    #would like this to work on the XLSearchOutput object, which probably needs
    #flag to show if it is residue pair reduced etc.
    require(ggplot2)
    datTab <- removeModule(datTab, removeMods)
    datTabR <- datTab
    tempXLinkAA1 <- datTab$XLink.AA.1
    tempAcc1 <- datTab$Acc.1
    datTabR$XLink.AA.1 <- datTab$XLink.AA.2
    datTabR$Acc.1 <- datTab$Acc.2
    datTabR$XLink.AA.2 <- tempXLinkAA1
    datTabR$Acc.2 <- tempAcc1
    datTab <- rbind(datTab, datTabR)
    if (!displayEmpty) {
        modTab <- modTab[names(modTab) %in% 
                             unique(c(datTab$Acc.1, datTab$Acc.2))]
    }
    modTab <- .rejiggerMods(modTab)
    gg <- ggplot(datTab, aes(XLink.AA.1, XLink.AA.2)) +
        geom_point(aes(size=num), col=col, alpha=0.6) +
        geom_vline(aes(xintercept=XLink.AA.1), modTab, color="palevioletred2",
                   linetype="21", size=0.75, alpha=0.5) +
        geom_hline(aes(yintercept=XLink.AA.2), modTab, color="palevioletred2",
                   linetype="21", size=0.75, alpha=0.5) +
        facet_grid(Acc.2 ~ Acc.1, space="free", scales="free", switch="both") + 
        theme(panel.background = element_rect(fill=rgb(0.98,0.98,0.98)),
               panel.spacing = unit(0.5, units="line"),
               panel.border = element_rect(color = "grey", fill = NA, size = 1),
               panel.grid.major = element_line(color = "grey", size=.5),
               panel.grid.minor = element_line(color = "grey", size=.25)) +
        scale_x_continuous(breaks=seq(0,max(datTab$XLink.AA.1),by=50),
                           minor_breaks=seq(0,max(datTab$XLink.AA.1),by=25)) +
        scale_y_continuous(breaks=seq(0,max(datTab$XLink.AA.2),by=50),
                           minor_breaks=seq(0,max(datTab$XLink.AA.2),by=25))
    plot(gg)
}

removeModule <- function(datTab, modules) {
    datTab <- datTab[!(datTab$Modul.1 %in% modules | datTab$Modul.2 %in% modules),]
    return(datTab)
}



###############################################################################
#                               testing area                                  #
###############################################################################

# test <- new(Class="PPsearchCompareXL",
#            directory="/Users/mtrnka/Projects-Collaboration/GrossLab/Serena/edc_new/",
#            dataFile="edc_dda_res3.txt",
#            modFile="../module_definitions.txt",
#            pdbDir="../pdb",
#            pdbFile="1kx5.pdb",
#            chainMapFile="chain_map.txt")
# 
# test <- renameProtein(test, "TMLA-H3K9", "H3")
# 
# sourceTable <- list("nucleosomes"=
#                         levels(factor(test@dataTable$Fraction))[1:6],
#                     "nucleosomes+swi6"=
#                         levels(factor(test@dataTable$Fraction))[7:12])
# 
# test2 <- splitSources(test, sourceTable)
# 
# test3 <- lapply(test2, bestResPair)
# 
# test3 <- lapply(test3, function(x) prepPairs(x, 10))
# 
# pairPlot(test3[[1]], "red")
#pairPlot(test3[[2]], "blue")


edc5555s <- new(Class="PPsearchCompareXL",
                directory="/Users/mtrnka/Projects-Collaboration/GrossLab/Serena/edc_new/edc_four/",
                dataFile="edc5555s.txt",
                modFile="../../module_definitions.txt",
                pdbDir="../../pdb",
                pdbFile="1kx5.pdb",
                chainMapFile="chain_map.txt",
                redState=factor("xlinkedResPairs"))




# 
# test3 <- reduceTo(test2, "xlinkedResPair", "Score.Diff")
# nrow(test3)
# test3 <- reduceTo(test3, "xlinkedResPair", "Score")
# nrow(test3)
# test3 <- reduceTo(test3, "xlinkedResPair", "Peptide.2")
# nrow(test3)
# test3 <- reduceTo(test3, "xlinkedResPair", "Peptide.1")
# nrow(test3)
# test3 <- reduceTo(test3, "xlinkedResPair", "MSMS.Info")
# nrow(test3)
# distanceDir <- "../pdb"
# chains <- readChainMap(paste(distanceDir,"/chain_map.txt",sep=""))
# nucl <- parsePDB(paste(distanceDir,"/1kx5.pdb",sep=""))
# test@dataTable <- measureDistances(test@dataTable, nucl, chains)

# setwd("/Users/mtrnka/Projects-Collaboration/GrossLab/Serena/edc_new/")
# dataTab <- readProspectorXLOutput("edc_dda_res3.txt")
# 
# test <- dataTab[15,]
# p1 <- new(Class="halfLink",
#           sequence=test$Peptide.1,
#           base_sequence=test$DB.Peptide.1,
#           length=nchar(p1@base_sequence),
#           accession=test$Acc.1,
#           xlinked_aa=as.integer(test$XLink.AA.1),
#           species=test$Species.1)
# 
# p2 <- new(Class="halfLink",
#           sequence=test$Peptide.2,
#           base_sequence=test$DB.Peptide.2,
#           length=nchar(p2@base_sequence),
#           accession=test$Acc.2,
#           xlinked_aa=as.integer(test$XLink.AA.2),
#           species=test$Species.2)
# 
# pxl <- new(Class="xLinkedPepPair",
#            pep1=p1,
#            pep2=p2)
# 

#library(xml2)
#xtest <- read_xml("~/prospector.5.19.22/web/cgi-bin/msprod_test1.xml")
#xlist <- as_list(xtest)
#testline <- test2[680,]
#xlist$parameters$data_filename = url_escape(paste(
#    "/Users/mtrnka/Projects-Collaboration/GrossLab/Serena/peaklists/edc_dda/",
#    testline$Fraction,sep=""))
#xlist$parameters$max_charge <- as.character(testline$z)
#xlist$parameters$msms_precursor_charge <- as.character(testline$z)
#xlist$parameters$sequence <- url_escape(testline$Peptide.1)
#xlist$parameters$sequence2 <- url_escape(testline$Peptide.2)
#xlist$parameters$spot_number <- as.character(testline$RT)
#xlist$parameters$spectrum_number <- as.character(testline$Spectrum)



# EvH <- new(Class="PPsearchCompareXL",
#            directory="/Users/mtrnka/Projects-Mine/LumosTest1/",
#            dataFile="EvH1.txt",
#            modFile="module_definitions.txt")
# evh2 <- EvH@dataTable
# table(evh2$Type)
# evh2$Res.1 <- paste(evh2$XLink.AA.1,evh2$Acc.1,sep=":")
# evh2$Res.2 <- paste(evh2$XLink.AA.2,evh2$Acc.2,sep=":")
# evh2$xlink <- ifelse(evh2$Res.1 <= evh2$Res.2,paste(evh2$Res.1,evh2$Res.2,sep="."),paste(evh2$Res.2,evh2$Res.1,sep="."))
# evh2$xlink <- as.factor(evh2$xlink)
# xtable <- table(evh2$xlink, evh2$Type)
# xtable2 <- as.data.frame(xtable)
# library(reshape)
# xtable2 <- cast(xtable2, Var1 ~ Var2, value="Freq", fun.aggregate = sum)
# sum(xtable2$ETD ==0 & xtable2$HCD>0)
# sum(xtable2$ETD >0 & xtable2$HCD==0)
# sum(xtable2$ETD >0 & xtable2$HCD>0)

















#######################################################################
#                   Outline for other possible classes                #
#######################################################################

# setClass(Class="xLinkedResPair",
#          slots=c(
#              pep1="halfLink",
#              pep2="halfLink",
#              distance="numeric",
#              num_respairs="integer")
# )
# 
# setClass(Class="xLinkedPepPair",
#          slots=c(
#              pep1="halfLink",
#              pep2="halfLink",
#              xl_reagent="factor",
#              decoy="logical",
#              num_spec="integer",
#              distance="numeric",
#              dval="numeric",
#              pvalue="numeric")
# )
# 
# setClass(Class="halfLink",
#          slots=c(
#              sequence="character",
#              base_sequence="character",
#              length="integer",
#              accession="character",
#              xlinked_aa="integer",
#              species="character",
#              protein="character",
#              modul="character",
#              chain="character",
#              decoy="logical",
#              intra="integer",
#              inter="integer"),
#          prototype(
#              sequence=NA_character_,
#              base_sequence=NA_character_,
#              length=NA_integer_,
#              accession=NA_character_,
#              xlinked_aa=NA_integer_,
#              species=NA_character_,
#              protein=NA_character_,
#              modul=NA_character_,
#              chain=NA_character_,
#              decoy=NA,
#              intra=NA_integer_,
#              inter=NA_integer_)
# )
# 
# setClass(Class="spectralMatchContainer",
#          slots=c(
#              mz="numeric",
#              z="integer",
#              ppm="numeric",
#              fraction="character",
#              rt="numeric",
#              rt_units="factor", #sec, min
#              spectrum="integer",
#              scan_no="integer",
#              score="numeric",
#              matched="numeric",
#              expect="numeric",
#              intensity="numeric",
#              ion_table="matrix")
# )
# 
# setClass(Class="xLinkedSpectralMatch",
#          slots=c(
#              xlink="xLinkedPepPair",
#              scorediff="numeric",
#              pep1.half_score="numeric",
#              pep1.rank="integer",
#              pep2.half_score="numeric",
#              pep2.rank="integer"),
#          contains="spectralMatchContainer")
