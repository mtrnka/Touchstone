setClass(Class="PPsearchCompareXL",
         representation=representation(
             parentDir="character",
             dataTable="data.frame",
             peaklistDir="character",
             fastaDB="character",
             modulFile="list",
             chainMap="character",
             pdbDir="character")
)

setMethod(f="initialize",
          signature="PPsearchCompareXL",
          definition=function(.Object, directory, dataFile, modFile){
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
              #Add validity check on .Object
              #validObject(.Object)
              return(.Object)
          }
)

setMethod(f="show",
          signature="PPsearchCompareXL",
          definition=function(object){
              cat("*** Class: PPsearchCompareXL *** \n")
              cat("*** Truncated Data table *** \n")
              print(head(object@dataTable, 5))
          }
)

setGeneric("getSearchTable", function(object) {
    standardGeneric("getSearchTable")
})

setMethod("getSearchTable", "PPsearchCompareXL",
          function(object){
              return(object@dataTable)
          }
)

setClass(Class="xLinkedResPair",
         representation=representation(
             pep1="halfLink",
             pep2="halfLink",
             distance="numeric",
             num_respairs="integer")
)

setClass(Class="xLinkedPepPair",
         representation=representation(
             pep1="halfLink",
             pep2="halfLink",
             xl_reagent="factor",
             decoy="logical",
             num_spec="integer",
             distance="numeric",
             dval="numeric",
             pvalue="numeric")
)

setClass(Class="halfLink",
         representation=representation(
             sequence="character",
             base_sequence="character",
             length="integer",
             accession="character",
             xlinked_aa="integer",
             species="character",
             protein="character",
             modul="character",
             chain="character",
             decoy="logical",
             intra="integer",
             inter="integer"),
         prototype(
             sequence=NA_character_,
             base_sequence=NA_character_,
             length=NA_integer_,
             accession=NA_character_,
             xlinked_aa=NA_integer_,
             species=NA_character_,
             protein=NA_character_,
             modul=NA_character_,
             chain=NA_character_,
             decoy=NA,
             intra=NA_integer_,
             inter=NA_integer_)
)

setClass(Class="spectralMatchContainer",
         representation=representation(
             mz="numeric",
             z="integer",
             ppm="numeric",
             fraction="character",
             rt="numeric",
             rt_units="factor", #sec, min
             spectrum="integer",
             scan_no="integer",
             score="numeric",
             matched="numeric",
             expect="numeric",
             intensity="numeric",
             ion_table="matrix")
)

setClass(Class="xLinkedSpectralMatch",
         representation=representation(
             xlink="xLinkedPepPair",
             scorediff="numeric",
             pep1.half_score="numeric",
             pep1.rank="integer",
             pep2.half_score="numeric",
             pep2.rank="integer"),
         contains="spectralMatchContainer")

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
                                  c("Range_low","Module")]
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

measureDistance <- function(searchTable, parsedPDB, chainMap) {
    chainLookup <- function(proteinName) {chainMap[[proteinName]]}
    pdbList <- list(parsedPDB)
    chains1 <- sapply(searchTable$Protein.1, chainLookup)
    chains2 <- sapply(searchTable$Protein.2, chainLookup)
    distance <- mapply(
        multiEuclideanDistance,
        pdbList,
        searchTable$XLink.AA.1,
        chains1,
        searchTable$Xlink.AA.2,
        chains2)
    searchTable$distance <- distance
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

#tests
test <- new(Class="PPsearchCompareXL",
            directory="/Users/mtrnka/Projects-Collaboration/GrossLab/Serena/edc_new/",
            dataFile="edc_dda_res3.txt",
            modFile="../module_definitions.txt")

chains <- readChainMap("/Users/mtrnka/Projects-Collaboration/GrossLab/Serena/pdb/chain_map.txt")
nucl <- parsePDB("1kx5")


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

