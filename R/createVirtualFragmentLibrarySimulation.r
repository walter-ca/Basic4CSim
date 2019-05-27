
.createVirtualFragmentLibrarySimulation_BSgenome <- function(chosenGenome, firstCutter, secondCutter, readLength, onlyNonBlind = TRUE, useOnlyIndex = TRUE, minSize = 0, maxSize = -1, minFragEndSize = 0, maxFragEndSize = 10000000, useAllData = TRUE, chromosomeName = "chr1", libraryName = "default") {

    chromosomeNames = seqnames(chosenGenome)
    
    totalFragments = NULL
    totalFragmentsRev = NULL
    
    ## for each chromosome: split current chromosome at the first cutter sequence
    ## and check for presence of the second cutter within the resulting fragments
    ## --> remove non-unique and blind fragments (if chosen; default == TRUE) for final output
    for (i in 1:length(chromosomeNames)) {
        
        chromosomeToSplit = chosenGenome[[i]]
        
        if (class(chromosomeToSplit) == "MaskedDNAString") {
            chromosomeToSplit = unmasked(chromosomeToSplit)
        }
        
        if (useAllData) {
            message(paste("analyzing ", chromosomeNames[i], "...", sep = ""))
            chromosomeToSplitRev = reverseComplement(chromosomeToSplit)
            currentChromosome = splitChromosome(firstCutter, secondCutter, chromosomeToSplit, chromosomeNames[i])
            currentChromosomeRev = splitChromosome(firstCutter, secondCutter, chromosomeToSplitRev, chromosomeNames[i])
            totalFragments = rbind(totalFragments, currentChromosome)
            totalFragmentsRev = rbind(currentChromosomeRev, totalFragmentsRev)
        } else {
            
            if (length(chromosomeToSplit) > 20000000) {
                message(paste("analyzing ", chromosomeNames[i], "...", sep = ""))
                chromosomeToSplitRev = reverseComplement(chromosomeToSplit)
                currentChromosome = splitChromosome(firstCutter, secondCutter, chromosomeToSplit, chromosomeNames[i])
                currentChromosomeRev = splitChromosome(firstCutter, secondCutter, chromosomeToSplitRev, chromosomeNames[i])
                totalFragments = rbind(totalFragments, currentChromosome)
                totalFragmentsRev = rbind(currentChromosomeRev, totalFragmentsRev)
            } else {
                message(paste("skipping ", chromosomeNames[i], " due to its length", sep = ""))
            }
        }
    }
    
    if (!(is.null(totalFragments))) {
        createVirtualFragmentLibrarySimulationMain(totalFragments, totalFragmentsRev, firstCutter, secondCutter, readLength, onlyNonBlind, useOnlyIndex, minSize, maxSize, minFragEndSize, maxFragEndSize, chromosomeName, libraryName)
    } else {
        message("No fragments created; please use 'useAllData = TRUE' for genomes with small chromosomes")
    }
}


.createVirtualFragmentLibrarySimulation_DNAString <- function(chosenGenome, firstCutter, secondCutter, readLength, onlyNonBlind = TRUE, useOnlyIndex = TRUE, minSize = 0, maxSize = -1, minFragEndSize = 0, maxFragEndSize = 10000000, useAllData = TRUE, chromosomeName = "chr1", libraryName = "default") {
    
    totalFragments = NULL
    totalFragmentsRev = NULL

    ## only one chromosome present: split current chromosome at the first cutter sequence
    ## and check for presence of the second cutter within the resulting fragments
    ## --> remove non-unique and blind fragments (if chosen; default == TRUE) for final output
    if (useAllData) {
        chromosomeToSplit = chosenGenome
        chromosomeToSplitRev = reverseComplement(chosenGenome)
        totalFragments = splitChromosome(firstCutter, secondCutter, chromosomeToSplit, chromosomeName)
        totalFragmentsRev = splitChromosome(firstCutter, secondCutter, chromosomeToSplitRev, chromosomeName)
        createVirtualFragmentLibrarySimulationMain(totalFragments, totalFragmentsRev, firstCutter, secondCutter, readLength, onlyNonBlind, useOnlyIndex, minSize, maxSize, minFragEndSize, maxFragEndSize, chromosomeName, libraryName)
    } else {
        chromosomeToSplit = chosenGenome
        
        if (length(chromosomeToSplit) > 20000000) {      
            chromosomeToSplitRev = reverseComplement(chosenGenome)
            totalFragments = splitChromosome(firstCutter, secondCutter, chromosomeToSplit, chromosomeName)
            totalFragmentsRev = splitChromosome(firstCutter, secondCutter, chromosomeToSplitRev, chromosomeName)
            createVirtualFragmentLibrarySimulationMain(totalFragments, totalFragmentsRev, firstCutter, secondCutter, readLength, onlyNonBlind, useOnlyIndex, minSize, maxSize, minFragEndSize, maxFragEndSize, chromosomeName, libraryName)
        } else {
            message("Chromosome's length is below threshold; use 'useAllData = TRUE' (default) to create a virtual fragment library for any chromosome length")
        }
    }
}



.splitChromosome_DNAString <- function(firstCutter, secondCutter, chromosomeToSplit, chromosomeName) {
    
    ## first step: get position of the cutter sequences and calculate start and end of fragment in between
    ## --> cutter sequence neglected for fragment to provide disjunct fragment intervals
    rawFragments = matchPattern(firstCutter, chromosomeToSplit)
    fragmentStart = c(1, end(rawFragments) + 1)
    fragmentEnd = c(start(rawFragments) - 1, length(chromosomeToSplit))

    ## second step: read chromosome as string and check which fragments are blind (no second cutter present
    ## or non-blind (second cutter sequence present within fragment sequence)
    chromosomeTotal = toString(chromosomeToSplit)
    fragmentSequences = unlist(strsplit(chromosomeTotal, split=toupper(firstCutter)))
    secondCutterPresent = grepl(toupper(secondCutter), fragmentSequences)

    ## sequences at start and end of chromosome not counted as non-blind fragments 
    ## --> valid non-blind fragments: only [FCE ... SCE ... FCE], not [1 ... (SCE) ... FCE] or [FCE ... (SCE) ... END] 
    secondCutterPresent[1] = FALSE
    secondCutterPresent[length(secondCutterPresent)] = FALSE
    
    emptyLastFrag = FALSE
    
    if (length(fragmentEnd) > length(secondCutterPresent)) {
        secondCutterPresent = c(secondCutterPresent, FALSE)
        fragmentSequences = c(fragmentSequences, "")
    emptyLastFrag = TRUE
    }
    
    ## fragments total
    fragmentTable = data.frame(chromosomeName, fragmentStart, fragmentEnd, secondCutterPresent, fragmentSequences)
    fragmentTable[,5] = as.vector(fragmentTable[,5])
  
    ## delete possible empty fragment if cutter sequence is at the end of the genome
    if (emptyLastFrag) {
        fragmentTable = fragmentTable[-nrow(fragmentTable),]
    }
    
    secondCutterPos = gregexpr(toupper(secondCutter), fragmentTable[,5])
    secondCutterFirst = NULL
    secondCutterLast = NULL
    
    for (i in 1:length(secondCutterPos)) {
        secondCutterFirst[i] = secondCutterPos[[i]][1]
        secondCutterLast[i] = secondCutterPos[[i]][length(secondCutterPos[[i]])]
    }
    
    fragmentTable["secondCutterFirst"] = secondCutterFirst
    fragmentTable["secondCutterLast"] = secondCutterLast
  
    ## return bed-like format: chromosome name, fragment start, fragment end plus fragment sequence and second cutter positions
    return(fragmentTable)
}




.createVirtualFragmentLibrarySimulationMain <- function(totalFragments, totalFragmentsRev, firstCutter, secondCutter, readLength, onlyNonBlind = TRUE, useOnlyIndex = TRUE, minSize = 0, maxSize = -1, minFragEndSize = 0, maxFragEndSize = 10000000, chromosomeName = "chr1", libraryName = "default") {
    
    ## if necessary, change chromosome name from "chr1"..."chrY" to "1"..."Y"
    if (useOnlyIndex) {
        totalFragments$chromosomeName = sub("chr", "", totalFragments$chromosomeName)
        totalFragmentsRev$chromosomeName = sub("chr", "", totalFragmentsRev$chromosomeName)
    }
    
    totalFragments["fragmentLength"] = nchar(totalFragments$fragmentSequences)
    totalFragmentsRev["fragmentLength"] = nchar(totalFragmentsRev$fragmentSequences)

    ## calculate frag end lengths
    leftFragEndLength = ifelse(totalFragments$secondCutterFirst == -1, totalFragments$fragmentLength, totalFragments$secondCutterFirst - 1)
    rightFragEndLength = ifelse(totalFragments$secondCutterFirst == -1, totalFragments$fragmentLength, totalFragments$fragmentLength - totalFragments$secondCutterLast + 1 - nchar(secondCutter))
    leftFragEndLengthRev = ifelse(totalFragmentsRev$secondCutterFirst == -1, totalFragmentsRev$fragmentLength, totalFragmentsRev$secondCutterFirst - 1)
    rightFragEndLengthRev = ifelse(totalFragmentsRev$secondCutterFirst == -1, totalFragmentsRev$fragmentLength, totalFragmentsRev$fragmentLength - totalFragmentsRev$secondCutterLast + 1 - nchar(secondCutter))
  
    ## valid 4C-seq reads map to one of the experiment's fragment ends and continue with the downstream sequence
    fragSeqStart = substr(totalFragments$fragmentSequences, start = 1, stop = readLength)
    fragSeqStartRev = substr(totalFragmentsRev$fragmentSequences, start = 1, stop = readLength)

    ## check fragment ends for uniqueness
    uniqueFromStart = !duplicated(c(fragSeqStart, fragSeqStartRev))
    uniqueFromEnd = !duplicated(c(fragSeqStart, fragSeqStartRev), fromLast = TRUE)
    uniqueSeq = uniqueFromStart & uniqueFromEnd

    uniqueFragStart = uniqueSeq[1:nrow(totalFragments)]
    uniqueFragEnd = rev(uniqueSeq[(nrow(totalFragments)+1):(nrow(totalFragments)*2)])
    
    rm(totalFragmentsRev)
  
    ## mark fragments with length < readLength as not usable, check minimum and maximum fragment end sizes
    minFragEndSize = max(readLength, minFragEndSize)
    leftLongEnough = ifelse((leftFragEndLength >= minFragEndSize & leftFragEndLength <= maxFragEndSize), TRUE, FALSE)
    rightLongEnough = ifelse((rightFragEndLength >= minFragEndSize & rightFragEndLength <= maxFragEndSize), TRUE, FALSE)
    
    leftFragEndValid = uniqueFragStart & leftLongEnough
    rightFragEndValid = uniqueFragEnd & rightLongEnough
  
    ## centre of fragment: middle between first and last second cutter site, if second cutter present
    fragmentCentre = ifelse(totalFragments$secondCutterFirst == -1, round((totalFragments$fragmentStart + totalFragments$fragmentEnd) / 2), round((totalFragments$secondCutterFirst + totalFragments$secondCutterLast) / 2) + totalFragments$fragmentStart + 1)
    
    isNonBlind = totalFragments$secondCutterPresent
    
    leftSequence = substr(totalFragments$fragmentSequences, start = 1, stop = readLength)  
    rightSequence = substr(totalFragments$fragmentSequences, start = nchar(totalFragments$fragmentSequences) - readLength + 1, stop = nchar(totalFragments$fragmentSequences))
    
    leftSequence = ifelse (leftSequence == "", "-", leftSequence)
    rightSequence = ifelse (rightSequence == "", "-", rightSequence)
    
    fragData = data.frame(totalFragments[,1:3], fragmentCentre, isNonBlind, "fragmentLength" = totalFragments$fragmentLength, leftFragEndLength, rightFragEndLength, leftFragEndValid, rightFragEndValid, leftSequence, rightSequence)

    if (onlyNonBlind) {
        fragData = subset(fragData, fragData$isNonBlind == TRUE)
    }
    
    ## if chosen, keep only fragments with a fragment length of more than X bp...
    fragData = subset(fragData, fragData$fragmentLength >= minSize)
    ## ... and less than Y bp
    if (maxSize != -1) {
        fragData = subset(fragData, fragData$fragmentLength <= maxSize)
    }
    
    ## output: fragment coordinates and flags for uniqueness
    if (libraryName == "default") {
        if (onlyNonBlind) {
            libraryName = paste("fragments_", firstCutter, "_", secondCutter, ".csv", sep = "")
        } else {
            libraryName = paste("fragments_", firstCutter, "_", secondCutter, "_with_blind_fragments", ".csv", sep = "")
        }
    }
  
    if (libraryName == "") {
        return(fragData)
    } else {
        write.table(fragData, file = libraryName, sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    }
}


setMethod("createVirtualFragmentLibrarySimulation",
          signature=signature(chosenGenome="BSgenome", firstCutter="character", secondCutter="character", readLength="numeric"),
          .createVirtualFragmentLibrarySimulation_BSgenome)

setMethod("createVirtualFragmentLibrarySimulation",
          signature=signature(chosenGenome="DNAString", firstCutter="character", secondCutter="character", readLength="numeric"),
          .createVirtualFragmentLibrarySimulation_DNAString)


setMethod("splitChromosome",
          signature=signature(firstCutter="character", secondCutter="character", chromosomeToSplit="DNAString", chromosomeName="character"),
          .splitChromosome_DNAString)


setMethod("createVirtualFragmentLibrarySimulationMain",
          signature=signature(totalFragments="data.frame", totalFragmentsRev="data.frame", firstCutter="character", secondCutter="character", readLength="numeric"),
          .createVirtualFragmentLibrarySimulationMain)
