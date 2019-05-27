
.simulateSample <-function(simFragLibPath, vpChr, vpStart, firstCutter, numberRandom = 3, writeTo = "simulatedSample", printControls = FALSE) {

    ## import simulated library
    simLib = read.table(simFragLibPath, header = TRUE)
    simLib = subset(simLib, simLib$chromosomeName == vpChr)
    simTable = data.frame(simLib[,1:3], "readsL" = 0, "readsR" = 0, simLib[,c(5,9,10)])
    colnames(simTable) = c("chromosomeName", "start", "end", "readsL", "readsR", "isNonBlind", "leftFragEndValid", "rightFragEndValid")

    indexVP = as.numeric(row.names(subset(simTable, simTable$start == vpStart)))
    randomData = makeRandomPeaks(peakNumber = numberRandom, vpStart)

    ## start the actual simulation with adding a power-law based background to fragment ends...
    simTable = addBackground(simTable, vpStart)

    ## .. add peaks / main interactions...
    simTable = addPeaks(simTable, randomData, vpStart)
    ## ... weaker noise peaks ...
    simTable = addRandomPeaks(simTable, vpStart)
    ## ... and adapt the fragments ...
    simTable = adaptBackground(simTable, vpStart, useNormDist = TRUE, vpRegionDist = 150000)
    ## ... and add viewpoint fragment
    simTable[indexVP,]$readsL = round(runif(1, 100000, 200000))

    if (printControls) {
        ## output fragments as control
        write.table(data.frame(simTable[,1:3], round((simTable[,4]+simTable[,5])/2)), file = paste(writeTo, "_fragments.csv", sep = ""), sep = "\t", row.names = FALSE, col.names = FALSE)
        ## plot viewpoint area as control
        nearCis = subset(simTable, simTable$start > (vpStart - 100000) & simTable$end < (vpStart + 100000))
        pdf(file = paste(writeTo, "_near-cis.pdf", sep = ""), width = 9, height = 5)
        plot(c(nearCis$start, nearCis$end), pmin(5000, c(nearCis$readsL, nearCis$readsR)), type = "h", xlab = "fragments", ylab = "read count")
        dev.off()
    }
    
    ## print fastq
    printFastq(simTable, simLib, writeTo, firstCutter)
}


setMethod("simulateSample",
          signature=signature(simFragLibPath="character", vpChr="character", vpStart="numeric", firstCutter="character"),
          .simulateSample)
