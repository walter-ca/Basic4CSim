
.printFastq <- function(simTable, simLib, ID, firstCutter, lengthRead = 26) {

    simFastq = subset(simTable, simTable$readsL > 0 | simTable$readsR > 0)

    sink(paste(ID, "_raw.fastq", sep = ""))

    for (i in 1:nrow(simFastq)) {

        readNumberL = simFastq$readsL[i]
        readNumberR = simFastq$readsR[i]
        readSeq1 = as.character(subset(simLib, simLib$fragmentStart == simFastq$start[i])$leftSequence)
        readSeq2 = as.character(subset(simLib, simLib$fragmentEnd == simFastq$end[i])$rightSequence)

        tempQuality = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
        realQuality = substr(tempQuality, 1, lengthRead + nchar(firstCutter))

        if (nchar(readSeq1) == lengthRead) {

            ## print left fragment end reads
            readSeq1 = paste(firstCutter, readSeq1, sep = "")
            for (j in 1:readNumberL) {
                cat(paste("@SimRead_Left_", i, "_", j, sep = ""))
                cat("\n")
                cat(readSeq1)
                cat("\n")
                cat("+")
                cat("\n")
                cat(realQuality)
                cat("\n")
            }
        }
        if (nchar(readSeq2) == lengthRead) {

            ## print right fragment end reads
            readSeq2 = paste(readSeq2, firstCutter, sep = "")
            readSeq2 = as.character(reverseComplement(DNAString(readSeq2)))
            for (j in 1:readNumberR) {
                cat(paste("@SimRead_Right_", i, "_", j, sep = ""))
                cat("\n")
                cat(readSeq2)
                cat("\n")
                cat("+")
                cat("\n")
                cat(realQuality)
                cat("\n")
            }
        }
    }
    
    sink()
}


setMethod("printFastq",
          signature=signature(simTable="data.frame", simLib="data.frame", ID="character", firstCutter="character"),
          .printFastq)
