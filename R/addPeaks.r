
.addPeaks <- function(simTable, randomPeaks, vpStart, maxVPArea = 7500, sdVP = 2000, vpRegionDist = 100000, minRandomVPPeak = -500, maxRandomVPPeak = 500, chanceRandom = 0.8, differenceRandom = 0.2) {

    rnVP = row.names(subset(simTable, simTable$start > (vpStart - vpRegionDist) & simTable$end < (vpStart + vpRegionDist)))
    randomVPPeaks = randomPeaks

    ## adapt peaks if necessary
    randomVPPeaks$max = randomVPPeaks$max * round(runif(nrow(randomPeaks), 1-differenceRandom, 1+differenceRandom))
    randomVPPeaks$sd = randomVPPeaks$sd * round(runif(nrow(randomPeaks), 1-differenceRandom, 1+differenceRandom))


    for (i in rnVP[1]:rnVP[length(rnVP)]) {
        ## VP area
        for (j in 1:nrow(randomVPPeaks)) {
            
            randomAdd = round(runif(1,minRandomVPPeak,maxRandomVPPeak))
            randomValue = runif(1,0,1)

            ## adapt fragments if necessary (i.e. if random adds are > 0)
            tempSimValue = round(randomVPPeaks$max[j] / (1/(randomVPPeaks$sd[j]*sqrt(2*pi))) * dnorm((simTable[i,2]), mean = randomVPPeaks$mean[j], sd = randomVPPeaks$sd[j]))
            if (tempSimValue > simTable[i,4]) {
                if (randomValue < chanceRandom) {
                    tempSimValue = tempSimValue + randomAdd
                }
                simTable[i,4] = max(0, tempSimValue)
            }    
            
            randomAdd = round(runif(1,minRandomVPPeak,maxRandomVPPeak))
            randomValue = runif(1,0,1)

            ## second fragmend end
            tempSimValue = round(randomVPPeaks$max[j] / (1/(randomVPPeaks$sd[j]*sqrt(2*pi))) * dnorm((simTable[i,3]), mean = randomVPPeaks$mean[j], sd = randomVPPeaks$sd[j]))
            if (tempSimValue > simTable[i,5]) {
                if (randomValue < chanceRandom) {
                    tempSimValue = tempSimValue + randomAdd
                }
                simTable[i,5] = max(0, tempSimValue)
            }
        }
        
        ## VP
        if (maxVPArea > 0) {
            simTable[i,4] = max(simTable[i,4], round(maxVPArea / (1/(sdVP*sqrt(2*pi))) * dnorm((simTable[i,2]), mean=vpStart, sd = sdVP)))
            simTable[i,5] = max(simTable[i,5], round(maxVPArea / (1/(sdVP*sqrt(2*pi))) * dnorm((simTable[i,3]), mean=vpStart, sd = sdVP)))  
        }
    }

    return(simTable)
}


setMethod("addPeaks",
          signature=signature(simTable="data.frame", randomPeaks="data.frame", vpStart="numeric"),
          .addPeaks)
