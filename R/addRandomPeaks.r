
.addRandomPeaks <- function(simTable, vpStart, maxVPArea = 1500, sdVP = 10000, vpRegionDist = 200000, otherPeaks = 5, rmaxRange = c(500, 1500), rsdRange = c(500, 1000), minRandomVPPeak = -200, maxRandomVPPeak = 500, chanceRandomVPPeak = 0.9) {

    ## vp peaks
    randomMeans = runif(otherPeaks,vpStart-vpRegionDist,vpStart+vpRegionDist)
    rsd = runif(otherPeaks,rsdRange[1],rsdRange[2])
    rmax = runif(otherPeaks, rmaxRange[1], rmaxRange[2])

    rnVP = row.names(subset(simTable, simTable$start > (vpStart - vpRegionDist) & simTable$end < (vpStart + vpRegionDist)))
    randomVPPeaks = data.frame("ID" = 1:otherPeaks, "mean" = randomMeans, "max" = round(rmax, 2), "sd" = round(rsd, 2))
    
    for (i in rnVP[1]:rnVP[length(rnVP)]) {

        ## VP
        if (maxVPArea > 0) {
            tempSimValue = max(simTable[i,4], round(maxVPArea / (1/(sdVP*sqrt(2*pi))) * dnorm((simTable[i,2]), mean=vpStart, sd = sdVP)))
            if (tempSimValue > simTable[i,4]) {
                if (randomValue < chanceRandomVPPeak) {
                    tempSimValue = tempSimValue + randomAdd
                }
                simTable[i,4] = max(0, tempSimValue)
            } 
            tempSimValue = max(simTable[i,5], round(maxVPArea / (1/(sdVP*sqrt(2*pi))) * dnorm((simTable[i,3]), mean=vpStart, sd = sdVP)))
            if (tempSimValue > simTable[i,5]) {
                if (randomValue < chanceRandomVPPeak) {
                    tempSimValue = tempSimValue + randomAdd
                }
                simTable[i,5] = max(0, tempSimValue)
            }  
        }

        ## VP area
        for (j in 1:otherPeaks) {
            
            randomAdd = round(runif(1,minRandomVPPeak,maxRandomVPPeak))
            randomValue = runif(1,0,1)
            
            tempSimValue = round(rmax[j] / (1/(rsd[j]*sqrt(2*pi))) * dnorm((simTable[i,2]), mean = randomMeans[j], sd = rsd[j]))
            if (tempSimValue > simTable[i,4]) {
                if (randomValue < chanceRandomVPPeak) {
                    tempSimValue = tempSimValue + randomAdd
                }
                simTable[i,4] = max(0, tempSimValue)
            }    
            
            randomAdd = round(runif(1,minRandomVPPeak,maxRandomVPPeak))
            randomValue = runif(1,0,1)
            
            tempSimValue = round(rmax[j] / (1/(rsd[j]*sqrt(2*pi))) * dnorm((simTable[i,3]), mean = randomMeans[j], sd = rsd[j]))
            if (tempSimValue > simTable[i,5]) {
                if (randomValue < chanceRandomVPPeak) {
                    tempSimValue = tempSimValue + randomAdd
                }
                simTable[i,5] = max(0, tempSimValue)
            }
        }
        
        randomVPPeaks = randomVPPeaks[order(randomVPPeaks$mean),]
    }
    
    return(simTable)
}


setMethod("addRandomPeaks",
          signature=signature(simTable="data.frame", vpStart="numeric"),
          .addRandomPeaks)
