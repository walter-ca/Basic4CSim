
.adaptBackground <- function(simTable, vpStart, vpRegionDist = 1000000, vpCovRate = 0.95, blindFactor = 0.2, useNormDist = FALSE, minVal = 1000) {

    rnVP = row.names(subset(simTable, simTable$start > (vpStart - vpRegionDist) & simTable$end < (vpStart + vpRegionDist)))

    ## additional step for viewpoint area: adapt data for chosen coverage rate around vp
    for (i in rnVP[1]:rnVP[length(rnVP)]) { 
        randomValue = runif(1,0,1)
        if (randomValue > vpCovRate) {
            simTable[i,4] = 0
        } 
        randomValue = runif(1,0,1)
        if (randomValue > vpCovRate) {
            simTable[i,5] = 0
        } 
        if (useNormDist) {
            ## manipulate read counts: less reads for fragments with higher distance from vp
            temp = simTable[i,]$readsL*dnorm(simTable[i,]$start, vpStart, sd = vpRegionDist*1/3)/dnorm(vpStart, vpStart, sd = vpRegionDist*1/3)
            if (simTable[i,4] > minVal) {
                simTable[i,4] = max(minVal, temp) - runif(1, 0, minVal/2)
            }
            temp = simTable[i,]$readsR*dnorm(simTable[i,]$start, vpStart, sd = vpRegionDist*1/3)/dnorm(vpStart, vpStart, sd = vpRegionDist*1/3)
            if (simTable[i,5] > minVal) {
                simTable[i,5] = max(minVal, temp) - runif(1, 0, minVal/2)
            }
        }
    }

    ## penalty for blind fragment read count
    simTable[simTable$isNonBlind == FALSE,]$readsL = simTable[simTable$isNonBlind == FALSE,]$readsL*blindFactor
    simTable[simTable$isNonBlind == FALSE,]$readsR = simTable[simTable$isNonBlind == FALSE,]$readsR*blindFactor
    
    ## round read number to integer
    simTable$readsL = round(simTable$readsL)
    simTable$readsR = round(simTable$readsR)

    return(simTable)
}


setMethod("adaptBackground",
          signature=signature(simTable="data.frame", vpStart="numeric"),
          .adaptBackground)

