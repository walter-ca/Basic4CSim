
.addBackground <- function(simTable, vpStart, powerlawAlpha = 1.35, imRegionDist = 1000000, bgFragCov = 0.01, imFragCov = 0.8, maxFragBG = 1000, maxFragIM = 1500) {

    ## power law distribution: xmin = 1, maximum value limited by maxFragBG
    powerLawBG = randomPowerLaw(round(bgFragCov*nrow(simTable))*100, powerlawAlpha, 1) 
    powerLawBG = round(subset(powerLawBG, (powerLawBG < maxFragBG)))
    simDistrBG = c(rep(0, 2*(round((1-bgFragCov)*nrow(simTable)))), sample(powerLawBG, 2*(round(bgFragCov*nrow(simTable)))))
    simTable[,4] = pmax(simTable[,4], sample(simDistrBG, nrow(simTable)))
    simTable[,5] = pmax(simTable[,5], sample(simDistrBG, nrow(simTable)))
    
    ## higher likelihood for fragment to be covered in intermediate region around vp
    simTableIM = subset(simTable, simTable$start > (vpStart - imRegionDist) & simTable$end < (vpStart + imRegionDist))
    powerLawIM = randomPowerLaw(round(imFragCov*nrow(simTableIM))*100, powerlawAlpha, 1)
    powerLawIM = round(subset(powerLawIM, (powerLawIM < maxFragIM)))
    simDistrIMBG = c(rep(0, 2*(round((1-imFragCov)*nrow(simTableIM)))), sample(powerLawIM, 2*(round(imFragCov*nrow(simTableIM))))) # background distribution in intermediate region
    simTableIM[,4] = sample(simDistrIMBG, nrow(simTableIM))
    simTableIM[,5] = sample(simDistrIMBG, nrow(simTableIM))

    ## put data together and add higher likelihood for fragments in IM area to be covered if the VP gets nearer
    rnIM = as.numeric(rownames(simTableIM))
    maxDNorm = dnorm(vpStart, mean = vpStart, sd = (1/3 * imRegionDist))
    ## left frag ends (IM)
    rVal = runif(nrow(simTableIM), 0, 1)
    testVal = simTableIM$start
    normVal = dnorm(testVal, mean = vpStart, sd = (1/3 * imRegionDist)) / maxDNorm
    factor = rep(0, nrow(simTableIM))
    factor[rVal < normVal] = 1
    simTable[rnIM,4] = (pmax(simTable[rnIM,4], simTableIM[,4]))*factor
    ## right frag ends (IM)
    rVal = runif(nrow(simTableIM), 0, 1)
    testVal = simTableIM$end
    normVal = dnorm(testVal, mean = vpStart, sd = (1/3 * imRegionDist)) / maxDNorm
    factor = rep(0, nrow(simTableIM))
    factor[rVal < normVal] = 1
    simTable[rnIM,5] = (pmax(simTable[rnIM,5], simTableIM[,5]))*factor
    
    return(simTable)
}


setMethod("addBackground",
          signature=signature(simTable="data.frame", vpStart="numeric"),
          .addBackground)
