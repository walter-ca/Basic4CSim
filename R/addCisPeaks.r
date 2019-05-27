
.addCisPeaks <- function(simTable, peaks, minRandomAdd = -500, maxRandomAdd = 500, chanceRandom = 0.8, differenceRandom = 0.2, zeroRate = 0, blindFactor = 0.2, block = FALSE) {

    ## make random peaks within specified range
    peaks$max = peaks$max * round(runif(nrow(peaks), 1-differenceRandom, 1+differenceRandom))
    peaks$sd = peaks$sd * round(runif(nrow(peaks), 1-differenceRandom, 1+differenceRandom))

    simTable$index = 1:nrow(simTable)

    for (j in 1:nrow(peaks)) {

        region = subset(simTable, simTable$start >= (peaks$mean[j] - 3*peaks$sd[j]) & simTable$end <= (peaks$mean[j] + 3*peaks$sd[j]))
        
        if (nrow(region) > 0) {
            
            for (k in 1:nrow(region)) {

                ## adapt interaction form if necessary (block or normal distribution curve)
                if (block == FALSE) {
                    simTable[region$index[k],4] = max(simTable[region$index[k],4],
                                round(peaks$max[j] / (1/(peaks$sd[j]*sqrt(2*pi))) * dnorm((simTable[region$index[k],2]),
                                                                                          mean=peaks$mean[j], sd = peaks$sd[j])))
                    simTable[region$index[k],5] = max(simTable[region$index[k],5],
                                round(peaks$max[j] / (1/(peaks$sd[j]*sqrt(2*pi))) * dnorm((simTable[region$index[k],3]),
                                                                                          mean=peaks$mean[j], sd = peaks$sd[j])))
                }

                ## add noise within specified range to fragments
                rn = runif(2,0,1)
                if (rn[1] < chanceRandom) {
                    simTable[region$index[k],4] = max(0, simTable[region$index[k],4] + runif(1, minRandomAdd, maxRandomAdd))
                }
                if (rn[2] < chanceRandom) {
                    simTable[region$index[k],5] = max(0, simTable[region$index[k],5] + runif(1, minRandomAdd, maxRandomAdd))
                }

                ## reduce simulated read count for blind fragments
                if (simTable[region$index[k],6] == FALSE) {
                    simTable[region$index[k],4] = simTable[region$index[k],4] * blindFactor
                    simTable[region$index[k],5] = simTable[region$index[k],5] * blindFactor
                }

                ## adapt fragment coverage rate (-> add additional '0 read fragments' if zero rate is > 0)
                rn = runif(2,0,1)
                if (rn[1] < zeroRate) {
                    simTable[region$index[k],4] = 0
                }
                if (rn[2] < zeroRate) {
                    simTable[region$index[k],5] = 0
                }            
                
            }
        }
    }
  
    return(simTable[,1:8])
}


setMethod("addCisPeaks",
          signature=signature(simTable="data.frame", peaks="data.frame"),
          .addCisPeaks)
