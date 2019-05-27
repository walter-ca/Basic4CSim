
.makeRandomPeaks <- function(vpStart, vpRegionDist = 100000, peakNumber = 6, rmaxRange = c(1500, 3000), rsdRange = c(500, 1000), fragmentLibrary = "none") {

    # create set of random peaks within a specified region around the viewpoint
    if (fragmentLibrary == "none") {
        randomMeans = c(runif(peakNumber, vpStart-vpRegionDist, vpStart-10000), runif(peakNumber, vpStart+10000, vpStart+vpRegionDist))
        randomMeans = sample(randomMeans, peakNumber)
    } else {
        # if a fragment library with valid fragments is provided, use those as position for random peaks
        # --> this guarantees that random peaks are not allocated to completely repetitive regions, but do contain at least one valid fragment
        fragLib = read.csv(fragmentLibrary, sep = "\t", header = TRUE)[,1:3]
        fragLib = subset(fragLib, fragLib[,2] >= vpStart-vpRegionDist & fragLib[,3] <= vpStart+vpRegionDist)
        fragLib = subset(fragLib, fragLib[,3] <= vpStart-10000 | fragLib[,2] >= vpStart+10000)
        randomMeans = sample(fragLib[,2], peakNumber)
    }
    rsd = runif(peakNumber, rsdRange[1], rsdRange[2])
    rmax = runif(peakNumber, rmaxRange[1], rmaxRange[2])

    randomVPPeaks = data.frame("ID" = 1:peakNumber, "mean" = randomMeans, "max" = round(rmax, 2), "sd" = round(rsd, 2))

    return(randomVPPeaks)
}


setMethod("makeRandomPeaks",
          signature=signature(vpStart="numeric"),
          .makeRandomPeaks)
