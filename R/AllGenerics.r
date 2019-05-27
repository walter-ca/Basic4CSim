
# generics for data simulation
setGeneric("adaptBackground", signature=c("simTable", "vpStart"),
    function(simTable, vpStart, vpRegionDist = 1000000, vpCovRate = 0.95, blindFactor = 0.2, useNormDist = FALSE, minVal = 1000)
        standardGeneric("adaptBackground"))

setGeneric("addBackground", signature=c("simTable", "vpStart"),
    function(simTable, vpStart, powerlawAlpha = 1.35, imRegionDist = 1000000, bgFragCov = 0.01, imFragCov = 0.8, maxFragBG = 1000, maxFragIM = 1500)
        standardGeneric("addBackground"))

setGeneric("addCisPeaks", signature=c("simTable", "peaks"),
    function(simTable, peaks, minRandomAdd = -500, maxRandomAdd = 500, chanceRandom = 0.8, differenceRandom = 0.2, zeroRate = 0, blindFactor = 0.2, block = FALSE)
        standardGeneric("addCisPeaks"))

setGeneric("addCisPeaks", signature=c("simTable", "peaks"),
    function(simTable, peaks, minRandomAdd = -500, maxRandomAdd = 500, chanceRandom = 0.8, differenceRandom = 0.2, zeroRate = 0, blindFactor = 0.2, block = FALSE)
        standardGeneric("addCisPeaks"))

setGeneric("addPeaks", signature=c("simTable", "randomPeaks", "vpStart"),
    function(simTable, randomPeaks, vpStart, maxVPArea = 7500, sdVP = 2000, vpRegionDist = 100000, minRandomVPPeak = -500, maxRandomVPPeak = 500, chanceRandom = 0.8, differenceRandom = 0.2)
        standardGeneric("addPeaks"))

setGeneric("addRandomPeaks", signature=c("simTable", "vpStart"),
    function(simTable, vpStart, maxVPArea = 1500, sdVP = 10000, vpRegionDist = 200000, otherPeaks = 5, rmaxRange = c(500, 1500), rsdRange = c(500, 1000), minRandomVPPeak = -200, maxRandomVPPeak = 500, chanceRandomVPPeak = 0.9)
        standardGeneric("addRandomPeaks"))

setGeneric("makeRandomPeaks", signature=c("vpStart"),
    function(vpStart, vpRegionDist = 100000, peakNumber = 6, rmaxRange = c(1500, 3000), rsdRange = c(500, 1000), fragmentLibrary = "none")
        standardGeneric("makeRandomPeaks"))

setGeneric("randomPowerLaw", signature=c("n"),
    function(n, alpha=2, xmin=1)
        standardGeneric("randomPowerLaw"))

setGeneric("simulateSample", signature=c("simFragLibPath", "vpChr", "vpStart", "firstCutter"),
    function(simFragLibPath, vpChr, vpStart, firstCutter, numberRandom = 3, writeTo = "simulatedSample", printControls = FALSE)
        standardGeneric("simulateSample"))


# generics for fragment library creation
setGeneric("createVirtualFragmentLibrarySimulation", signature=c("chosenGenome", "firstCutter", "secondCutter", "readLength"),
    function(chosenGenome, firstCutter, secondCutter, readLength, onlyNonBlind = TRUE, useOnlyIndex = FALSE, minSize = 0, maxSize = -1, minFragEndSize = 0, maxFragEndSize = 10000000, useAllData = TRUE, chromosomeName = "chr1", libraryName = "default")
        standardGeneric("createVirtualFragmentLibrarySimulation"))

setGeneric("splitChromosome", signature=c("firstCutter", "secondCutter", "chromosomeToSplit", "chromosomeName"),
    function(firstCutter, secondCutter, chromosomeToSplit, chromosomeName, onlyNonBlind = TRUE)
        standardGeneric("splitChromosome"))

setGeneric("createVirtualFragmentLibrarySimulationMain", signature=c("totalFragments", "totalFragmentsRev", "firstCutter", "secondCutter", "readLength"),
    function(totalFragments, totalFragmentsRev, firstCutter, secondCutter, readLength, onlyNonBlind = TRUE, useOnlyIndex = FALSE, minSize = 0, maxSize = -1, minFragEndSize = 0, maxFragEndSize = 10000000, chromosomeName = "chr1", libraryName = "default")
        standardGeneric("createVirtualFragmentLibrarySimulationMain"))


# generics for data export
setGeneric("printFastq", signature=c("simTable", "simLib", "ID", "firstCutter"),
    function(simTable, simLib, ID, firstCutter, lengthRead = 26)
        standardGeneric("printFastq"))

