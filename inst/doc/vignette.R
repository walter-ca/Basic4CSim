### R code from vignette source 'vignette.Rnw'

###################################################
### code chunk number 1: vignette.Rnw:36-38
###################################################
options(width=90)
options(continue=" ")


###################################################
### code chunk number 2: preliminaries
###################################################
library(Basic4CSim)


###################################################
### code chunk number 3: createVirtualFragmentLibrarySimulation
###################################################
testGenome = DNAString("TCATGAAAGATCGGGGCATGTT")
fragmentData = createVirtualFragmentLibrarySimulation(
    chosenGenome = testGenome, firstCutter = "catg", 
    secondCutter = "gatc", readLength = 3,  chromosomeName = 
    "test", libraryName = "")
head(fragmentData)


###################################################
### code chunk number 4: makeRandomPeaks
###################################################
vpStart = 69999869
makeRandomPeaks(peakNumber=6, vpStart)


###################################################
### code chunk number 5: addBackground
###################################################
simLibFile <- system.file("extdata", "simLib_short.csv", 
    package="Basic4CSim")
# load virtual simulation fragment library
simLib = read.csv(simLibFile, sep = "\t", header = TRUE)
vpStart = 69999869
# create empty simulation fragment count table
simTable = data.frame(simLib[,1:3], "readsL" = 0, "readsR" = 0, 
    simLib[,c(5,9,10)])
colnames(simTable) = c("chromosomeName", "start", "end", 
    "readsL", "readsR", "isNonBlind", "leftFragEndValid", "rightFragEndValid")
# add actual background
set.seed(42)
simTable = addBackground(simTable, vpStart)
head(simTable)


###################################################
### code chunk number 6: addPeaks
###################################################
simTableFile <- system.file("extdata", "simTable_bg.csv", 
    package="Basic4CSim")
simTable = read.csv(simTableFile, sep = "\t", header = TRUE)
vpStart = 69999869
set.seed(42)
randomPeaks = makeRandomPeaks(peakNumber=6, vpStart)
simTable = addPeaks(simTable, randomPeaks, vpStart)
head(simTable)


###################################################
### code chunk number 7: addRandomPeaks
###################################################
simTableFile <- system.file("extdata", "simTable_peaks.csv", 
    package="Basic4CSim")
simTable = read.csv(simTableFile, sep = "\t", header = TRUE)
vpStart = 69999869
set.seed(42)
simTable = addRandomPeaks(simTable, vpStart)
head(simTable)


###################################################
### code chunk number 8: adaptBackground
###################################################
simTableFile <- system.file("extdata", "simTable_noise.csv", 
    package="Basic4CSim")
simTable = read.csv(simTableFile, sep = "\t", header = TRUE)
vpStart = 69999869
set.seed(42)
simTable = adaptBackground(simTable, vpStart, useNormDist = TRUE,
    vpRegionDist = 150000)
head(simTable)


###################################################
### code chunk number 9: printFastq (eval = FALSE)
###################################################
## simLibFile <- system.file("extdata", "simLib_short.csv", 
##     package="Basic4CSim")
## simLib = read.csv(simLibFile, sep = "\t", header = TRUE)
## simTableFile <- system.file("extdata", "simTable_adapted.csv", 
##     package="Basic4CSim")
## simTable = read.csv(simTableFile, sep = "\t", header = TRUE)
## printFastq(simTable, simLib, "test", "CATG")


###################################################
### code chunk number 10: sessionInfo
###################################################
sessionInfo()


