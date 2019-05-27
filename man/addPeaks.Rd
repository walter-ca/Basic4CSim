\name{addPeaks}

\alias{addPeaks}
\alias{addPeaks,data.frame,data.frame,numeric-method}

\title{Add predefined peaks to a simulated 4C-seq sample}

\description{
This helper function adds predefined regions of higher signal (peaks) in near-cis to a 4C-seq sample. Random noise as well as modifiers for blind fragments can be added. 
}

\usage{addPeaks(simTable, randomPeaks, vpStart, maxVPArea = 7500, sdVP = 2000, 
vpRegionDist = 100000, minRandomVPPeak = -500, maxRandomVPPeak = 500, 
chanceRandom = 0.8, differenceRandom = 0.2)}

\arguments{
    \item{simTable}{virtual simulation fragment library}
    \item{randomPeaks}{list of peaks to be simulated: ID (e.g. number), mean (genomic coordinate), max (maximum signal for the peak in question), sd (standard deviation, "width")}
    \item{vpStart}{viewpoint position}
    \item{maxVPArea}{height of viewpoint peak (with possible exception of actual viewpoint fragment)}
    \item{sdVP}{width of viewpoint peak}
    \item{vpRegionDist}{length of near-cis region where randomly generated peaks are simulated}
    \item{minRandomVPPeak}{minimum for randomly chosen summand to fragment read count (can be negative)}
    \item{maxRandomVPPeak}{maximum for randomly chosen summand to fragment read count (can be negative)}
    \item{chanceRandom}{chance for occurrence of random summands}
    \item{differenceRandom}{parameter to manipulate peak heights per peak: (1-dR) <= max <= (1+dR)}
}

\value{Data frame with added intervals of higher signal (peaks) in near-cis}

\examples{
    if(interactive()) {
        simTableFile <- system.file("extdata", "simTable_bg.csv", 
            package="Basic4CSim")
        simTable = read.csv(simTableFile, sep = "\t", header = TRUE)
        vpStart = 69999869
        set.seed(42)
        randomPeaks = makeRandomPeaks(peakNumber=6, vpStart)
        simTable = addPeaks(simTable, randomPeaks, vpStart)
        head(simTable)
    }
}

\author{Carolin Walter}

\keyword{addPeaks}

