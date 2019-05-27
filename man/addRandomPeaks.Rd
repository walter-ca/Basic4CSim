\name{addRandomPeaks}

\alias{addRandomPeaks}
\alias{addRandomPeaks,data.frame,numeric-method}

\title{Add random near-cis peaks to a simulated 4C-seq sample}

\description{
This helper function adds random regions of higher signal (peaks) in near-cis to a 4C-seq sample. Random noise as well as modifiers for blind fragments can be added. The default shape of peaks resembles a normal distribution, however block shaped peaks can be created as well. 
}

\usage{addRandomPeaks(simTable, vpStart, maxVPArea = 1500, sdVP = 10000, 
vpRegionDist = 200000, otherPeaks = 5, rmaxRange = c(500, 1500), 
rsdRange = c(500, 1000), minRandomVPPeak = -200, maxRandomVPPeak = 500, 
chanceRandomVPPeak = 0.9)}

\arguments{
    \item{simTable}{virtual simulation fragment library}
    \item{vpStart}{viewpoint position}
    \item{maxVPArea}{height of viewpoint peak (with possible exception of actual viewpoint fragment)}
    \item{sdVP}{width of viewpoint peak}
    \item{vpRegionDist}{length of near-cis region where randomly generated peaks are simulated}
    \item{otherPeaks}{number of random peaks to generate}
    \item{rmaxRange}{minimum and maximum height for random peaks}
    \item{rsdRange}{minimum and maximum width for random peaks}
    \item{minRandomVPPeak}{minimum for randomly chosen summand to fragment read count (can be negative)}
    \item{maxRandomVPPeak}{maximum for randomly chosen summand to fragment read count (can be negative)}
    \item{chanceRandomVPPeak}{chance for occurrence of random summands}
}

\value{Data frame with additional random peaks and viewpoint modification}

\examples{
    if(interactive()) {
        simTableFile <- system.file("extdata", "simTable_peaks.csv", 
            package="Basic4CSim")
        simTable = read.csv(simTableFile, sep = "\t", header = TRUE)
        vpStart = 69999869
        set.seed(42)
        simTable = addRandomPeaks(simTable, vpStart)
        head(simTable)
    }
}

\author{Carolin Walter}

\keyword{addRandomPeaks}

