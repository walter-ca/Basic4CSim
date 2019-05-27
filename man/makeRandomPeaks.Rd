\name{makeRandomPeaks}

\alias{makeRandomPeaks}
\alias{makeRandomPeaks,numeric-method}

\title{Make a set of randomly positioned peaks within a predefined near-cis region}

\description{
This helper function creates a set of random peaks with a predefined range for height and width in near-cis. If a fragment library is provided, peak centers are located on a fragment position. 
}

\usage{makeRandomPeaks(vpStart, vpRegionDist = 100000, peakNumber = 6, 
rmaxRange = c(1500, 3000), rsdRange = c(500, 1000), fragmentLibrary = "none")}

\arguments{
    \item{vpStart}{viewpoint position}
    \item{vpRegionDist}{length of viewpoint region to either side of the viewpoint where peaks can be simulated}
    \item{peakNumber}{number of peaks to simulate}
    \item{rmaxRange}{range of simulated peaks' heights (maximum)}
    \item{rsdRange}{range of simulated peaks' widths (standard deviation)}
    \item{fragmentLibrary}{optional fragment library with valid fragment positions; if one is provided, peak centers are located at a valid fragment position}
}

\value{Data frame with random peak data for use in \code{addPeak} when creating similar, but not identical simulated replicate data}

\examples{
    if(interactive()) {
        vpStart = 69999869
        makeRandomPeaks(peakNumber=6, vpStart)
    }
}

\author{Carolin Walter}

\keyword{makeRandomPeaks}

