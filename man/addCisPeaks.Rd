\name{addCisPeaks}

\alias{addCisPeaks}
\alias{addCisPeaks,data.frame,data.frame-method}

\title{Add peaks in cis for a simulated 4C-seq sample}

\description{
This helper function adds regions of higher signal (peaks) in cis to a 4C-seq sample. Random noise as well as modifiers for blind fragments can be added. The default shape of peaks resembles a normal distribution, however block shaped peaks can be created as well. 
}

\usage{addCisPeaks(simTable, peaks, minRandomAdd = -500, maxRandomAdd = 500, 
chanceRandom = 0.8, differenceRandom = 0.2, zeroRate = 0, blindFactor = 0.2, block = FALSE)}

\arguments{
    \item{simTable}{virtual simulation fragment library}
    \item{peaks}{list of peaks to be simulated: ID (e.g. number), mean (genomic coordinate), max (maximum signal for the peak in question), sd (standard deviation, "width")}
    \item{minRandomAdd}{minimum for randomly chosen summand to fragment read count (can be negative)}
    \item{maxRandomAdd}{maximum for randomly chosen summand to fragment read count (can be negative)}
    \item{chanceRandom}{chance for occurrence of random summands}
    \item{differenceRandom}{parameter to manipulate peak heights per peak: (1-dR) <= max <= (1+dR)}
    \item{zeroRate}{percentage of fragments in the peak region which are set to 0}
    \item{blindFactor}{factor applied to read counts on blind fragments within peak regions}
    \item{block}{parameter that determines the shape of simulated peaks (default: normal distribution / bell shape)}
}

\value{Data frame with added intervals of higher signal (peaks) in cis}

\author{Carolin Walter}

\keyword{addCisPeaks}

