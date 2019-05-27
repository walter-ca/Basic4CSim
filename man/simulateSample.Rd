\name{simulateSample}

\alias{simulateSample}
\alias{simulateSample,character,character,numeric,character-method}

\title{Add random near-cis peaks to a simulated 4C-seq sample}

\description{
This function is a wrapper function that creates a simple simulated 4C-seq sample with a default background, some random noise, and near-cis peak regions. 
}

\usage{simulateSample(simFragLibPath, vpChr, vpStart, firstCutter, 
numberRandom = 3, writeTo = "simulatedSample", printControls = FALSE)}

\arguments{
    \item{simFragLibPath}{path for a corresponding virtual simulated fragment library containing fragment end sequences}
    \item{vpChr}{viewpoint chromosome}
    \item{vpStart}{viewpoint position}
    \item{firstCutter}{first restriction enzyme sequence}
    \item{numberRandom}{number of randomly simulated near-cis peaks}
    \item{writeTo}{name for the simulated 4C-seq sample's fastq file}
    \item{printControls}{if TRUE, print control graphics for background and near-cis read distribution}
}

\value{An exported fastq file for a basic simulated 4C-seq sample}

\examples{
    if(interactive()) {
        simulateSample()
    }
}

\author{Carolin Walter}

\keyword{simulateSample}

