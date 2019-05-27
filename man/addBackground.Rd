\name{addBackground}

\alias{addBackground}
\alias{addBackground,data.frame,numeric-method}

\title{Add background for a simulated 4C-seq sample}

\description{
This helper function adds background noise to an empty 4C-seq sample. If the virtual fragment library of the 4C-seq sample does already include counts, those may be overwritten. 
}

\usage{addBackground(simTable, vpStart, powerlawAlpha = 1.35, imRegionDist = 1000000, 
bgFragCov = 0.01, imFragCov = 0.8, maxFragBG = 1000, maxFragIM = 1500)}

\arguments{
    \item{simTable}{virtual simulation fragment library (should be empty)}
    \item{vpStart}{position of the experiment's viewpoint}
    \item{powerlawAlpha}{alpha for powerlaw background distribution (default: 1.35)}
    \item{imRegionDist}{length of intermediate region between viewpoint and far-cis region}
    \item{bgFragCov}{coverage rate per fragment for background region fragments (typically low)}
    \item{imFragCov}{coverage rate per fragment for intermediate region fragments}
    \item{maxFragBG}{maximum signal strength / read number for background fragments}
    \item{maxFragIM}{maximum signal strength / read number for fragments in the intermediate region between viewpoint and pure background}
}

\value{Data frame with background reads}

\examples{
    if(interactive()) {
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
    }
}

\author{Carolin Walter}

\keyword{addBackground}

