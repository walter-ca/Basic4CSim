\name{adaptBackground}

\alias{adaptBackground}
\alias{adaptBackground,data.frame,numeric-method}

\title{Adapt background for a simulated 4C-seq sample}

\description{
This helper function adapts the already simulated background for a simulated 4C-seq sample. The read count for blind fragments is reduced by a given factor (default: 0.2), coverage around the viewpoint is reduced to a certain rate (optional), and the general background read count per fragment is reduced with higher distance from the experiment's viewpoint (optional). 
}

\usage{adaptBackground(simTable, vpStart, vpRegionDist = 1000000, vpCovRate = 0.95,
 blindFactor = 0.2, useNormDist = FALSE, minVal = 1000)}

\arguments{
    \item{simTable}{virtual simulation fragment library, with background}
    \item{vpStart}{position of the experiment's viewpoint}
    \item{vpRegionDist}{length of adapted region}
    \item{vpCovRate}{coverage rate of background around viewpoint (default: 0.95)}
    \item{blindFactor}{factor for blind fragment read count (default: 0.2)}
    \item{useNormDist}{if TRUE, fragments with higher distance from viewpoint will receive fewer reads (to blend transition from near-cis to cis)}
    \item{minVal}{minimum read count of fragments to be considered for near-cis adaption}
}

\value{Data frame with adapted simulated background reads}

\examples{
    if(interactive()) {
        simTableFile <- system.file("extdata", "simTable_noise.csv", 
            package="Basic4CSim")
        simTable = read.csv(simTableFile, sep = "\t", header = TRUE)
        vpStart = 69999869
        set.seed(42)
        simTable = adaptBackground(simTable, vpStart, useNormDist = TRUE,
            vpRegionDist = 150000)
        head(simTable)
    }
}

\author{Carolin Walter}

\keyword{adaptBackground}

