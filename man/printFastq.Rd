\name{printFastq}

\alias{printFastq}
\alias{printFastq,data.frame,data.frame,character,character-method}

\title{Add random near-cis peaks to a simulated 4C-seq sample}

\description{
This function collects actual sequence data per fragment for a simulated 4C-seq sample (in a fragment table format), and exports them with a high quality score per read as a fastq file. 
}

\usage{printFastq(simTable, simLib, ID, firstCutter, lengthRead = 26)}

\arguments{
    \item{simTable}{virtual simulation fragment library with fragment read counts representing a simulated 4C-seq sample}
    \item{simLib}{original virtual simulation fragment library, with fragment end sequences for simulation purposes}
    \item{ID}{sample name}
    \item{firstCutter}{first restriction enzyme sequence}
    \item{lengthRead}{read length (to check for too short fragments which can cause errors in other programs when exported)}
}

\value{An exported fastq file for a simulated 4C-seq sample}

\examples{
    if(interactive()) {
        simLibFile <- system.file("extdata", "simLib_short.csv", 
            package="Basic4CSim")
        simLib = read.csv(simLibFile, sep = "\t", header = TRUE)
        simTableFile <- system.file("extdata", "simTable_adapted.csv", 
            package="Basic4CSim")
        simTable = read.csv(simTableFile, sep = "\t", header = TRUE)
        printFastq(simTable, simLib, "test", "CATG")
    }
}

\author{Carolin Walter}

\keyword{printFastq}

