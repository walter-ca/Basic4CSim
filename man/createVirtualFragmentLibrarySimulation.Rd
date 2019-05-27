\name{createVirtualFragmentLibrarySimulation}

\alias{createVirtualFragmentLibrarySimulation}
\alias{createVirtualFragmentLibrarySimulation,BSgenome,character,character,numeric-method}
\alias{createVirtualFragmentLibrarySimulation,DNAString,character,character,numeric-method}

\title{Create a virtual simulation fragment library with two fragment end sequences per fragment from a provided genome and two restriction enzymes}

\description{\pkg{Basic4CSim} adapts \pkg{Basic4Cseq}'s \code{createVirtualFragmentLibrary} function to provide sequence data for fragment end data simulation.}

\usage{
createVirtualFragmentLibrarySimulation(chosenGenome, firstCutter, secondCutter, 
readLength, onlyNonBlind = TRUE, useOnlyIndex = FALSE, minSize = 0, maxSize = -1, 
minFragEndSize = 0, maxFragEndSize = 10000000, useAllData = TRUE, 
chromosomeName = "chr1", libraryName = "default")
}

\arguments{
    \item{chosenGenome}{The genome that is to be digested in silico with the provided enzymes; can be an instance of BSgenome or DNAString (\pkg{Basic4Cseq})}
    \item{firstCutter}{First of two restriction enzymes (\pkg{Basic4Cseq})}
    \item{secondCutter}{Second of two restriction enzymes (\pkg{Basic4Cseq})}
    \item{readLength}{Read length for the experiment (\pkg{Basic4Cseq})}
    \item{onlyNonBlind}{Variable that is TRUE (default) if only non-blind fragments are considered (i.e. all blind fragments are removed) (\pkg{Basic4Cseq})}
    \item{useOnlyIndex}{Convenience function to adapt the annotation style of the chromosomes ("chr1", \dots "chrY" or "1", \dots, "Y"); parameter has to be set to match the BAM file in question (\pkg{Basic4Cseq})}
    \item{minSize}{Filter option that allows to delete fragments below a certain size (in bp) (\pkg{Basic4Cseq})}
    \item{maxSize}{Filter option that allows to delete fragments above a certain size (in bp) (\pkg{Basic4Cseq})}
    \item{minFragEndSize}{Filter option that allows to delete fragment ends below a certain size (in bp) (\pkg{Basic4Cseq})}
    \item{maxFragEndSize}{Filter option that allows to delete fragment ends above a certain size (in bp) (\pkg{Basic4Cseq})}
    \item{useAllData}{Variable that indicates if all data of a BSgenome package is to be used. If FALSE, chromosome names including a "_" are removed, reducing the set of chromosomes to (1 \dots 19, X, Y, MT) for the mouse genome or (1 \dots 22, X, Y, MT) for the human genome (\pkg{Basic4Cseq})}
    \item{chromosomeName}{Chromosome name for the virtual fragment library if a \code{DNAString} object is used instead of a \code{BSgenome} object (\pkg{Basic4Cseq})}
    \item{libraryName}{Name of the file the created virtual fragment library is written to. Per default the file is called "fragments_firstCutter_secondCutter.csv". The fragment data is returned as a data frame if and only if an empty character string is chosen as \code{libraryName}. (\pkg{Basic4Cseq})}
}

\value{A tab-separated file with the specified virtual fragment library (cp. \pkg{Basic4Cseq}), with added fragment end sequences}

\examples{
    if(interactive()) {
        library(BSgenome.Ecoli.NCBI.20080805)
        fragmentData = createVirtualFragmentLibrarySimulation(chosenGenome = Ecoli$NC_002655, firstCutter = "catg", secondCutter = "gtac", readLength = 30,  onlyNonBlind = TRUE, chromosomeName = "NC_002655", libraryName = "fragments_Ecoli.csv")
    }
}

\author{Carolin Walter}

\keyword{createVirtualFragmentLibrarySimulation}
