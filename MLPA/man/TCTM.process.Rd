\name{TCTM.process}
\alias{TCTM.process}
\title{
  LD-RTPCR fusion identification by Sanger
}
\description{
  Automatically interpret gene fusions found by Sanger sequencing using the Ligation-Dependent PCR protocol. 
}
\usage{
  TCTM.process(input, design, sheet = NA, output = ".", cores = NA, ...)
}

\arguments{
  \item{input}{
    Single character value, the path to the directory containing the AB1 files to process.
}
  \item{design}{
    Data.frame describing all possible fusions.
}
  \item{sheet}{
    Single character value, the name and path of a CSV file describing the files to process. 3 columns are expected: \code{ID} which gives a simpler sample name to use in outputs, \code{way} which defines if sequencing was 'forward' or 'reverse', and \code{file} which gives the file name and path relative to the \code{input} argument.
}
  \item{output}{
    Single character value, the path to a directory in which to produce output files (will be created if doesn't yet exists).
}
  \item{cores}{
    Single integer value, the amount of CPUs to use on the local machine to parallelize the computation. If \code{NA}, a guess will be made. If 1, computation will not use the \code{parallel} package at all but only loop over samples.
}
  \item{\dots}{
    Further arguments are passed to \code{TCTM.process.one}.
}
}

\value{
  Invisibly returns the aggregated table of top results for all samples.
  
  Various files are produced, in location set by the \code{output} argument :
  \describe{
    \item{Top.csv}{The aggregated table of top results for all samples.}
    \item{*.pdf}{One plot for each sample, showing the sequencing profile and the best alignments found.}
  }
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{MLPA.process}}
}

