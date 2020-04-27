\name{fusions.process}
\alias{fusions.process}
\title{
  LD-RTPCR fusion identification by Sanger
}
\description{
  Automatically interpret gene fusions found by Sanger sequencing using the Ligation-Dependent PCR protocol. 
}
\usage{
  fusions.process(input, design, sheet = NA, output = ".", cores = NA, ...)
}

\arguments{
  \item{input}{
    Single character value, the path to the directory containing the AB1 files to process.
}
  \item{design}{
    Data.frame describing all possible fusions (see Details).
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
    Further arguments are passed to \code{fusions.process.core}.
}
}

\details{
  \code{design} must contain one row for each possible combination of a left primer with a right primer, whether this fusion is expected and relevant or not.
  
  Expected columns in \code{design} are (excluding extra columns required with \code{extra}) :
  \describe{
    \item{left.name}{Character, the name of the left primer.}
    \item{left.seq}{Character (uppercase), the sequence of the left primer (gene-specific part only).}
    \item{left.unileft}{Character (uppercase), the sequence of the left universal primer used for amplification.}
    \item{left.symbol}{Character, the symbol of the gene targeted by the left primer.}
    \item{left.GRCh38}{Character, the genomic coordinates of the last base of the left primer (chromosome:position:strand).}
    \item{left.GRCh38_band}{Character, the cytogenetic location of the gene targeted by the left primer.}
    \item{right.name}{Character, the name of the right primer.}
    \item{right.seq}{Character (uppercase), the sequence of the right primer (gene-specific part only).}
    \item{right.uniright}{Character (uppercase), the sequence of the right universal primer used for amplification.}
    \item{right.symbol}{Character, the symbol of the gene targeted by the right primer.}
    \item{right.GRCh38}{Character, the genomic coordinates of the last base of the right primer (chromosome:position:strand).}
    \item{right.GRCh38_band}{Character, the cytogenetic location of the gene targeted by the right primer.}
    \item{seq_forward}{Character (uppercase), the complete sequence expected in forward sequencing (concatenation of \code{left.unileft}, \code{left.seq}, \code{right.seq}, \code{right.uniright} and the right tail, if any).}
    \item{seq_reverse}{Character (uppercase), the complete sequence expected in reverse sequencing (reverse complement of a concatenation of the left tail, if any, \code{left.unileft}, \code{left.seq}, \code{right.seq}, \code{right.uniright}).}
  }
  
  Please contact the authors to obtain a relevant design object.
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
  \code{\link{GEP.process}}
}

