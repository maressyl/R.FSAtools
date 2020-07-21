\name{genotype.fsa}
\alias{genotype.closest.fsa}
\alias{genotype.ratio.fsa}
\alias{genotype.N1.fsa}
\title{
  Calls alleles for a SNP genotyping experiment
}

\description{
  Calls alleles in experiments where fragments of different sizes are expected according to the allele.
  
  \code{genotype.closest.fsa} selects for each allele the closest expected normalized peak height among the 0, 1 and 2 copy values provided in the design.
  
  \code{genotype.ratio.fsa} calls an allele as absent, heterozygous or homozygous using fixed thresholds on the proportion of signal from a locus.
  
  \code{genotype.N1.fsa} calls an allele as present if it exceeds a given proportion of the expected signal for 1 copy, as provided in the design. When a single allele is called present for a locus, it is considered homozygous.
}

\usage{
  genotype.closest.fsa(x)
  genotype.ratio.fsa(x, homo = 0.85, hetero = c(0.3, 0.7))
  genotype.N1.fsa(x, threshold = 0.1)
}

\arguments{
  \item{x}{
    The \code{fsa} object to use, which must have a \code{peaks} attribute (see \code{\link{peaks.fsa}}). Peaks are expected to be named according to the \code{LOCUS_ALLELE} pattern. \code{genotype.closest.fsa()} expected optional N0, N1 and N2 columns with expected normalized heights for 0, 1 and 2 copies of the allele.
}
  \item{homo}{
    Single numeric value, the ratio of the considered allele signal over all signal from this locus to call a homozygous allele. Similarly alleles below \code{1 - homo} will be called "absent".
}
  \item{hetero}{
    Numeric vector of length two, minimum and maximum ratios of the considered allele signal over all signal from this locus to call a heterozygous allele.
}
  \item{threshold}{
    Single numeric value, the ratio of the considered allele signal over the expected value for 1n to call an allele.
}
}

\note{
  All three functions assume peaks to be named according to the following convention : "ALLELE - LOCUS".
}

\value{
  Returns \code{x}, with a new or updated \code{genotypes} attribute, a \code{data.frame} with a row for each locus :
  \item{call}{The concatenation of the two called alleles, possibly with \code{?}.}
  \item{alleles}{The comma-separated list of ratios observed for all alleles (\code{genotype.ratio.fsa()} only).}
  
  A \code{calls} vector attribute is also set, corresponding to the \code{call} column only.
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{generic.process}}
}

\examples{
  # Using a design file
  design <- designFile(system.file("extdata/design_SNP.conf", package="FSAtools"))
  
  # Example FSA file provided
  fsa <- read.fsa(system.file("extdata/fsa_SNP/A7840.fsa", package="FSAtools"))
  fsa <- align.fsa(
    fsa,
    channel = design$align.fsa$channel,
    outThreshold = design$align.fsa$outThreshold,
    useLadder = design$align.fsa$useLadder
  )
  fsa <- peaks.fsa(fsa, peaks=design$GLOBALS$PEAKS)
  
  # Genotype
  fsa <- genotype.ratio.fsa(fsa)
  print(attr(fsa, "genotypes"))
  print(attr(fsa, "calls"))
}

