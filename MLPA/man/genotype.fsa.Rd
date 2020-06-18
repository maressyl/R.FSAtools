\name{genotype.fsa}
\alias{genotype.fsa}
\title{
  Calls alleles for a SNP genotyping experiment
}

\description{
  Calls alleles using fixed thresholds on ratios in experiments where fragments of different sizes are expected according to the allele.
  
  Assumes all loci are diploid.
}

\usage{
  genotype.fsa(x, homo = 0.85, hetero = c(0.3, 0.7))
}

\arguments{
  \item{x}{
    The \code{fsa} object to use, which must have a \code{peaks} attribute (see \code{\link{peaks.fsa}}). Peaks are expected to be named according to the \code{LOCUS_ALLELE} pattern.
}
  \item{homo}{
    Single numeric value, the ratio of the considered allele signal over all signal from this locus to call a homozygous allele. Similarly alleles below \code{1 - homo} will be called "absent".
}
  \item{hetero}{
    Numeric vector of length two, minimum and maximum ratios of the considered allele signal over all signal from this locus to call a heterozygous allele.
}
}

\value{
  Returns \code{x}, with a new or updated \code{genotypes} attribute, a \code{data.frame} with a row for each locus :
  \item{call}{The concatenation of the two called alleles, possibly with \code{?}.}
  \item{alleles}{The comma-separated list of ratios observed for all alleles.}
  
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
  design <- designFile(system.file("extdata/design_SNP.conf", package="MLPA"))
  
  # Example FSA file provided
  fsa <- read.fsa(system.file("extdata/fsa_SNP/A7840.fsa", package="MLPA"))
  fsa <- align.fsa(
    fsa,
    channel = design$align.fsa$channel,
    outThreshold = design$align.fsa$outThreshold,
    useLadder = design$align.fsa$useLadder
  )
  fsa <- peaks.fsa(fsa, peaks=design$GLOBALS$PEAKS)
  
  # Genotype
  fsa <- genotype.fsa(fsa)
  print(attr(fsa, "genotypes"))
  print(attr(fsa, "calls"))
}

