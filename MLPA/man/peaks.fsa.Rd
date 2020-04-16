\name{peaks.fsa}
\alias{peaks.fsa}
\title{
  Get maximal value in ranges
}
\description{
  Look for the maximal value in one or many ranges, typically for peak detection.
}

\usage{
  peaks.fsa(x, ranges, logTransform = FALSE, lowThreshold = 1000, channels = "6-FAM",
    noiseRange = c(-10, 0), primerRange = c(35, 45))
}

\arguments{
  \item{x}{
    An aligned object of class \code{fsa}, as returned by \code{\link{align.fsa}}.
}
  \item{ranges}{
    A named list of ranges, numeric vectors of length two (minimal and maximal values). They defines the sizes (in base pairs) for which a maximum is required.
}
  \item{logTransform}{
    Single logical value, whether to apply log transformation (base 2) to normalized values (previously floored to 0 and summed with 1) or not.
}
  \item{lowThreshold}{
    Single numeric value, threshold for which "low profile" warnings are called if all peaks are lower.
}
  \item{channels}{
    Single character value, the name of the \code{x} channel to browse.
}
  \item{noiseRange}{
    Numeric vector of length 2, defining the range (relative to the starting range of the first peak defined in \code{ranges}) in which measure the noise (in bp). If the noise peak is 20 percent greater than the first peak, a warning is raised as the accuracy of the measure may be compromised.
}
  \item{primerRange}{
    Numeric vector of length 2, defining the range in which measure primer signals (in bp). This is implemented for QC experimentation and may not be useful in current practice.
}
}

\value{
  Returns a \code{data.frame} with a row for each range :
  \item{gene}{The name of the range described, extracted from \code{ranges}.}
  \item{size.min}{Minimal size of the range described.}
  \item{size.max}{Maximal size of the range described.}
  \item{peak.size}{Size at which the maximum was found, in base pairs.}
  \item{peak.height}{Maximum found, in fluorescence units.}
  \item{peak.offScale}{Is there any off-scale value in the range ?}
  \item{normalized}{Current peak's height divided by the mean of all peak heights.}
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{GEP.process}}
}

\examples{
  # Example FSA file provided
  fsa <- read.fsa(system.file("extdata/fsa/A5918.fsa", package="MLPA"))
  fsa <- align.fsa(fsa)
  
  # Single interval
  print(peaks.fsa(fsa, ranges=list(IRF4=c(86.2, 87.5))))
  
  # Using a design file
  design <- designFile(system.file("extdata/design.conf", package="MLPA"))
  print(peaks.fsa(fsa, ranges=design$PEAKS$ranges))
}

