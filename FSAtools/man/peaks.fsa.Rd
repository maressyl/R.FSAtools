\name{peaks.fsa}
\alias{peaks.fsa}
\title{
  Get maximal value in ranges
}
\description{
  Look for the maximal value in one or many ranges, typically for peak detection.
}

\usage{
  peaks.fsa(x, peaks, names, size.min, size.max, channels, colors,
    logTransform = FALSE, lowThreshold = 1000, noiseRange = c(-10, 0))
}

\arguments{
  \item{x}{
    An aligned object of class \code{fsa}, as returned by \code{\link{align.fsa}}.
}
  \item{peaks}{
    A \code{data.frame} with one row for each peak to consider. Can be provided instead of individual \code{names} (corresponding to \code{peaks} row names instead of a column), \code{size.min}, \code{size.max}, \code{channels} and \code{colors} argument vectors.
}
  \item{names}{
    Character vector, the names to give to the peaks (duplicated values should be avoided). If \code{peaks} is provided, this vector will be used instead.
}
  \item{size.min}{
    Numeric vector, the minimal size (in base pairs) to look for the corresponding peak. If \code{peaks} is provided, this vector will be used instead.
}
  \item{size.max}{
    Numeric vector, the maximal size (in base pairs) to look for the corresponding peak. If \code{peaks} is provided, this vector will be used instead.
}
  \item{channels}{
    Character vector, the name of the channel in \code{x} in which to look for the corresponding peak. If \code{peaks} is provided, this vector will be used instead.
}
  \item{colors}{
    Vector defining the color to use in future plots to highlight the corresponding peak. If \code{peaks} is provided, this vector will be used instead.
}
  \item{logTransform}{
    Single logical value, whether to apply log transformation (base 2) to normalized values (previously floored to 0 and summed with 1) or not.
}
  \item{lowThreshold}{
    Single numeric value, threshold for which "low profile" warnings are called if all peaks are lower.
}
  \item{noiseRange}{
    Numeric vector of length 2, defining the range (relative to the starting range of the first peak defined in \code{ranges}) in which measure the noise (in bp). If the noise peak is 20 percent greater than the first peak, a warning is raised as the accuracy of the measure may be compromised.
}
}

\value{
  Returns \code{x}, with a new or updated \code{peaks} attribute, a \code{data.frame} with a row for each range :
  \item{size.min}{User-provided argument.}
  \item{size.max}{User-provided argument.}
  \item{channels}{User-provided argument.}
  \item{colors}{User-provided argument.}
  \item{size}{Size at which the maximum was found, in base pairs.}
  \item{height}{Maximum found, in fluorescence units.}
  \item{offScale}{Is there any off-scale value in the range ?}
  \item{normalized}{Current peak's height divided by the mean of all peak heights.}
  
  A \code{normalized} vector attribute is also set, corresponding to the \code{normalized} column only.
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{generic.process}}
}

\examples{
  # Example FSA file provided
  fsa <- read.fsa(system.file("extdata/fsa_GEP/A5918.fsa", package="FSAtools"))
  fsa <- align.fsa(fsa)
  
  # Single custom interval
  fsa <- peaks.fsa(
    fsa,
    names = "IRF4",
    size.min = 86.2,
    size.max = 87.5,
    channels = "6-FAM",
    colors = "blue"
  )
  print(attr(fsa, "peaks"))
  
  # Using a design file
  design <- designFile(system.file("extdata/design_GEP.conf", package="FSAtools"))
  fsa <- peaks.fsa(fsa, peaks=design$GLOBALS$PEAKS)
  print(attr(fsa, "peaks"))
}

