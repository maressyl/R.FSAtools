\name{filter.fsa}
\alias{filter.fsa}
\title{
  Applies filter() to a "fsa" object
}

\description{
  Replaces the requested column of values by the output of \code{\link{filter}}, possibly after masking values out of a specified index or bp range.
}

\usage{
  filter.fsa(x, channel, ..., from = NA, to = NA, units = "bp")
}

\arguments{
  \item{x}{
    An object of class \code{fsa}, as returned by \code{\link{read.fsa}}
}
  \item{channel}{
    Single character value, the name of the channel used for size markers.
}
  \item{\dots}{
    Further arguments to be passed to \code{\link{filter}}.
}
  \item{from}{
    Single numeric value, the starting offset (integer index or numeric bp) to consider. No subsetting will be applied if \code{NA}.
}
  \item{to}{
    Single numeric value, the last offset (integer index or numeric bp) to consider. No subsetting will be applied if \code{NA}.
}
  \item{units}{
    Either "index" or "bp", defining the unit of \code{from} and \code{to}. Notice \code{x} must have been processed by \code{\link{align.fsa}} to use "bp".
}
}

\value{
  Returns \code{x}, with updated content.
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{read.fsa}}
}

\examples{
  # Example FSA file provided
  fsa <- read.fsa(system.file("extdata/fsa_GEP/A5918.fsa", package="FSAtools"))
  fsa <- align.fsa(fsa)
  
  # Profile before filtering
  plot(fsa)
  
  # Plot subset of the profile (base pairs)
  fsa <- filter.fsa(fsa, channel="ROX", filter=20, from=40, to=140, units="bp")
  
  # Profile after filtering
  plot(fsa)
}
