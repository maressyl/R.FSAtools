\name{print.fsa}
\alias{print.fsa}
\title{
  Print method for "fsa" objects
}

\description{
  Prints a short summary of an \code{fsa} object.
}

\usage{
\method{print}{fsa}(x, ...)
}

\arguments{
  \item{x}{
    The \code{fsa} object to print.
}
  \item{\dots}{
    Currently ignored.
}
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{read.fsa}}
}

\examples{
  # Example FSA file provided
  fsa <- read.fsa(system.file("extdata/fsa/A5918.fsa", package="MLPA"))
  print(fsa)
  
  # Aligned version
  fsa <- align.fsa(fsa)
  print(fsa)
}
