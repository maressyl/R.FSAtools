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

\value{
  Invisibly returns \code{TRUE} on success.
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
  print(fsa)
  
  # Aligned version
  fsa <- align.fsa(fsa)
  print(fsa)
}
