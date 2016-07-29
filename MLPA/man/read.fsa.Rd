\name{read.fsa}
\alias{read.fsa}
\title{
  Imports a .fsa file from Applied Biosystems
}

\description{
  This function parses a FSA file holding fragment analysis data, using \pkg{seqinr} package's \code{\link{read.abif}}.
}

\usage{
  read.fsa(file, applyLowess = TRUE, processed = FALSE, ...)
}

\arguments{
  \item{file}{
    Single character value, the name and path of the file to parse.
}
  \item{applyLowess}{
    Single logical value, whether to apply \code{\link{lowess}} on intensities to smooth time-related biases or not.
}
  \item{processed}{
    Single logical value, whether to use processed DATA values (as stored in sets 9 to 12, not always available) rather than raw values (sets 1 to 4). If \code{NA}, processed ones will be used as long as they are available, else raw ones will be used instead.
}
  \item{\dots}{
    Further arguments to be passed to \code{\link{read.abif}}.
}
}

\value{
  A S3 object of class \code{fsa}
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{read.abif}}, \code{\link{MLPA.process}}, \code{\link{plot.fsa}}
}

\examples{
  # Example FSA file provided
  fsa <- read.fsa(system.file("extdata/fsa/A5918.fsa", package="MLPA"))
  print(fsa)
}
