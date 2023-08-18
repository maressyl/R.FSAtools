\name{plot.fsaModel}
\alias{plot.fsaModel}
\title{
  Plot method for "fsaModel" objects
}

\description{
  Plots a \code{fsaModel} object.
}

\usage{
\method{plot}{fsaModel}(x, xlab = "Score", lwd = 3, ...)
}

\arguments{
  \item{x}{
    The \code{fsaModel} object to plot.
}
  \item{xlab}{
    To be passed to \code{\link{plot}}.
}
  \item{lwd}{
    To be passed to \code{\link{plot}}.
}
  \item{\dots}{
    Further arguments to be passed to \code{\link{plot}}.
}
}

\value{
  Invisibly returns \code{TRUE} on success.
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{train}}
}

\examples{
  # Example FSA file provided
  fsa <- read.fsa(system.file("extdata/fsa_GEP/A5918.fsa", package="FSAtools"))
  
  # Add model from design file
  design <- designFile(system.file("extdata/design_GEP.conf", package="FSAtools"))
  fsa <- add.model(fsa, model=design$GLOBALS$MODEL)
  
  # Plot model
  plot(attr(fsa, "model"))
}

