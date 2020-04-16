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

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{train}}
}

\examples{
  # Build model from design file
  design <- designFile(system.file("extdata/design.conf", package="MLPA"))
  design$model$disable <- NULL
  model <- do.call("model", design$model)
  
  # Plot model
  plot(model)
}

