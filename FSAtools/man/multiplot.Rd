\name{multiplot}
\alias{multiplot}
\title{
  Wrapper to layout
}

\description{
  Calls \code{\link{layout}} using atomic arguments and convenient defaults, mainly to be included in \code{\link{generic.process}} design files.
}

\usage{
  multiplot(nrow, ncol, widths = rep.int(1, ncol), heights = rep.int(1, nrow),
    indexes = 1:(nrow * ncol), byrow = FALSE, respect = FALSE, cex = 1)
}

\arguments{
  \item{nrow}{
    Single integer value, the amount of rows into which divide the screen for multiple plots.
}
  \item{ncol}{
    Single integer value, the amount of columns into which divide the screen for multiple plots.
}
  \item{widths}{
    To be passed to \code{\link{layout}}.
}
  \item{heights}{
    To be passed to \code{\link{layout}}.
}
  \item{indexes}{
    Integer vector, the ordering of plots while building the \code{mat} matrix for \code{\link{layout}}.
}
  \item{byrow}{
    Single logical value, whether to fill the \code{mat} matrix for \code{\link{layout}} with \code{indexes} by row or by column.
}
  \item{respect}{
    To be passed to \code{\link{layout}}.
}
  \item{cex}{
    Single numeric value, if not \code{NA} \code{\link{par}} will be called to force this value, as large layout automatically change it.
}
}

\value{
  Invisibly returns \code{TRUE} on success.
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{generic.process}}
}

\examples{
  multiplot(nrow=1, ncol=3, widths=c(1,2,2))
  plot(1:5)
  plot(1:5)
  plot(1:5)
}

