\name{plot.fsa}
\alias{plot.fsa}
\title{
  Plot method for "fsa" objects
}

\description{
  Plots a \code{fsa} object. For each selected channel, a line is drawn bewteen measured fluorescence intensities (y axis) along the electrophoresis time (x axis).
}

\usage{
\method{plot}{fsa}(x, units = NA, channels = NA, chanColors = NA, ladder = TRUE,
    offScaleCol = "#FF0000", offScalePch = "+", offScaleCex = 0.4, bg = "white",
    fg = "black", title = "", xlab = NA, ylab = "Intensity", xlim = NA, ylim = NA,
    xaxt = "s", yaxt = "s", bty = "o", xaxp = NA, nticks = 5, all.bp = TRUE, ...)
}

\arguments{
  \item{x}{
    The \code{fsa} object to plot.
}
  \item{units}{
    Single character value, the unit to use on x axis. "index" uses the raw index contained in files, "bp" usess base pair estimations but needs the object to be aligned first using \code{align.fsa}. \code{NA} will select "bp" if \code{x} is aligned, "index" elsewhere.
}
  \item{channels}{
    Character or integer vector, the channels to plot. If \code{NA}, all channels are selected.
}
  \item{chanColors}{
    Character vector defining colors to use to plot channels. Can be named according to channel names stored in \code{x}, or parallel with \code{channels} (first color for first channel, etc, no recycling). If \code{NA}, colors stored in \code{x} are used. See the \code{col} argument in \code{\link{par}} for further details on allowed values.
}
  \item{ladder}{
    Single logical value, whether to add an x axis with size ladder peaks or not. Raises a warning if \code{x} was not aligned before plotting.
}
  \item{offScaleCol}{
    To be passed to \code{\link{points}} for off-scale value plot (see \code{\link{par}} for allowed values).
}
  \item{offScalePch}{
    To be passed to \code{\link{points}} for off-scale value plot (see \code{\link{par}} for allowed values).
}
  \item{offScaleCex}{
    To be passed to \code{\link{points}} for off-scale value plot (see \code{\link{par}} for allowed values).
}
  \item{bg}{
    See \code{\link{par}} for further details.
}
  \item{fg}{
    See \code{\link{par}} for further details. This value is also used for col.axis, col.lab, col.main and col.sub graphical parameters.
}
  \item{title}{
    Single character value, the main title to print on the plot.
}
  \item{xlab}{
    See \code{\link{plot}} for further details. If \code{NA}, \code{units} is used.
}
  \item{ylab}{
    See \code{\link{plot}} for further details.
}
  \item{xlim}{
    See \code{\link{plot}} for further details. If \code{NA}, \code{x} range is used.
}
  \item{ylim}{
    See \code{\link{plot}} for further details. If \code{NA}, \code{x} range is used.
}
  \item{xaxt}{
    See \code{\link{par}} for further details.
}
  \item{yaxt}{
    See \code{\link{par}} for further details.
}
  \item{bty}{
    See \code{\link{par}} for further details.
}
  \item{xaxp}{
    See \code{\link{par}} for further details. If \code{NA}, a suitable value is computed.
}
  \item{nticks}{
    Single integer value. When \code{xaxp} is \code{NA} and \code{units} is "bp", this values fixes the interval between X axis labels.
}
  \item{all.bp}{
    Single logical value, whether to force an unlabeled axis tick at each bp when \code{units} is "bp" or not.
}
  \item{\dots}{
    Further arguments to be passed to \code{\link{plot}}.
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
  
  # Plot whole profile
  plot(fsa)
  
  # Plot subset of the profile (time index)
  plot(fsa, units="index", xlim=c(4000,5000))
  
  # Plot subset of the profile (base pairs)
  fsa <- align.fsa(fsa)
  plot(fsa, units="bp", xlim=c(80,130))
}

