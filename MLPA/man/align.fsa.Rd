\name{align.fsa}
\alias{align.fsa}
\title{
  Aligns peaks using size ladder
}

\description{
  This function adds to a \code{fsa} object a linear regression model allowing the raw time indexes to be converted into base pair sizes, using a known size markers ladder.
}

\usage{
  align.fsa(x, channel = "ROX", fullLadder = c(50, 60, 90, 100, 120, 150, 160, 180, 190,
    200, 220, 240, 260, 280, 290, 300, 320, 340, 360, 380, 400), useLadder = c(50, 60, 90,
    100, 120), outThreshold = 0.15, noiseLevel = 10, surePeaks = 5,
    trim = c("forward", "backward", "none"), maskOffScale = FALSE, rMin = 0.999,
    rescue = FALSE, ylim = NULL, ...)
}

\arguments{
  \item{x}{
    An object of class \code{fsa}, as returned by \code{\link{read.fsa}}
}
  \item{channel}{
    Single character value, the name of the channel used for size markers.
}
  \item{fullLadder}{
    Integer vector, the size markers used in the assay (in base pairs).
}
  \item{useLadder}{
    Integer vector, the size markers to use for the alignment (using only size markers nearing the expected size for the experimental peaks may achieve a more precise alignment). They must be present in \code{fullLadder}. If \code{NULL}, \code{fullLadder} will be used entirely.
}
  \item{outThreshold}{
    Single numeric value, maximal distance from the computed size-marker intensity for a peak to be considered as a size-marker. If lower than 1, it is considered as a proportion of the size-marker intensity computed from sure peaks.
}
  \item{noiseLevel}{
    Single numeric value, minimal intensity for a local maximum to be considered as a peak.
}
  \item{maskOffScale}{
    Single logical value, whether to mask indexes with off-scale values in any channel to limit side-effects or not.
}
  \item{surePeaks}{
    Single integer value, amount of peaks to use to compute size-marker intensity. They are selected at the end of the profile, as most artefacts are observed ate the beginning. Consider to reduce this value if your assay was prematurely ended.
}
  \item{trim}{
    Single character value, defining how to behave when more/less peaks than expected are read. "forward" will keep first peaks and adjust discarding the last ones, "backward" will keep last peaks and adjust discarding the first ones, and "none" will generate an error.
}
  \item{rMin}{
    Single numeric value, minimum adjusted r squared value (see \code{\link{summary.lm}}) to consider an alignment as "good". Poor alignments raise a warning, and may be due to artefactual peaks in the size-marker channel or errors in \code{fullLadder} definition. Consider lowering \code{outThreshold} and raising \code{noiseLevel} to minimize artefact selection.
}
  \item{rescue}{
    Single logical value, whether to plot a "rescue" profile or not. Rescue profiles are calls to \code{\link{plot.fsa}} on which diverse additionnal data is drawn to help diagnose alignment problems.
}
  \item{ylim}{
    To be passed to \code{\link{plot.fsa}} for the alignment rescue plot, if enabled (see \code{rescue}).
}
  \item{\dots}{
    Further arguments to be passed to \code{\link{plot.fsa}} for the alignment rescue plot, if enabled (see \code{rescue}).
}
}

\value{
  Returns the object of class \code{fsa} provided with updated \code{\link{attributes}} :
  \item{ladderModel}{A numeric vector of linear regression coefficients to use to convert raw indexes into base pairs.}
  \item{ladderExact}{A named numeric vector of raw indexes at which size markers were detected.}
}

\author{
  Sylvain Mareschal, Philippe Ruminy
}

\examples{
  # Example FSA file provided
  fsa <- read.fsa(system.file("extdata/fsa/A5918.fsa", package="MLPA"))
  
  # Plot subset of the profile (time index)
  plot(fsa, units="index", xlim=c(4000,5000))
  
  # Align using full ladder
  fullLadder <- c(
    50, 60, 90, 100, 120, 150, 160, 180, 190, 200, 220,
    240, 260, 280, 290, 300, 320, 340, 360, 380, 400
  )
  fsa <- align.fsa(fsa, fullLadder=fullLadder)
  
  # Plot subset of the profile (base pairs)
  plot(fsa, units="bp", xlim=c(80,130))
}
