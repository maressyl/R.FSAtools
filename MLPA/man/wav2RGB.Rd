\name{wav2RGB}
\alias{wav2RGB}
\title{
  Converts light wavelengths to RGB colors
}

\description{
  Converts wavelengths in nanometers into corresponding visible colors.
}

\usage{
  wav2RGB(wav)
}

\arguments{
  \item{wav}{
    Numeric vector of wavelengths (in nanometers) to convert into colors.
}
}

\value{
  Returns a character vector of the same length as \code{wav}, with an RGB color for each wavelength. Wavelengths out of visible ranges return black.
}

\references{
  \url{http://codingmess.blogspot.fr/2009/05/conversion-of-wavelength-in-nanometers.html}
}

\author{
  Sylvain Mareschal
}

\examples{
  wv <- seq(from=300, to=800, by=10)
  plot(x=wv, y=rep(1, length(wv)), col=wav2RGB(wv), pch=19)
}
