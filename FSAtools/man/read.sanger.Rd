\name{read.sanger}
\alias{read.sanger}
\title{
  Imports a .ab1 file from Applied Biosystems corresponding to Sanger sequencing
}

\description{
  This function parses a FSA/AB1 file using \code{\link{read.fsa}}, with few adjustments for Sanger sequencing experiments.
}

\usage{
  read.sanger(file, channelOrder = NULL, guess.threshold = 0.3, processed = NA,
     lowess = FALSE, ...)
}

\arguments{
  \item{file}{
    Single character value, the name and path of the file to parse.
}
  \item{channelOrder}{
    Character vector, providing 'A', 'C', 'G' and 'T' in the order of the used channels. If \code{NULL}, a guess will be attempted based on the called sequence.
}
  \item{guess.threshold}{
    Single numeric value, setting the tolerance to use for channel guessing validation. Lower values mean higher chances to get an error for channel guessing failure.
}
  \item{processed}{
    To be passed to \code{\link{read.fsa}}.
}
  \item{lowess}{
    To be passed to \code{\link{read.fsa}}.
}
  \item{\dots}{
    To be passed to \code{\link{read.fsa}}.
}
}

\value{
  A S3 object of class \code{fsa}
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{read.fsa}}, \code{\link{read.abif}}
}

