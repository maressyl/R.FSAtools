\name{read.fsa}
\alias{read.fsa}
\title{
  Imports a .fsa file from Applied Biosystems
}

\description{
  This function parses a FSA file holding fragment analysis data, using \pkg{seqinr} package's \code{\link{read.abif}}.
}

\usage{
  read.fsa(file, lowess = TRUE, lowess.top = 200, processed = FALSE, meta.extra = NULL,
    quiet = FALSE, ...)
}

\arguments{
  \item{file}{
    Single character value, the name and path of the file to parse.
}
  \item{lowess}{
    Single logical value, whether to apply \code{\link{lowess}} on intensities to smooth time-related biases or not.
}
  \item{lowess.top}{
    Single numeric value, values flagged as "off-scale" or above this threshold will be replaced by \code{lowess.top} to compute the lowess smooth, in order to limit the impact of high and wide peaks.
}
  \item{processed}{
    Single logical value, whether to use processed DATA values (as stored in sets 9 to 12, not always available) rather than raw values (sets 1 to 4). If \code{NA}, processed ones will be used as long as they are available, else raw ones will be used instead.
}
  \item{meta.extra}{
    Named character vector, defining which extra fields to extract to populate the \code{runMetaData} attribute. The vector names define the human-readable names to use in output, the vector values provide the 4 uppercase letter code to extract (all values will be gathered in a vector if the code is used several times). See the reference provided in \code{\link{read.abif}} for existing codes in the ABIF file format.
}
  \item{quiet}{
    Single logical value, whether to print FSA meta-data read from the file or not.
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
  \code{\link{read.abif}}, \code{\link{generic.process}}, \code{\link{plot.fsa}}, \code{\link{read.sanger}}
}

\examples{
  # Example FSA file provided
  fsa <- read.fsa(system.file("extdata/fsa_GEP/A5918.fsa", package="FSAtools"))
  print(fsa)
}

