\name{generic.log}
\alias{generic.log}
\title{
  Wrapper for \code{generic.process}
}
\description{
  This function is mainly a wrapper for \code{\link{generic.process}}, diverting messages, warnings and errors to a more readable log file.
}

\usage{
  generic.log(..., logFile)
}

\arguments{
  \item{\dots}{
    Arguments to be passed to \code{\link{generic.process}}.
}
  \item{logFile}{
    Single character value, the path and name of the log file where to divert output.
}
}

\value{
  Either an \code{error} object if one occured, an integer number of warnings which happened during the (otherwise successfull) processing or \code{TRUE} if everything went fine.
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{generic.process}}, \code{\link{generic.interface}}
}

\examples{
  # Working in temporary directory
  output <- sprintf("\%s/GEP", tempdir())
  logFile <- sprintf("\%s.log", output)
  
  # Direct analysis
  generic.process(
    input = system.file("extdata/fsa_GEP", package="FSAtools"),
    design = system.file("extdata/design_GEP.conf", package="FSAtools"),
    output = output
  )
  
  # Logged analysis (check logFile)
  generic.log(
    input = system.file("extdata/fsa_GEP", package="FSAtools"),
    design = system.file("extdata/design_GEP.conf", package="FSAtools"),
    output = output,
    logFile = logFile
  )
}
