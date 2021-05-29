\name{generic.process}
\alias{generic.process}
\alias{generic.interface}
\title{
  Processing multiple FSA files
}
\description{
  \code{generic.process} handles the whole analysis of a series of .fsa files according to the pipeline described in the user-provided design file, generating tabular and graphical profiles.
  
  \code{generic.interface} summons a Tcl-Tk interface to call \code{generic.process} interactively.
}

\usage{
  generic.process(input, design, output, include = NULL, exclude = NULL,
     progressBar = NULL)
  generic.interface()
}

\arguments{
  \item{input}{
    Single character value, the path to a directory containing .fsa files to analyse. Notice it will be explored recursively, so sub-directories are allowed.
}
  \item{design}{
    Single character value, the path to a design file, as handled by \code{\link{designFile}}.
}
  \item{output}{
    Single character value, the path to a ".pdf" or ".log" file that will be created during the analysis.
}
  \item{include}{
    Single character value, a regular expression files (with relative path) in \code{input} must match to be processed (ignored if \code{NULL}).
}
  \item{exclude}{
    Single character value, a regular expression files (with relative path) in \code{input} must not match to be processed (ignored if \code{NULL}).
}
  \item{progressBar}{
    A \code{ttkprogressbar} to increment during the processing, or \code{NULL}. This argument is only provided to connect \code{GEP.interface} and \code{GEP.process}, thus it should be ignored.
}
}

\details{
  The content of the analysis pipeline is fully controlled by the design file, see \code{\link{designFile}} for details and the \code{examples} section below for two working examples provided in the package.
  
  More generally, \code{generic.process} loops over the list of .fsa files in the \code{input} directory and calls the requested functions one after the other, updating the \code{fsa} object at each step.
}

\value{
  Return nothing. \code{generic.process} raise errors, warnings and messages which are intercepted by \code{generic.interface} and redirected to the log file (\code{output}.log).
}

\references{
  Mareschal, Ruminy et al (2015) <doi:10.1016/j.jmoldx.2015.01.007> "Accurate Classification of Germinal Center B-Cell-Like/Activated B-Cell-Like Diffuse Large B-Cell Lymphoma Using a Simple and Rapid Reverse Transcriptase-Multiplex Ligation-Dependent Probe Amplification Assay: A CALYM Study"
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{designFile}}
}

\examples{
  ### EXAMPLE 1 : Gene expression (RT-MLPA) ###
  
  # Working in temporary directory
  output <- sprintf("\%s/GEP", tempdir())
  
  # See files before analysis
  dir(system.file("extdata", package="FSAtools"))
  
  # Launch analysis in package directory
  generic.process(
    input = system.file("extdata/fsa_GEP", package="FSAtools"),
    design = system.file("extdata/design_GEP.conf", package="FSAtools"),
    output = output
  )
  
  # List resulting files
  dir(dirname(output), full.names=TRUE)
  
  
  ### EXAMPLE 2 : Genotyping ###
  
  # Working in temporary directory
  output <- sprintf("\%s/SNP", tempdir())
  
  # See files before analysis
  dir(system.file("extdata", package="FSAtools"))
  
  # Launch analysis in package directory
  generic.process(
    input = system.file("extdata/fsa_SNP", package="FSAtools"),
    design = system.file("extdata/design_SNP.conf", package="FSAtools"),
    output = output
  )
  
  # List resulting files
  dir(dirname(output), full.names=TRUE)
}
