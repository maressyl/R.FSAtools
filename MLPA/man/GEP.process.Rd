\name{GEP.process}
\alias{GEP.process}
\alias{GEP.interface}
\title{
  MLPA peak detection
}
\description{
  \code{GEP.process} handles the whole analysis from .fsa files, generating tables of expression values and graphical profiles.
  
  \code{GEP.interface} summons a Tcl-Tk interface to call \code{GEP.process} interactively.
}

\usage{
  GEP.process(input, design, output, overwrite = list(), gene.cex = 1.3,
    file.line = 3, mar = c(5,4,5,1), progressBar = NULL)
  GEP.interface()
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
  \item{overwrite}{
    Named list, to be passed to \code{designFile()}.
}
  \item{gene.cex}{
    Single numeric value, the character expanding factor for gene names.
}
  \item{file.line}{
    Single numeric value, the line on which print the file name on plots.
}
  \item{mar}{
    Numeric vector with 4 values, the margin sizes on bottom, left, top and right sides respectively.
}
  \item{progressBar}{
    A \code{ttkprogressbar} to increment during the processing, or \code{NULL}. This argument is only provided to connect \code{GEP.interface} and \code{GEP.process}, thus it should be ignored.
}
}

\value{
  Return nothing. \code{GEP.process} raise errors, warnings and messages which are intercepted by \code{GEP.interface} and redirected to the log file.
  
  Various files are produced, in location set by the \code{output} argument :
  \describe{
    \item{~.expr.tsv}{Numeric matrix of normalized expressions (each sample is divided by its means).}
    \item{~.peaks.tsv}{Table collecting all the peaks called during the analysis, with their size, intensity, sample, annotation and off-scale status.}
    \item{~.pdf}{Profiles of the samples analysed. See \code{\link{designFile}} for customisation.}
    \item{~.log}{Log file of errors, warnings and message (only \code{GEP.interface} produces it).}
  }
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
  # Working in temporary directory
  output <- sprintf("\%s/test.log", tempdir())
  
  # See files before analysis
  dir(system.file("extdata", package="MLPA"))
  
  # Launch analysis in package directory
  GEP.process(
    input = system.file("extdata/fsa", package="MLPA"),
    design = system.file("extdata/design.conf", package="MLPA"),
    output = output
  )
  
  # List resulting files
  dir(dirname(output), full.names=TRUE)
}
