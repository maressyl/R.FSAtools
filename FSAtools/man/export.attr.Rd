\name{export.attr}
\alias{export.attr}
\title{
  Print an attribute of a 'fsa' object to a file
}

\description{
  Adds the content of an attribute of a 'fsa' object to a CSV file, either appending new rows or columns.
  
  Typically used in \code{\link{generic.process}} via the design file to export numeric data during the processing.
}

\usage{
  export.attr(x, attr, file, meta = character(0), sep = "\t", dec = ".", quote = TRUE)
}

\arguments{
  \item{x}{
    The \code{fsa} object whose attribute is to be printed.
}
  \item{attr}{
    Single character value, the name of the attribute to print.
}
  \item{file}{
    Single character value, the path and name to the file to create or update.
}
  \item{meta}{
    Character vector, the names of \code{x} meta-data fields to export as extra columns.
}
  \item{sep}{
    To be passed to \code{\link{write.table}}.
}
  \item{dec}{
    To be passed to \code{\link{write.table}}.
}
  \item{quote}{
    To be passed to \code{\link{write.table}}.
}
}

\details{
  New data will be appended below the content of the file with extra 'meta' columns, row names will be added only if the file was empty.
}

\value{
  Invisibly returns \code{TRUE} on success.
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{read.fsa}}, \code{\link{peaks.fsa}}, \code{\link{genotype.closest.fsa}}, \code{\link{genotype.ratio.fsa}}, \code{\link{classify}}
}

