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
  export.attr(x, attr, file, sample, sep = "\t", dec = ".", quote = TRUE)
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
  \item{sample}{
    Single character value, the sample name to add for disambiguation when appending to a file.
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
  If the attribute is a \code{data.frame}, it will be appended below the content of the file with an extra 'sample' column, row names will be added only if the file was empty.
  
  If the attribute is \code{atomic}, it will be appended as a new row to the existing matrix, using 'sample' as row name.
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{read.fsa}}, \code{\link{peaks.fsa}}, \code{\link{genotype.fsa}}, \code{\link{classify}}
}

