\name{designFile}
\alias{designFile}
\title{
  Process interface's design file
}

\description{
  This function is a slave for \code{\link{GEP.interface}}. It process a design file and returns its processed elements as a list.
}

\usage{
  designFile(fileName, overwrite = list())
}

\arguments{
  \item{fileName}{
    Single character value, the path and name of a design file to process.
}
  \item{overwrite}{
    Two-level named list as returned by \code{designFile}, values from this list will replace parameters parsed from the design file.
}
}
\details{
  Design files are tab-separated files, with a header line and dots as decimal separators. It can be quoted, with double or simple quotes. "#" are considered as comment markup.
  
  Design files are plain text files, separated in multiple sections. Sections begin with a "[NAME]" line and end when the next section begin. Comments (line beggining with the "#" character) and blank lines can be added everywhere, except into the PEAKS table. Quotes should be avoided, and multiple values for a single parameter (vectors) can be obtained by concatenating multiple values separated by commas (no spacing).
  
  The "[DESIGN]" section contain name/value pairs, separated by tabulations :
  \describe{
    \item{author}{The name of the design author (for human readers only).}
    \item{purpose}{The description of the design (for human readers only).}
    \item{MLPA}{Version of the MLPA package for which the design was created (separated with dots).}
    \item{updated}{Date of the last design update (YYYY-MM-DD).}
  }
  
  The "[PEAKS]" section contain a tab-separated table, with a header line and dots as decimal separators. It should contain one row for each peak that is to be measured in the assay.
  \describe{
    \item{name}{Character, the name of the gene described (mist be unique in the table).}
    \item{channel}{Character, the name of the channel in which a measure has to be done.}
    \item{size.min}{Numeric, the start of the range in which a measure has to be done. For size markers, the range is defined in time index units, for genes in base pairs (after an \code{\link{align.fsa}} call).}
    \item{size.max}{Numeric, the end of the range in which a measure has to be done. For size markers, the range is defined in time index units, for genes in base pairs (after an \code{\link{align.fsa}} call).}
    \item{color}{Character, the color tu use for the range on the plot. This can be an english color name, or an hexadecimal color specification as "#000000". Notice transparency will be added.}
  }
  
  Function-specific sections may also be present, to define R function arguments. Currently \code{\link{read.fsa}}, \code{\link{align.fsa}}, \code{\link{plot.fsa}}, \code{\link{model}} and \code{\link{classify}} functions are handled. Section names are supposed to respect function name case, surrounded by square brackets. These sections should contain name/value pairs, separated by tabulations, using only valid function arguments as names (please refer to the corresponding help page). Arguments not defined in the design file will be set to the function's default (raising a warning), and arguments not used by the function will be ignored (raising a warning too). An additionnal \code{disable} argument is handled (single logical value), to disable the call to the corresponding function.
}

\value{
  Returns a multi-level \code{list}, with a direct children per section.
  
  "DESIGN" directly transcribes as a named and typed \code{list} the elements described above.
  
  "PEAKS" contains a list with the following elements :
  \item{ranges}{The concatenated \code{size.min} and \code{size.max} columns, for genes.}
  \item{channels}{The \code{channel} column, for genes.}
  \item{weights}{The \code{weigh} column, for genes.}
  \item{colors}{The plain color from \code{color} column, for genes.}
  \item{backgrounds}{The transparent version of the \code{color} column, for genes.}
  
  Function-specific sections directly transcribe arguments as a named and typed \code{list}.
}

\note{
  It is highly recommended to new users to build their design modifying "extdata/design.conf", updating the "[model]" and "[PEAKS]" sections.
  
  Replacing the content of the "[model]" and "[classify]" sections by "disable TRUE" is also a good way to start with designs, as the classification model is optionnal and can be added later once the technique is developped.
  
  The "extdata/recue.conf" design can prove particularly useful when processing fails with a non-obvious message (generally because the alignment failed). Processing files with this design disables almost everything and provides a raw profile that can help disgnose the problem.
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{GEP.interface}}
}

\examples{
  # Example file provided
  file <- system.file("extdata/design.conf", package="MLPA")
  design <- designFile(file)
  
  # Alignment rescue design provided
  file <- system.file("extdata/rescue.conf", package="MLPA")
  design <- designFile(file)
}

