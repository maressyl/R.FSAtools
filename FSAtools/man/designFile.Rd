\name{designFile}
\alias{designFile}
\title{
  Process interface's design file
}

\description{
  This function is a slave for \code{\link{generic.process}}. It process a design file and returns its processed elements as a list.
}

\usage{
  designFile(fileName)
}

\arguments{
  \item{fileName}{
    Single character value, the path and name of a design file to process.
}
}
\details{
  Design files are text file split in multiple sections. Each section starts with a "[NAME]" line and ends when the next section begins. Lines starting with a \code{#} sign are ignored, as well as blank lines.
  
  Standard sections refer to an existing R function (the section name is expected to match the function name in a case-sensitive manner), each line in the section setting an argument to call this function : the first value is the argument name, then after a tabulation come one or many values separated by tabulations. Multiple values will be aggregated into a vector, and \code{type.convert} will try to guess the correct type.
  
  One or many modifiers can be added in the section name after the function name followed by the \code{:} sign and separated by commas (e.g. \code{file.remove:first,nowarn}). Modifiers have the following effects :
  \describe{
	\item{first}{The function will be called only once, while processing the first .fsa file.}
	\item{last}{The function will be called only once, while processing the last .fsa file.}
	\item{nowarn}{The function will be embedded inside \code{suppressWarnings} to silently ignore warnings.}
  }
  
  \code{$NAME} and \code{@NAME} can be used while setting an argument to refer to global variable \code{NAME}. Use \code{$} reference for character variables, which will be replaced with \code{gsub} and can thus be combined (e.g. \code{$NAME.txt}). Use \code{@} reference to obtain the raw R variable, regardless of its type. Default globals are :
  \describe{
	\item{FILE_PATH}{The full name and path of the .fsa file currently considered in the loop.}
	\item{FILE_DIR}{The parent directory of the .fsa file currently considered in the loop.}
	\item{FILE_NAME}{The base name (without path) of the .fsa file currently considered in the loop.}
	\item{OBJECT}{The last object of class 'fsa' returned by any function called in the pipeline.}
	\item{OUTPUT_PATH}{The full name and path defined by \code{output} while calling \code{\link{generic.process}}.}
	\item{OUTPUT_DIR}{The parent directory of \code{output} while calling \code{\link{generic.process}}.}
	\item{OUTPUT_NAME}{The base name (without path) defined by \code{output} while calling \code{\link{generic.process}}.}
  }
  
  Sections with full uppercase names (only letters and \code{_} are allowed) will define a new global variable with matching name. The global will be a named list, each line in the section defining a vector (the first value being the name for the vector in the list). The \code{table} modifier can be used while defining globals (e.g. \code{PEAKS:table}), to request instead the section to be parsed as a TSV file (one row in the section is one row in the table, columns are separated by tabulations, see \code{\link{read.table}}). The first row will be used as column names and the first column as row names.
  
  A DESIGN global (\code{[DESIGN]} section) is strongly recommended to keep design file self-explained, the following elements are suggested :
  \describe{
    \item{author}{The name of the design author (for human readers only).}
    \item{purpose}{The description of the design (for human readers only).}
    \item{FSAtools}{Version of the FSAtools package for which the design was created (separated with dots).}
    \item{updated}{Date of the last design update (YYYY-MM-DD).}
  }
  
  Please refer to the two provided working examples to help building your own designs. Full processing with the examples are described in \code{\link{generic.process}} \code{examples} section, direct access to the example design files is shown below.
}

\value{
  Returns a multi-level \code{list}, with a direct children per function to call and an extra \code{GLOBALS} element.
  
  Children are named according to the function to call, thus multiple children can have the same name.
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{generic.process}}, \code{\link{generic.interface}}
}

\examples{
  # Example file provided
  file <- system.file("extdata/design_GEP.conf", package="FSAtools")
  design <- designFile(file)
  
  # Alignment rescue design provided
  file <- system.file("extdata/design_SNP.conf", package="FSAtools")
  design <- designFile(file)
}

