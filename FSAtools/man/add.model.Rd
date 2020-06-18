\name{add.model}
\alias{add.model}
\title{
  Object constructor for binary predictors
}

\description{
  This function aggregates the data required to predict class in \code{\link{classify}}.
}

\usage{
  add.model(x, model, groupMeans, groupSDs, groupNames, groupColors, threshold,
    geneNames, geneTs, geneMs)
}

\arguments{
  \item{x}{
    An object of class \code{fsa} to which add a model.
}
  \item{model}{
    A \code{list} aggregating the model parameters. Can be provided instead of the argument vectors below individually.
}
  \item{groupMeans}{
    Numeric vector of length 2, the means of the scores in each group as computed on a training series. If \code{model} is provided, this vector will be used instead.
}
  \item{groupSDs}{
    Numeric vector of length 2, the standard deviations of the scores in each group as computed on a training series. If \code{model} is provided, this vector will be used instead.
}
  \item{groupNames}{
    Character vector of length 2, the names of the group described in \code{groupMeans}, \code{groupSDs} and \code{groupColors}. If \code{model} is provided, this vector will be used instead.
}
  \item{groupColors}{
    Character vector of length 2, the colors to use to plot each group (see \code{\link{par}} for allowed values). If \code{model} is provided, this vector will be used instead.
}
  \item{threshold}{
    Single numeric vector, the confidence threshold to use for prediction (a call will be made only if it is at least at this level of certainty). If \code{model} is provided, this vector will be used instead.
}
  \item{geneNames}{
    Character vector, the names of the genes whose expression is to be used. If \code{model} is provided, this vector will be used instead.
}
  \item{geneTs}{
    Numeric vector, for each gene in \code{geneNames}, the statistic of a \code{\link{t.test}} comparing its expression between the two groups in a training series. If \code{model} is provided, this vector will be used instead.
}
  \item{geneMs}{
    Numeric vector, for each gene in \code{geneNames}, the mean expression in the whole training series. If \code{model} is provided, this vector will be used instead.
}
}

\value{
  Returns an S3 object of class \code{fsaModel}.
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{classify}}
}

\examples{
  # Example FSA file provided
  fsa <- read.fsa(system.file("extdata/fsa_GEP/A5918.fsa", package="FSAtools"))
  
  # Add model from design file
  design <- designFile(system.file("extdata/design_GEP.conf", package="FSAtools"))
  fsa <- add.model(fsa, model=design$GLOBALS$MODEL)
  
  # Observe model
  print(attr(fsa, "model"))
  plot(attr(fsa, "model"))
}
