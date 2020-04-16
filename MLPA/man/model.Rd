\name{model}
\alias{model}
\title{
  Object constructor for binary predictors
}

\description{
  This function aggregates the data required to predict class in \code{\link{classify}}.
}

\usage{
  model(groupMeans, groupSDs, groupNames, groupColors = c("blue", "red"),
    threshold = 0.9, geneNames, geneTs, geneMs)
}

\arguments{
  \item{groupMeans}{
    Numeric vector of length 2, the means of the scores in each group as computed on a training series.
}
  \item{groupSDs}{
    Numeric vector of length 2, the standard deviations of the scores in each group as computed on a training series.
}
  \item{groupNames}{
    Character vector of length 2, the names of the group described in \code{groupMeans}, \code{groupSDs} and \code{groupColors}.
}
  \item{groupColors}{
    Character vector of length 2, the colors to use to plot each group (see \code{\link{par}} for allowed values).
}
  \item{threshold}{
    Single numeric vector, the confidence threshold to use for prediction (a call will be made only if it is at least at this level of certainty).
}
  \item{geneNames}{
    Character vector, the names of the genes whose expression is to be used.
}
  \item{geneTs}{
    Numeric vector, for each gene in \code{geneNames}, the statistic of a \code{\link{t.test}} comparing its expression between the two groups in a training series.
}
  \item{geneMs}{
    Numeric vector, for each gene in \code{geneNames}, the mean expression in the whole training series.
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
  # Build from design file
  design <- designFile(system.file("extdata/design.conf", package="MLPA"))
  design$model$disable <- NULL
  model <- do.call("model", design$model)
  
  # Observe model
  print(model)
  plot(model)
}
