\name{train}
\alias{train}
\title{
  Training function for binary predictors
}

\description{
  This function build a model from data to predict class in \code{\link{classify}}.
}

\usage{
  train(peakMatrix, group, filter.p = 0.05, ...)
}

\arguments{
  \item{peakMatrix}{
    Numeric matrix of normalized peak heights with samples in rows and peaks in columns.
}
  \item{group}{
    Two-level factor defining the group of every samples in \code{peaks}.
}
  \item{filter.p}{
    Single numeric value, if not \code{NA} only genes for which the t-test p is lower than this will be used in the model.
}
  \item{\dots}{
    Further arguments to be passed to \code{\link{model}}.
}
}

\value{
  Returns an S3 object of class \code{fsaModel}.
}

\author{
  Sylvain Mareschal
}

\seealso{
  \code{\link{model}}, \code{\link{classify}}
}

