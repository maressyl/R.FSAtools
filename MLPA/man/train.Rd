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

\examples{
  # Underlying truth for pseudo-data (10 genes)
  geneNames <- paste("gene", LETTERS[1:10], sep=".")
  geneMean <- abs(rnorm(10))
  groupShift <- rnorm(10, sd=0.1)
  
  # Generate pseudo-data for 50 samples
  mtx <- NULL
  for(g in 1:10) {
    x <- rnorm(n=50, mean=geneMean[g], sd=0.1)
    x[1:25] <- x[1:25] + groupShift[g]
    x[26:50] <- x[26:50] - groupShift[g]
    mtx <- cbind(mtx, x)
  }
  colnames(mtx) <- geneNames
  rownames(mtx) <- c(
    paste("group1", 1:25, sep="."),
    paste("group2", 26:50, sep=".")
  )
  
  # Train model
  group <- c(
    rep("group1", 25),
    rep("group2", 25)
  )
  model <- train(mtx, group)
  plot(model)
  
  # Compare model to truth
  i <- match(geneNames, model$geneNames)
  out <- data.frame(
    gene = geneNames,
    true.M = geneMean,
    model.M = model$geneMs[i],
    true.shift = groupShift,
    model.T = model$geneTs[i]
  )
  print(out)
}

