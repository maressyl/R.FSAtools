# Build a model object for classify()
# Author : Sylvain Mareschal <maressyl@gmail.com>
model <- function(
			groupMeans,
			groupSDs,
			groupNames,
			groupColors = c("blue", "red"),
			threshold	= 0.9,
			geneNames,
			geneTs,
			geneMs
		)
	{
	object <- list(
		groupMeans = groupMeans,
		groupSDs = groupSDs,
		groupNames = groupNames,
		groupColors = groupColors,
		threshold	= threshold,
		geneNames = geneNames,
		geneTs = geneTs,
		geneMs = geneMs
	)
	class(object) <- "fsaModel"
	return(object)	
}

