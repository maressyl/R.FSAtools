# Add a LPS model data to a fsa object
add.model <- function(
			x,
			model,
			groupMeans,
			groupSDs,
			groupNames,
			groupColors,
			threshold,
			geneNames,
			geneTs,
			geneMs
		)
	{
	# Checks
	if(!is(x, "fsa")) stop("'x' must be a 'fsa' object")
	
	# Model parameters
	if(missing(model)) {
		# Create list from elements
		model <- list(
			groupMeans = groupMeans,
			groupSDs = groupSDs,
			groupNames = groupNames,
			groupColors = groupColors,
			threshold = threshold,
			geneNames = geneNames,
			geneTs = geneTs,
			geneMs = geneMs
		)
	} else {
		# Check
		if(!is.list(model)) stop("'model' must be a list")
		
		# Overwrite
		if(!missing(groupMeans))  model$groupMeans <- groupMeans
		if(!missing(groupSDs))    model$groupSDs <- groupSDs
		if(!missing(groupNames))  model$groupNames <- groupNames
		if(!missing(groupColors)) model$groupColors <- groupColors
		if(!missing(threshold))   model$threshold <- threshold
		if(!missing(geneNames))   model$geneNames <- geneNames
		if(!missing(geneTs))      model$geneTs <- geneTs
		if(!missing(geneMs))      model$geneMs <- geneMs
		
		# Check
		names <- c("groupMeans", "groupSDs", "groupNames", "groupColors", "threshold", "geneNames", "geneTs", "geneMs")
		if(!setequal(names(model), names)) stop("'model' must have the following components : ", paste(names, collapse=", "))
	}
	
	# Store in object
	class(model) <- "fsaModel"
	attr(x, "model") <- model
	
	return(x)	
}

