# Train a model object for classify()
# Author : Sylvain Mareschal <maressyl@gmail.com>
train <- function(
		peakMatrix,
		group,
		filter.p = 0.05,
		...
		)
	{
	# Coercions
	if(!is.factor(group)) group <- factor(group)
	
	# Checks
	if(!is.matrix(peakMatrix) || !is.numeric(peakMatrix)) stop("'peakMatrix' must be a numeric matrix")
	if(length(levels(group)) != 2L)                       stop("'group' must be a two-level factor")
	if(any(table(group) == 0L))                           stop("All two levels of 'group' must be represented")
	if(nrow(peakMatrix) != length(group))                 stop("'peakMatrix' must have as many rows as 'group' has values")
	
	# Initialize
	geneNames <- colnames(peakMatrix)
	geneTs <- rep(as.numeric(NA), ncol(peakMatrix))
	geneMs <- rep(as.numeric(NA), ncol(peakMatrix))
	genePs <- rep(as.numeric(NA), ncol(peakMatrix))
	
	# Test all genes
	for(i in 1:ncol(peakMatrix)) {
		# Perform a t-test
		test <- t.test(peakMatrix[,i] ~ group)
		
		# Collect statistics
		geneTs[i] <- test$statistic
		genePs[i] <- test$p.value
		geneMs[i] <- mean(peakMatrix[,i])
	}
	
	# Keep only significant genes
	if(!is.na(filter.p)) {
		keep <- genePs < filter.p
		geneNames <- geneNames[keep]
		geneTs <- geneTs[keep]
		geneMs <- geneMs[keep]
		genePs <- genePs[keep]
		peakMatrix <- peakMatrix[,keep]
	}

	# Compute score
	score <- apply(geneTs * (t(peakMatrix) / geneMs), 2, sum)
	groupMeans <- as.double(tapply(X=score, INDEX=group, FUN=mean))
	groupSDs <- as.double(tapply(X=score, INDEX=group, FUN=sd))
	groupNames <- levels(group)
	
	# Gather as a model object
	out <- model(
		groupMeans = groupMeans,
		groupSDs = groupSDs,
		groupNames = groupNames,
		geneNames = geneNames,
		geneTs = geneTs,
		geneMs = geneMs,
		...
	)
	
	return(out)
}

