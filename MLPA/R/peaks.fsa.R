# Get maximum over windows
# Author : Sylvain Mareschal <maressyl@gmail.com>
peaks.fsa <- function(
		x,
		ranges,
		logTransform = FALSE,
		lowThreshold = 1000,
		channels = "6-FAM",
		noiseRange = c(-10, 0),
		primerRange = c(35, 45)
		)
	{
	# Ranges check
	if(!is(x, "fsa"))                    stop("'x' must be an 'fsa' object")
	if(!is.list(ranges))                 stop("'ranges' must be a list or a deisgn file name")
	if(length(ranges) == 0)              stop("'ranges' must not be empty")
	if(any(sapply(ranges, length) != 2)) stop("'ranges' must contain length 2 vectors")
	if(any(!sapply(ranges, is.numeric))) stop("'ranges' must contain numeric vectors")
	
	# Points sizes (in bp)
	if(is.null(attr(x, "ladderModel"))) stop("Can't use 'bp' unit without alignment ('ladderModel' attribute)")
	size <- attr(x, "ladderModel")[2] * 1:nrow(x) + attr(x, "ladderModel")[1]
	
	# Identify first peak
	firstPeak <- which.min(sapply(ranges, "[", 1))
	
	# Add primer peak
	ranges$"<primers>" <- primerRange
	
	# Channel recycling
	channels <- rep(channels, length.out=length(ranges))
	names(channels) <- names(ranges)
	
	out <- NULL
	for(geneName in names(ranges)) {
		# Scans in range
		w <- size >= ranges[[geneName]][1] & size <= ranges[[geneName]][2]
		if(!any(w)) stop("No data point in range for ", geneName, " (check alignment)")
		y <- x[ w , channels[geneName] ]
		
		# Off-scale value in range (any channel)
		isOffScale <- any(attr(x, "offScale") %in% which(w))
		
		# Data collection
		out <- rbind(out,
			data.frame(
				gene = geneName,
				size.min = ranges[[geneName]][1],
				size.max = ranges[[geneName]][2],
				peak.size = size[w][ which.max(y) ],
				peak.height = max(y),
				peak.offScale = isOffScale,
				stringsAsFactors = FALSE
			)
		)
	}
	
	# First peak warning
	w <- size >= ranges[[firstPeak]][1] + noiseRange[1] & size <= ranges[[firstPeak]][1] + noiseRange[2]
	y <- x[ w , channels[geneName] ]
	firstPeakValue <- out[ out$gene == names(firstPeak) , "peak.height" ]
	if(firstPeakValue < max(y) * 1.2) warning("First peak value compromised (value: ", round(firstPeakValue, 1), ", noise: ", round(max(y), 1), ")")
	
	# Warn about low values
	maxValue <- max(out[ !grepl("^<.+>$", out$gene) , "peak.height" ])
	if(maxValue < lowThreshold) warning("Low profile (maximal height: ", round(maxValue, 1), ")")
	
	# Warn about off-scale values
	isOffScale <- out$peak.offScale & !grepl("^<.+>$", out$gene)
	if(any(isOffScale)) warning(sum(isOffScale), " off-scale measures (", paste(out[ out$peak.offScale , "gene" ], collapse=", "), ")")
	
	# Normalize
	out$normalized <- out$peak.height / mean(out[ !grepl("^<.+>$", out$gene) , "peak.height" ])
	
	# log-transformation
	if(isTRUE(logTransform)) {
		out[ out$normalized < 0 , "normalized" ] <- 0
		out$normalized <- log(1 + out$normalized, 2)
	}
	
	return(out)
}

