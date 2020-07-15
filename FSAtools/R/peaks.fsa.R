# Get maximum over windows
peaks.fsa <- function(
		x,
		peaks,
		names,
		size.min,
		size.max,
		channels,
		colors,
		logTransform = FALSE,
		lowThreshold = 1000,
		noiseRange = c(-10, 0)
		)
	{
	# Peak table
	if(missing(peaks)) {
		# Recycle if necessary
		peaks <- data.frame(
			row.names = names,
			size.min = size.min,
			size.max = size.max,
			channel = channels,
			color = colors,
			stringsAsFactors = FALSE
		)
	} else {
		# Check
		if(!is.data.frame(peaks)) stop("'peaks' must be a data.frame")
		
		# Overwrite
		if(!missing(names))    rownames(peaks) <- names
		if(!missing(size.min)) peaks$size.min <- size.min
		if(!missing(size.max)) peaks$size.max <- size.max
		if(!missing(channels)) peaks$channels <- channels
		if(!missing(colors))   peaks$colors   <- colors
		
		# Check
		if(!all(c("size.min", "size.max", "channel", "color") %in% colnames(peaks))) stop("'peaks' must have 'size.min', 'size.max', 'channel' and 'color' columns")
	}
	
	# Points sizes (in bp)
	if(is.null(attr(x, "ladderModel"))) stop("Can't use 'bp' unit without alignment ('ladderModel' attribute)")
	size <- attr(x, "ladderModel")[2] * 1:nrow(x) + attr(x, "ladderModel")[1]
	
	# Identify first peak
	firstPeak <- which.min(peaks$size.min)
	
	for(i in 1:nrow(peaks)) {
		# Scans in range
		w <- size >= peaks[i,"size.min"] & size <= peaks[i,"size.max"]
		if(!any(w)) stop("No data point in range for ", rownames(peaks)[i], " (check alignment)")
		y <- x[ w , peaks[i,"channel"] ]
		
		# Off-scale value in range (any channel)
		isOffScale <- any(attr(x, "offScale") %in% which(w))
		
		# Data collection
		peaks[i,"size"] <- size[w][ which.max(y) ]
		peaks[i,"height"] <- max(y)
		peaks[i,"offScale"] <- isOffScale
	}
	
	# First peak warning
	noiseIndex <- size >= peaks[firstPeak,"size.min"] + noiseRange[1] & size <= peaks[firstPeak,"size.min"] + noiseRange[2]
	noiseValue <- x[ noiseIndex , peaks[firstPeak,"channel"] ]
	firstPeakValue <- peaks[firstPeak,"height"]
	if(firstPeakValue < max(noiseValue) * 1.2) warning("First peak value compromised (value: ", round(firstPeakValue, 1), ", noise: ", round(max(noiseValue), 1), ")")
	
	# Warn about low values
	maxValue <- max(peaks$height)
	if(maxValue < lowThreshold) warning("Low profile (maximal height: ", round(maxValue, 1), ")")
	
	# Warn about off-scale values
	if(any(peaks$offScale)) warning(sum(peaks$offScale), " off-scale measures (", paste(rownames(peaks)[peaks$offScale], collapse=", "), ")")
	
	# Normalize
	peaks$normalized <- peaks$height / mean(peaks$height)
	
	# log2-transformation
	if(isTRUE(logTransform)) {
		peaks[ peaks$normalized < 0L , "normalized" ] <- 0L
		peaks$normalized <- log(1L + peaks$normalized, 2L)
	}
	
	# Store in object
	attr(x, "peaks") <- peaks
	attr(x, "normalized") <- peaks$normalized
	names(attr(x, "normalized")) <- rownames(peaks)
	
	return(x)
}

