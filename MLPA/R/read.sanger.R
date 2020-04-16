# Imports a .ab1 / .fsa file from Applied Biosystems, using the provided converter
# Author : Sylvain Mareschal <maressyl@gmail.com>
read.sanger <- function(file, channelOrder=NULL, guess.threshold=0.3) {
	# Parse the ABIF file
	obj <- read.fsa(file, processed=TRUE, applyLowess=FALSE, meta.extra=c(sequence="PBAS", quality="PCON", peaks="PLOC"))

	# Called sequence
	seq <- strsplit(attr(obj, "runMetaData")$sequence["PBAS.1"], split="")[[1]]
	filter <- seq %in% c("A", "C", "G", "T")
	seq <- seq[filter]
	
	# Quality
	phred <- as.integer(charToRaw(attr(obj, "runMetaData")$quality["PCON.1"]))
	phred <- phred[filter]
	
	# Peak positions
	peaks <- attr(obj, "runMetaData")$peaks
	peaks <- as.integer(peaks[ grep("^PLOC\\.1", names(peaks)) ])
	
	# Enforce size consistency
	size <- min(length(seq), length(phred), length(peaks))
	if(length(seq) > size || length(phred) > size || length(peaks) > size) {
		warning("Base call (", length(seq), "), Phred score (", length(phred), ") and peak (", length(peaks), ") lengths differ, trimming assuming common start")
		seq <- seq[1:size]
		phred <- phred[1:size]
		peaks <- peaks[1:size]
	}
	
	# Store
	attr(obj, "seq") <- seq
	attr(obj, "phred") <- phred
	attr(obj, "peaks") <- peaks
	
	# Mask off-scale values (different X axis)
	attr(obj, "offScale") <- NULL
	
	if(is.null(channelOrder)) {
		# Guess channels from called sequence
		column <- NULL
		i <- tapply(X=1:length(seq), INDEX=seq, FUN=c)
		for(nt in c("A", "C", "G", "T")) {
			j <- i[[nt]]
			j <- j[ j <= length(peaks) & j <= length(phred) ]
			column[nt] <- weighted.mean(apply(obj[ peaks[j] ,], 1, which.max), phred[j])
		}
		
		# Reshape
		column <- sort(column)
		if(any(abs(column - 1:4) > guess.threshold)) stop("Channel guessing probably failed (", max(abs(column - 1:4)), ")")
		channels <- names(column)
		chanColors <- c("G"="black", "A"="darkgreen", "T"="darkred", "C"="darkblue")[ names(column) ]
		
		# Store
		colnames(obj) <- channels
		attr(obj, "colors") <- chanColors
	} else if(setequal(channelOrder, c("A","C","G","T"))) {
		# Custom ordering
		colnames(obj) <- channelOrder
		attr(obj, "colors") <- c("G"="black", "A"="darkgreen", "T"="darkred", "C"="darkblue")[ channelOrder ]
	} else stop("'channelOrder' must either be an ordered vector of 'A', 'C', 'G' and 'T' or NULL")
	
	# S3 class
	class(obj) <- "sanger"
	
	return(obj)
}

