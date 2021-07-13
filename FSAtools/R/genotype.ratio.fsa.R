# Call alleles from peak sets
genotype.ratio.fsa <- function(
		x,
		homo = 0.85,
		hetero = c(0.3, 0.7)
		)
	{
	# Checks
	if(!is(x, "fsa")) stop("'x' must be a 'fsa' object")
	
	# Peak table, with values
	peaks <- attr(x, "peaks")
	if(!is.data.frame(peaks)) stop("'x' must have been processed with peaks.fsa()")
	
	# Loci
	peaks$locus <- sub("^(.+) - (.+)$", "\\2", rownames(peaks))
	peaks$locus <- factor(peaks$locus, levels=unique(peaks$locus))
	peaks$allele <- sub("^(.+) - (.+)$", "\\1", rownames(peaks))
	peaks <- peaks[ order(peaks$locus, peaks$allele) ,]
	
	# Locus storage
	loci <- data.frame(
		row.names = levels(peaks$locus),
		stringsAsFactors = FALSE
	)
	loci$call <- ""
	
	# Compute allelic ratio
	for(locus in rownames(loci)) {
		i <- as.character(peaks$locus) == locus
		peaks[i,"ratio"] <- (peaks[i,"height"] - min(peaks$height)) / sum(peaks[i,"height"] - min(peaks$height))
	}
	
	# Per allele calls
	peaks[ peaks$ratio > homo , "call" ] <- "homozygous"
	peaks[ peaks$ratio < hetero[2] & peaks$ratio > hetero[1] , "call" ] <- "heterozygous"
	peaks[ peaks$ratio < 1L - homo , "call" ] <- "absent"
	
	# Homozygous
	indexes <- which(!is.na(peaks$call) & peaks$call == "homozygous")
	for(i in indexes) loci[ as.character(peaks[i,"locus"]) , "call" ] <- paste(loci[ as.character(peaks[i,"locus"]) , "call" ], peaks[i,"allele"], peaks[i,"allele"], sep="")
	
	# Heterozygous
	indexes <- which(!is.na(peaks$call) & peaks$call == "heterozygous")
	for(i in indexes) loci[ as.character(peaks[i,"locus"]) , "call" ] <- paste(loci[ as.character(peaks[i,"locus"]) , "call" ], peaks[i,"allele"], sep="")
	
	# Missing one
	i <- nchar(loci$call) == 1L
	loci[i,"call"] <- paste(loci[i,"call"], "?", sep="")
	
	# Missing two
	i <- nchar(loci$call) == 0L
	loci[i,"call"] <- "??"
	
	# Ratios
	ratios <- with(
		peaks[ order(peaks$ratio, decreasing=TRUE) ,],
		tapply(
			X = sprintf("%s (%.0f%%)", allele, ratio*100),
			INDEX = locus,
			FUN = paste, collapse = ", "
		)
	)
	loci[ names(ratios) , "alleles" ] <- as.character(ratios)
	
	# Store in object
	attr(x, "genotypes") <- loci
	attr(x, "calls") <- loci$call
	names(attr(x, "calls")) <- rownames(loci)
	
	return(x)
}

