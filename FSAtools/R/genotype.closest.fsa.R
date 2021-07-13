# Call alleles from peak sets
genotype.closest.fsa <- function(
		x
		)
	{
	# Checks
	if(!is(x, "fsa")) stop("'x' must be a 'fsa' object")
	
	# Peak table, with values
	peaks <- attr(x, "peaks")
	if(!is.data.frame(peaks))                        stop("'x' must have been processed with peaks.fsa()")
	if(!all(c("N0","N1","N2") %in% colnames(peaks))) stop("'x' peak table must have N0, N1 and N2 optional columns")
	
	# Loci
	peaks$locus <- sub("^(.+) - (.+)$", "\\2", rownames(peaks))
	peaks$locus <- factor(peaks$locus, levels=unique(peaks$locus))
	peaks$allele <- sub("^(.+) - (.+)$", "\\1", rownames(peaks))
	peaks <- peaks[ order(peaks$locus, peaks$allele) ,]
	
	# Closest key value
	peaks$N <- apply(abs(peaks[,c("N0","N1","N2")] - peaks$normalized), 1, which.min) - 1L
	
	# Call per locus
	peaks$alleles <- mapply(function(...) { paste(rep(...), collapse="") }, peaks$allele, peaks$N)
	loci <- as.data.frame(tapply(X=peaks$alleles, INDEX=peaks$locus, FUN=paste, collapse=""))
	colnames(loci)[1] <- "call"
	
	# Store in object
	attr(x, "genotypes") <- loci
	attr(x, "calls") <- loci$call
	names(attr(x, "calls")) <- rownames(loci)
	
	return(x)
}

