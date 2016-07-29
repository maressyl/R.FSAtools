# Print method for "fsa" S3 class
# Author : Sylvain Mareschal <maressyl@gmail.com>
print.fsa <- function(x, ...) {
	cat("Object of class MLPA::fsa\n\n")
	cat(nrow(x), " values on ", ncol(x), " channels (", paste(colnames(x), collapse=", "), ")\n", sep="")
	
	if(length(attr(x, "offScale")) > 0) cat(length(attr(x, "offScale")), "off-scale values masked\n")
	else                                cat("No off-scale value\n")
	
	if(attr(x, "lowess")) cat("Consistency along time axis enforced by LOWESS\n")
	
	if(is.null(attr(x, "ladderExact"))) cat("No time index to size model available\n")
	else                                cat("Aligned on ", length(attr(x, "ladderExact")), " size markers (", paste(names(attr(x, "ladderExact")), collapse=", "), ")\n", sep="")
}
