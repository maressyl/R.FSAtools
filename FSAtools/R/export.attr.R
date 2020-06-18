# Append the content of an fsa object attribute to a CSV file
export.attr <- function(
		x,
		attr,
		file,
		sample,
		sep = "\t",
		dec = ".",
		quote = TRUE
		)
	{
	# Checks
	if(!is(x, "fsa")) stop("'x' must be a 'fsa' object")
	
	# Attribute to export
	val <- attr(x, attr)
	if(is.null(val)) stop("'x' has no \"", attr, "\" attribute to export")
	
	if(is.data.frame(val)) {
		# Append with an extra 'sample' column
		val$sample <- sample
		write.table(val, file=file, sep=sep, dec=dec, quote=quote, append=TRUE, row.names=TRUE, col.names=ifelse(file.exists(file), FALSE, NA))
	} else if(is.atomic(val)) {
		# Get current matrix
		if(file.exists(file)) { mtx <- as.matrix(read.table(file, sep="\t", header=TRUE, row.names=1, stringsAsFactors=FALSE))
		} else                { mtx <- NULL
		}
		
		# Append a new column
		mtx <- rbind(mtx, val)
		rownames(mtx)[ nrow(mtx) ] <- sample
		
		# Export
		write.table(mtx, file=file, sep=sep, dec=dec, quote=quote, row.names=TRUE, col.names=NA)
	} else stop("Attribute to export must be a data.frame or atomic vector")
}

