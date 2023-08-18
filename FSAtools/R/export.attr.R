# Append the content of an fsa object attribute to a CSV file
export.attr <- function(
		x,
		attr,
		file,
		meta = character(0),
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
	
	# Transform vectors into 1-row data.frame
	row.names <- TRUE
	if(is.atomic(val)) {
		if(is.null(names(val))) names(val) <- sprintf("X%i", 1:length(val))
		val <- as.data.frame(as.list(val))
		row.names <- FALSE
	}
	
	if(is.data.frame(val)) {
		# Append metadata columns
		for(m in meta) {
			v <- attr(x, "metaData")[[m]]
			if(length(v) == 0L) { val[[m]] <- NA
			} else              { val[[m]] <- paste(v, collapse=", ")
			}
		}
		
		# Append with an extra 'sample' column
		if(file.exists(file)) { col.names <- FALSE
		} else if(row.names)  { col.names <- NA
		} else                { col.names <- TRUE
		}		
		write.table(val, file=file, sep=sep, dec=dec, quote=quote, append=TRUE, row.names=row.names, col.names=col.names)
	} else stop("Attribute to export must be a data.frame or atomic vector")
	
	invisible(TRUE)
}
