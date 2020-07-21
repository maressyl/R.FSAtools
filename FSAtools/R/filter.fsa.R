# Applies filter() to a channel of a FSA object
filter.fsa <- function(x, channel, ..., from=NA, to=NA, units="bp") {
	if(!is.na(from) && !is.na(to)) {
		# X coordinate
		if(units == "index") {
			if(nrow(x) > 0) { xcor <- 1:nrow(x)
			} else          { xcor <- integer(0)
			}
		} else if(units == "bp") {
			if(is.null(attr(x, "ladderModel"))) stop("Can't crop an unaligned object in 'bp' units")
			if(nrow(x) > 0) { xcor <- attr(x, "ladderModel")[2] * 1:nrow(x) + attr(x, "ladderModel")[1]
			} else          { xcor <- integer(0)
			}
		}
		
		# Values to use
		i <- xcor >= from & xcor <= to
		
		# Mask out-of-bound values
		x[!i,channel] <- NA
	} else {
		# Use all values
		i <- TRUE
	}
	
	# Apply filter
	x[i,channel] <- filter(x[i,channel], ...)
	
	return(x)
}
