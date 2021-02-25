# Generic .fsa collection processing with log file handler
generic.log <- function(..., logFile) {
	# Initialize log file
	cat("", file=logFile)
	
	# Follow-up
	warnCount <- 0L
	error <- NULL
	
	# Run under handler
	handle(
		expr = {
			generic.process(...)
		},
		messageHandler = function(m) {
			mes <- conditionMessage(m)
			if(grepl("^Processing", mes) || grepl("^Partial output", mes) || grepl("^All done", mes)) {
			         cat(sprintf("\n%-10s%s", "", mes), file=logFile, append=TRUE)
			} else { cat(sprintf("%-10s%s", "", mes), file=logFile, append=TRUE)
			}
		},
		warningHandler = function(w) {
			cat(sprintf("%-10s%s\n", "[WARNING]", conditionMessage(w)), file=logFile, append=TRUE)
			warnCount <<- warnCount + 1L
		},
		errorHandler = function(e) {
			cat(sprintf("%-10s%s\n", "[ERROR]", conditionMessage(e)), file=logFile, append=TRUE)
			error <<- e
			if(dev.cur() > 1L) dev.off()
		}
	)
	
	# Return
	if(!is.null(error)) {
		# Failed
		out <- error
	} else if(warnCount > 0L) {
		# Success with warnings
		out <- warnCount
	} else {
		# Success without warnings
		out <- TRUE
	}
	
	return(out)
}

