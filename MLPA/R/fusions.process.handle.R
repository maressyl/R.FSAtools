fusions.process.handle <- function(sampleName, ...) {
	out <- NULL
	dialogs <- NULL
	
	# To log an event
	newRow <- function(condition, type, sampleName) {
		data.frame(
			time = Sys.time(),
			PID = Sys.getpid(),
			sample = sampleName,
			type = type,
			text = sub("\n+$", "", conditionMessage(condition)),
			lineBreaks = nchar(sub("^.*[^\n](\n*)$", "\\1", conditionMessage(condition))),
			stringsAsFactors = FALSE
		)
	}
	
	# Process under monitoring
	handle(
		{ out <- fusions.process.core(sampleName, ...) },
		messageHandler = function(m) {
			dialogs <<- rbind(dialogs, newRow(m, "message", sampleName))
		},
		warningHandler = function(w) {
			dialogs <<- rbind(dialogs, newRow(w, "warning", sampleName))
		},
		errorHandler = function(e) {
			dialogs <<- rbind(dialogs, newRow(e, "error", sampleName))
		}
	)
	
	return(list(dialogs=dialogs, output=out))
}

