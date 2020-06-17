# Generic .fsa collection processing pipeline
# Author : Sylvain Mareschal <maressyl@gmail.com>
generic.process <- function(
		input,
		design,
		output,
		overwrite = list(),
		progressBar = NULL
	) {
	# Version
	message(date())
	message("MLPA ", as.character(packageVersion("MLPA")), "\n", sep="")
	
	# Checks
	if(length(dir(input, pattern="\\.fsa$", recursive=TRUE)) == 0) stop("Select an input directory containing at least one .fsa file to process.")
	if(design == "" || !grepl("\\.(conf)$", design))               stop("Select a design file with a \".conf\" extension.")
	if(!file.exists(design))                                       stop("Select an existing design file.")
	if(output == "")                                               stop("Select an output file name.")
	
	# Design file used
	message("Design file : ", normalizePath(design), " [", tools::md5sum(design), "]", sep="")
	
	# Design file processing
	design <- designFile(design, overwrite=overwrite)
	
	# File list
	toProcess <- dir(input, full.names=TRUE, recursive=TRUE, pattern="\\.fsa$")
	if(!is.null(progressBar)) {
		tcltk::tcl(progressBar, "configure", maximum=length(toProcess))
		tcltk::tcl(progressBar, "configure", value=0)
	}
	
	# Initialize globals
	globals <- list(
		FILE_PATH = NULL,
		FILE_NAME = NULL,
		OBJECT = NULL,
		RANGES = design$PEAKS$ranges,
		CHANNELS = design$PEAKS$channels,
		COLORS = design$PEAKS$colors,
		BACKGROUNDS = design$PEAKS$backgrounds,
		OUTPUT_PATH = output,
		OUTPUT_NAME = basename(output)
	)
	
	for(file in toProcess) {
		
		message("\nProcessing ", file)
		globals$FILE_PATH <- file
		globals$FILE_NAME <- basename(file)
		
		# Call functions listed in design order
		for(i in 1:length(design)) {
			if(!names(design)[i] %in% c("DESIGN", "PEAKS")) {
				# Loop modifiers
				modifier <- attr(design[[i]], "modifier")
				if(identical(modifier, "first") && file != head(toProcess, 1)) next
				if(identical(modifier, "last")  && file != tail(toProcess, 1)) next
				
				message("- ", names(design)[i])
				
				# Collect arguments from the design
				args <- design[[i]]
				
				# Replace globals
				if(length(args) > 0L) for(a in 1:length(args)) {
					for(g in names(globals)) {
						if(length(args[[a]]) == 1L && !is.na(args[[a]])) {
							# String replacement
							regex <- sprintf("\\$%s", g)
							if(grepl(regex, args[[a]])) {
								args[[a]] <- gsub(regex, globals[[g]], args[[a]])
							}
							
							# Variable replacement
							pattern <- sprintf("@%s", g)
							if(args[[a]] == pattern) {
								args[[a]] <- globals[[g]]
							}
						}
					}
				}
				
				# Call function
				out <- do.call(names(design)[i], args)
				
				# OBJECT global
				if(is(out, "fsa")) globals$OBJECT <- out
			}
		}
		
		# Progress bar
		if(!is.null(progressBar)) {
			tcltk::tcl(progressBar, "step", 1)
			tcltk::tcl("update")
		}
	}
	
	# Done
	message("All done")
}

