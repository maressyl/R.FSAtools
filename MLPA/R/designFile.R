# Parses the design file
# Author : Sylvain Mareschal <maressyl@gmail.com>
designFile <- function(fileName, overwrite=list()) {
	# Design sections
	rawDesign <- readLines(fileName)
	sectionStarts <- grep("^\\s*\\[([A-Za-z\\._:]+)\\]\\s*$", rawDesign)
	sectionNames <- sub("^\\s*\\[([A-Za-z\\._:]+)\\]\\s*$", "\\1", rawDesign[sectionStarts])
	boundaries <- c(sectionStarts, length(rawDesign)+1L)
	
	# Returned list
	conf <- list(GLOBALS=NULL)
	
	for(sectionName in sectionNames) {
		# Preserve design file order
		n <- length(conf) + 1L
		
		# Modifiers
		if(grepl("^(.+):(first|last|table)$", sectionName)) {
			modifier <- sub("^(.+):(first|last|table)$", "\\2", sectionName)
			section <- sub("^(.+):(first|last|table)$", "\\1", sectionName)
		} else {
			modifier <- NULL
			section <- sectionName
		}
		
		# Type
		if(grepl("^[A-Z_]+$", section) && !exists(section, mode="function")) { type <- "global"
		} else if(exists(section, mode="function"))                          { type <- "function"
		} else                                                               { stop("Design section \"", section, "\" is neither a R function or a valid global declaration")
		}
		
		# Subtype
		if(type == "global") {
			if(identical(modifier, "table")) { type <- "global.table"
			} else                           { type <- "global.list"
			}
		}
		
		# Section boundaries
		skip <- sectionStarts[ sectionNames == sectionName ]
		nrows <- boundaries[ boundaries > skip ][1] - skip - 1L
		
		# Extract from whole design
		content <- rawDesign[skip+(1:nrows)]
		
		# Discard comments and empty lines
		content <- grep("^[^#\t]+\t.+\\s*$", content, value=TRUE)
		
		if(length(content) > 0L) {
			if(type %in% c("function", "global.list")) {
				# Split name and values
				content <- strsplit(content, split="\t")
			
				# Reformate as a named list, separating vectors
				processed <- lapply(content, "[", -1)
				names(processed) <- sapply(content, "[", 1)
				
				# Guess types
				for(i in 1:length(processed)) processed[[i]] <- type.convert(processed[[i]], as.is=TRUE)
			} else if(type == "global.table") {
				# Parse data.frame
				processed <- read.table(textConnection(content), sep="\t", dec=".", row.names=1, header=TRUE, stringsAsFactors=FALSE)
			} else stop("Unknown type")
		} else {
			# Allocate an empty list anyway
			processed <- list()
		}
		
		if(type == "function") {
			# Store
			conf[[n]] <- processed
			names(conf)[n] <- section
			attr(conf[[n]], "modifier") <- modifier
			
			# Function arguments
			arguments <- formals(section)
			
			# Defined in design file but not used by this function
			notInFun <- setdiff(names(conf[[n]]), names(arguments))
			if(length(notInFun) > 0L && ! "..." %in% names(arguments)) message("Design setting(s) ignored for function ", section, " : ", paste(notInFun, collapse=", "))
			
			# Used by function but not defined in design
			notInDesign <- setdiff(names(arguments), c(names(conf[[n]]), "..."))
			if(length(notInDesign) > 0L) message("Design setting(s) set to default for function ", section, " : ", paste(notInDesign, collapse=", "))
		} else if(type %in% c("global.list", "global.table")) {
			# Store
			conf$GLOBALS[[section]] <- processed
		} else stop("Unknown type")
	}
	
	# Overwrite some options from CLI
	if(is.list(overwrite) && length(overwrite) > 0L) {
		for(section in names(overwrite)) {
			for(argument in names(overwrite[[section]])) {
				conf[[section]][[argument]] <- overwrite[[section]][[argument]]
			}
		}
	}
	
	# Version
	if("DESIGN" %in% names(conf$GLOBALS)) {
		softVersion <- as.character(packageVersion("MLPA"))
		if(softVersion != conf$GLOBALS$DESIGN$MLPA) warning(sprintf("Design (%s) and software (%s) version differ", conf$GLOBALS$DESIGN$MLPA, softVersion))
	} else warning("Design file contains no DESIGN section")
	
	return(conf)
}

