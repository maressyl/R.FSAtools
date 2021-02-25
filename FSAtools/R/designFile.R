# Parses the design file
designFile <- function(fileName) {
	# Design sections
	rawDesign <- readLines(fileName)
	regex <- "^\\s*\\[(.+)\\]\\s*$"
	sectionStarts <- grep(regex, rawDesign)
	sectionNames <- sub(regex, "\\1", rawDesign[sectionStarts])
	boundaries <- c(sectionStarts, length(rawDesign)+1L)
	
	# Returned list
	conf <- list(GLOBALS=NULL)
	
	for(s in 1:length(sectionNames)) {
		# Preserve design file order
		n <- length(conf) + 1L
		
		# Modifiers
		if(grepl(":", sectionNames[s])) {
			# Extract elements
			section <- sub(":.+$", "", sectionNames[s])
			modifiers <- sub("^.+:", "", sectionNames[s])
			modifiers <- strsplit(modifiers, split=",", fixed=TRUE)[[1]]
			
			# Validate modifiers
			unknown <- setdiff(modifiers, c("first", "last", "table", "nowarn"))
			if(length(unknown) > 0L) stop("Unknown modifier(s) : ", paste(unknown, collapse=", "))
		} else {
			# No modifier
			modifiers <- NULL
			section <- sectionNames[s]
		}
		
		# Type
		if(grepl("^[A-Z_]+$", section) && !exists(section, mode="function")) { type <- "global"
		} else if(exists(section, mode="function"))                          { type <- "function"
		} else                                                               { stop("Design section \"", section, "\" is neither a R function or a valid global declaration")
		}
		
		# Subtype
		if(type == "global") {
			if("table" %in% modifiers) { type <- "global.table"
			} else                     { type <- "global.list"
			}
		}
		
		# Section boundaries
		skip <- sectionStarts[s]
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
			attr(conf[[n]], "modifiers") <- modifiers
			
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
	
	# Version
	if("DESIGN" %in% names(conf$GLOBALS)) {
		if("FSAtools" %in% names(conf$GLOBALS$DESIGN)) {
			softVersion <- as.character(packageVersion("FSAtools"))
			if(softVersion != conf$GLOBALS$DESIGN$FSAtools) {
				warning(sprintf("Design (%s) and software (%s) version differ", conf$GLOBALS$DESIGN$FSAtools, softVersion))
			}
		} else if("MLPA" %in% names(conf$GLOBALS$DESIGN))   {
			stop("Design file pretends to be from the previous (uncompatible) version of FSAtools 'MLPA'")
		} else {
			warning("DESIGN section of the design file contains no FSAtools element with the expected package version")
		}
	} else {
		warning("Design file contains no DESIGN section")
	}
	
	return(conf)
}

