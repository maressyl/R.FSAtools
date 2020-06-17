# Parses the design file
# Author : Sylvain Mareschal <maressyl@gmail.com>
designFile <- function(fileName, overwrite=list()) {
	# Design sections
	rawDesign <- readLines(fileName)
	sectionStarts <- grep("^\\s*\\[([A-Za-z\\._:]+)\\]\\s*$", rawDesign)
	sectionNames <- sub("^\\s*\\[([A-Za-z\\._:]+)\\]\\s*$", "\\1", rawDesign[sectionStarts])
	boundaries <- c(sectionStarts, length(rawDesign)+1L)
	
	# Returned list
	conf <- list()
	
	
	
	### name:value sections (DESIGN and functions)
	
	sections <- setdiff(sectionNames, "PEAKS")
	for(section in sections) {
		# Preserve design file order
		n <- length(conf) + 1L
		
		# Section boundaries
		skip <- sectionStarts[ sectionNames == section ]
		nrows <- boundaries[ boundaries > skip ][1] - skip - 1L
		
		# Extract from whole design
		content <- rawDesign[skip+(1:nrows)]
		
		# Discard comments and empty lines
		content <- grep("^[^#\t]+\t.+\\s*$", content, value=TRUE)
		if(length(content) > 0) {
			# Split name and values
			content <- strsplit(content, split="\t")
		
			# Reformate as a named list, separating vectors
			conf[[n]] <- lapply(content, "[", -1)
			names(conf[[n]]) <- sapply(content, "[", 1)
			
			# Guess types
			for(i in 1:length(conf[[n]])) conf[[n]][[i]] <- type.convert(conf[[n]][[i]], as.is=TRUE)
		} else {
			# Allocate an empty list anyway
			conf[[n]] <- list()
		}
		
		# Modifiers
		if(grepl("^(.+):(first|last)$", section)) {
			modifier <- sub("^(.+):(first|last)$", "\\2", section)
			section <- sub("^(.+):(first|last)$", "\\1", section)
			attr(conf[[n]], "modifier") <- modifier
		}
		
		# Section name
		names(conf)[n] <- section
		
		if(section != "DESIGN") {
			# Section name must be an existing function
			if(!exists(section, mode="function")) stop("Design section \"", section, "\" refers to an unknown R function")
			
			# Function arguments
			arguments <- formals(section)
			
			# Defined in design file but not used by this function
			notInFun <- setdiff(names(conf[[n]]), names(arguments))
			if(length(notInFun) > 0L) message("Design setting(s) ignored for function ", section, " : ", paste(notInFun, collapse=", "))
			
			# Used by function but not defined in design
			notInDesign <- setdiff(names(arguments), names(conf[[n]]))
			if(length(notInDesign) > 0L) message("Design setting(s) set to default for function ", section, " : ", paste(notInDesign, collapse=", "))
		}
	}
	
	
	
	### PEAKS (FIXME)
	
	# Section boundaries
	skip <- sectionStarts[ sectionNames == "PEAKS" ]
	nrows <- boundaries[ boundaries > skip ][1] - skip - 1L
	
	# Check if content is empty
	content <- rawDesign[skip+(1:nrows)]
	content <- grep("^[^#\t]+\t.+\\s*$", content, value=TRUE)
	if(length(content) > 0) {
		# Parse
		peaks <- read.table(fileName, skip=skip, nrows=nrows, sep="\t", dec=".", row.names="name", comment.char="#", header=TRUE, stringsAsFactors=FALSE)
		conf$PEAKS <- list()
	
		# Genes
		conf$PEAKS$ranges <- as.list(as.data.frame(t(peaks[,c("size.min", "size.max")])))
		conf$PEAKS$channels <- peaks$channel
		conf$PEAKS$colors <- peaks$color
	
		# Add transparency to colors, turn "" into NA
		conf$PEAKS$colors[ conf$PEAKS$colors == "" ] <- NA
		trueColors <- !is.na(conf$PEAKS$colors)
		conf$PEAKS$backgrounds <- rep(as.character(NA), length(conf$PEAKS$colors))
		conf$PEAKS$backgrounds[ trueColors ] <- sprintf("#%s", apply(as.character(as.hexmode(rbind(col2rgb(conf$PEAKS$colors[ trueColors ]), 48L))), 2, paste, collapse=""))
	} else {
		# No peak
		conf$PEAKS <- list()
		conf$PEAKS$ranges <- integer(0)
		conf$PEAKS$channels <- character(0)
		conf$PEAKS$colors <- character(0)
		conf$PEAKS$backgrounds <- character(0)
	}
	
	
	
	### CHECKS
	
	# Overwrite some options from CLI
	if(is.list(overwrite) && length(overwrite) > 0L) {
		for(section in names(overwrite)) {
			for(argument in names(overwrite[[section]])) {
				conf[[section]][[argument]] <- overwrite[[section]][[argument]]
			}
		}
	}
	
	# Version
	softVersion <- as.character(packageVersion("MLPA"))
	if(softVersion != conf$DESIGN$MLPA) warning(sprintf("Design (%s) and software (%s) version differ", conf$DESIGN$MLPA, softVersion))
	
	# Can't disable read.fsa()
	if(isTRUE(conf$read.fsa$disable)) stop("Can't disable read.fsa()")
	
	
	
	return(conf)
}

