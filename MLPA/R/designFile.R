# Parses the design file
# Author : Sylvain Mareschal <maressyl@gmail.com>
designFile <- function(fileName, overwrite=list()) {
	# Design sections
	rawDesign <- readLines(fileName)
	sectionStarts <- grep("^\\s*\\[([A-Za-z\\._]+)\\]\\s*$", rawDesign)
	sectionNames <- sub("^\\s*\\[([A-Za-z\\._]+)\\]\\s*$", "\\1", rawDesign[sectionStarts])
	boundaries <- c(sectionStarts, length(rawDesign)+1L)
	
	# Returned list
	conf <- list()
	
	
	
	### name:value sections (DESIGN and functions)
	
	sections <- setdiff(sectionNames, "PEAKS")
	for(section in sections) {
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
			conf[[section]] <- lapply(content, "[", -1)
			names(conf[[section]]) <- sapply(content, "[", 1)
			
			# Guess types
			for(i in 1:length(conf[[section]])) conf[[section]][[i]] <- type.convert(conf[[section]][[i]], as.is=TRUE)
		}
	}
	
	# Functions handled
	for(fun in c("GEP.process", "read.fsa", "align.fsa", "peaks.fsa", "plot.fsa", "model", "classify")) {
		# Create missing section
		if(! fun %in% names(conf)) conf[[fun]] <- list()
		
		# Function arguments
		arguments <- formals(fun)
		arguments$disable <- FALSE
		
		# Defined in design file but not used by this function
		notInFun <- setdiff(names(conf[[fun]]), names(arguments))
		notInFun <- setdiff(notInFun, "disable")
		if(length(notInFun) > 0) message("Design setting(s) ignored for function ", fun, " : ", paste(notInFun, collapse=", "))
		
		# Used by function but not defined in design
		notInDesign <- setdiff(names(arguments), names(conf[[fun]]))
		setToDefault <- character(0)
		for(a in notInDesign) {
			# "disable" pseudo-argument
			if(a == "disable") conf[[fun]][[ a ]] <- FALSE
			
			# Do not eval missing arguments
			if(!is.name(arguments[[a]]) || as.character(arguments[[a]]) != "") {
				conf[[fun]][[ a ]] <- eval(arguments[[ a ]])
				setToDefault <- c(setToDefault, a)
			}
		}
		if(length(setToDefault) > 0) message("Design setting set to default for function ", fun, " : ", paste(setToDefault, collapse=", "))
	}
	
	
	
	### PEAKS
	
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

