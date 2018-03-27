# Imports a .fsa file from Applied Biosystems, using the provided converter
# Author : Sylvain Mareschal <maressyl@gmail.com>
read.fsa <- function(
		file,
		applyLowess = TRUE,
		processed = FALSE,
		meta.extra = NULL,
		...
	) {
	# Parse ABIF
	fsa <- read.abif(file, ...)
	
	# Scan dimensions
	channelCount <- fsa$Data$`Dye#.1`
	if("SCAN.1" %in% names(fsa$Data)) { scanCount <- fsa$Data$SCAN.1
	} else                            { scanCount <- fsa$Data$Scan.1 ### Legacy
	}
	if(scanCount == 0)    warning("No scan stored in this file")
	if(channelCount == 0) warning("No channel stored in this file")
	
	# Extract data
	dyeNames <- character(0)
	dyeWavelengths <- integer(0)
	dyeColors <- character(0)
	
	# DATA indexes
	if(is.na(processed))  { processed <- any(sprintf("DATA.%i", c(9:12, 205:299)) %in% names(fsa$Data))
	}
	if(isTRUE(processed)) { dataTracks <- c(9:12, 205:299)
	} else                { dataTracks <- c(1:4, 105:199)
	}
	
	channelValues <- list()
	for(i in 1:channelCount) {
		# Get raw DATA
		values <- fsa$Data[[ sprintf("DATA.%i", dataTracks[i]) ]]
		
		# Apply lowess to reduce time bias
		if(isTRUE(applyLowess) && scanCount > 0) values <- values - lowess(x=1:length(values), y=values)$y
		
		# Update value matrix
		channelValues[[i]] <- values
		
		# Get dye name
		if(sprintf("DyeN.%i", i) %in% names(fsa$Data))                      dyeNames[i] <- fsa$Data[[ sprintf("DyeN.%i", i) ]]
		if(dyeNames[i] == "" && sprintf("DyeS.%i", i) %in% names(fsa$Data)) dyeNames[i] <- fsa$Data[[ sprintf("DyeS.%i", i) ]]
		
		# Get dye wavelength
		if(sprintf("DyeW.%i", i) %in% names(fsa$Data)) {
			dyeWavelengths[i] <- fsa$Data[[ sprintf("DyeW.%i", i) ]]
			dyeColors[i] <- wav2RGB(dyeWavelengths[i])
		} else {
			dyeWavelengths[i] <- NA
			dyeColors[i] <- NA
		}
	}
	
	# Channel size consistency
	channelLengths <- sapply(channelValues, length)
	if(isTRUE(processed)) scanCount <- channelLengths[1]
	if(any(channelLengths != scanCount)) stop("Data length inconsistency: ", paste(channelLengths, collapse=", "))
	
	# Store channels
	x <- matrix(as.double(NA), ncol=channelCount, nrow=scanCount)
	for(i in 1:channelCount) x[,i] <- channelValues[[i]]
	
	# Automatic colors
	if(all(is.na(dyeColors))) {
		dyeColors <- rainbow(channelCount)
		warning("No dye wavelength found, using arbitrary colors")
	}
	
	# Use dye names
	names(dyeWavelengths) <- dyeNames
	names(dyeColors) <- dyeNames
	colnames(x) <- dyeNames
	
	# Attributes
	attr(x, "lowess") <- isTRUE(applyLowess)
	attr(x, "wavelengths") <- dyeWavelengths
	attr(x, "colors") <- dyeColors
	
	# Metadata to collect
	collect <- c(
		user = "User",
		machine = "MCHN",
		runModule.name = "RMdN",
		runModule.version = "RMdV",
		runProtocole.name = "RPrN",
		runProtocole.version = "RPrV",
		runDate = "RUND",
		runTime = "RUNT",
		meta.extra
	)
	
	# Start collection
	meta <- list()
	for(metaName in names(collect)) {
		values <- fsa$Data[ grep(sprintf("^%s\\.[0-9]+$", collect[ metaName ]), names(fsa$Data)) ]
		if(length(values) > 0L) {
			if(all(sapply(values, is.atomic))) {
				meta[[ metaName ]] <- unlist(values)
				if(length(values) == 1L) names(meta[[ metaName ]]) <- NULL
			} else {
				meta[[ metaName ]] <- as.matrix(as.data.frame(lapply(values, unlist)))
			}
		}
	}
	
	# Reshape dates
	if("runDate" %in% names(meta) && "runTime" %in% names(meta)) {
		dates <- sprintf("%04i-%02i-%02i %02i:%02i:%02i", meta$runDate["year",], meta$runDate["month",], meta$runDate["day",], meta$runTime["hour",], meta$runTime["minute",], meta$runTime["second",])
		names(dates) <- c("runStart", "runStop", "collectionStart", "collectionStop")[ 1:length(dates) ]
		meta$runDate <- NULL
		meta$runTime <- NULL
		meta <- c(meta, dates)
	}
	
	# Injection time (not portable)
	regex <- "^.+<Token>DC_Injection_Time</Token><Value>(.+?)</Value>.+$"
	if("RMdX.1" %in% names(fsa$Data) && grepl(regex, fsa$Data$RMdX.1)) meta$injectionTime <- sub(regex, "\\1", fsa$Data$RMdX.1)
	
	# Store metadata
	attr(x, "runMetaData") <- meta
	
	# Off scale values (if any)
	attr(x, "offScale") <- as.integer(fsa$Data$OfSc.1) + 1L
	
	# S3 class
	class(x) <- "fsa"
	
	return(x)
}
