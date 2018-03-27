# Imports a .fsa file from Applied Biosystems, using the provided converter
# Author : Sylvain Mareschal <maressyl@gmail.com>
read.fsa <- function(
		file,
		applyLowess = TRUE,
		processed = FALSE,
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
	
	# Collect metadata
	meta <- list()
	
	# Metadata (various)
	if("User.1" %in% names(ab1$Data))               meta$user <- ab1$Data$User.1
	if("MCHN.1" %in% names(ab1$Data))               meta$machine <- ab1$Data$MCHN.1
	if(all(c("RMdN", "RMdV") %in% names(ab1$Data))) meta$runModule <- paste(ab1$Data$RMdN, ab1$Data$RMdV, sep=" ")
	if(all(c("RPrN", "RPrV") %in% names(ab1$Data))) meta$runProtocol <- paste(ab1$Data$RPrN, ab1$Data$RPrV, sep=" ")
	
	# Metadata (injection time)
	regex <- "^.+<Token>DC_Injection_Time</Token><Value>(.+?)</Value>.+$"
	if("RMdX.1" %in% names(ab1$Data) && grepl(regex, ab1$Data$RMdX.1)) meta$injectionTime <- sub(regex, "\\1", ab1$Data$RMdX.1)
	
	# Metadata (date)
	if(all(c("RUND.1", "RUNT.1") %in% names(fsa$Data)) && is.list(fsa$Data$RUND.1) && is.list(fsa$Data$RUNT.1)) {
		meta$startDate <- sprintf(
			"%04i-%02i-%02i %02i:%02i:%02i",
			fsa$Data$RUND.1$year,
			fsa$Data$RUND.1$month,
			fsa$Data$RUND.1$day,
			fsa$Data$RUNT.1$hour,
			fsa$Data$RUNT.1$minute,
			fsa$Data$RUNT.1$second
		)
	}
	
	# Store metadata
	attr(x, "runMetaData") <- meta
	
	# Off scale values (if any)
	attr(x, "offScale") <- as.integer(fsa$Data$OfSc.1) + 1L
	
	# S3 class
	class(x) <- "fsa"
	
	return(x)
}
