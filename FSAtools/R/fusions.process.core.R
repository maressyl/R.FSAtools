fusions.process.core <- function(
		sampleName,
		samples,
		design,
		input,
		output,
		channelOrder = NULL,
		top = 5,
		dev.fun = cairo_pdf,
		dev.args = list(width=15, height=10),
		dev.close = TRUE
	) {
	# Check the file pair
	forwardFile <- sprintf("%s/%s", input, samples[ samples$ID == sampleName & samples$way == "forward" , "file" ])
	reverseFile <- sprintf("%s/%s", input, samples[ samples$ID == sampleName & samples$way == "reverse" , "file" ])
	if(length(forwardFile) != 1L) stop("No or multiple forward files found for sample \"", sampleName, "\"")
	if(length(reverseFile) != 1L) stop("No or multiple reverse files found for sample \"", sampleName, "\"")
	
	# Parse the file pair
	message("Parsing files...")
	forward <- read.sanger(forwardFile, channelOrder=channelOrder)
	reverse <- read.sanger(reverseFile, channelOrder=channelOrder)
	
	# 5 best alignments
	message("Identifying fusions...")
	idn <- identify.fusions(forward=forward, reverse=reverse, design=design, top=top)
	
	# Annotate output table
	tab <- idn$top
	tab$sample <- rep(sampleName, nrow(tab))
	tab$forwardFile <- rep(forwardFile, nrow(tab))
	tab$reverseFile <- rep(reverseFile, nrow(tab))
	
	# Plot identified fusions
	message("Plotting the alignment...")
	dev.args$file <- sprintf("%s/%s.pdf", output, gsub("/", "_", sampleName))
	do.call(dev.fun, args=dev.args)
	plot.fusions(sampleName=sampleName, forward=forward, reverse=reverse, forwardFile=forwardFile, reverseFile=reverseFile, identified=idn, design=design)
	if(isTRUE(dev.close)) dev.off()
	
	message("Done")
	
	return(tab)
}
