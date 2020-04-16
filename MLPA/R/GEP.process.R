# MLPA Gene Expression Profiling TCL-TK worker
# Author : Sylvain Mareschal <maressyl@gmail.com>
GEP.process <- function(
		input,
		design,
		output,
		overwrite = list(),
		gene.cex = 1.3,
		file.line = 3,
		mar = c(5,4,5,1),
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
	
	# GEP.process arguments
	process.args <- design$GEP.process
	for(varName in names(process.args)) {
		assign(varName, process.args[[varName]])
	}
	
	# Switchs
	alignEnabled        <- ! isTRUE(design$align.fsa$disable)
	rescueEnabled       <- ! isTRUE(design$align.fsa$disable) && isTRUE(design$align.fsa$rescue)
	peaksAnnotated      <- ! isTRUE(design$align.fsa$disable) && length(design$PEAKS$ranges) > 0
	plotEnabled         <- ! isTRUE(design$plot.fsa$disable)
	modelEnabled        <- ! isTRUE(design$model$disable)
	classifyEnabled     <- ! isTRUE(design$model$disable) && ! isTRUE(design$classify$disable)
	classifyPlotEnabled <- ! isTRUE(design$model$disable) && ! isTRUE(design$classify$disable) && isTRUE(design$classify$plot)
	
	# Open plotting device
	if(plotEnabled || rescueEnabled) {
		pdf(file=sub("\\.[^\\.]+$", ".pdf", output), title=sub("\\.[^\\.]+$", "", basename(output)), width=24, height=6, pointsize=14)
		par(mar=mar)
	}
	
	# File list
	toProcess <- dir(input, full.names=TRUE, recursive=TRUE, pattern="\\.fsa$")
	if(!is.null(progressBar)) {
		tcltk::tcl(progressBar, "configure", maximum=length(toProcess))
		tcltk::tcl(progressBar, "configure", value=0)
	}
	
	# Model object
	if(modelEnabled) {
		model.args <- design$model
		model.args$disable <- NULL
		model <- do.call("model", model.args)
	}
	
	peakTable <- NULL
	peakMatrix <- NULL
	sampleTable <- NULL
	
	for(f in toProcess) {
		# Log
		message("Processing ", f)
		
		# Import
		read.fsa.args <- design$read.fsa
		read.fsa.args$file <- f
		read.fsa.args$disable <- NULL
		x <- do.call("read.fsa", read.fsa.args)
		
		# Meta data
		meta <- attr(x, "runMetaData")
		message("FSA metadata : ", paste(sprintf("%s=\"%s\"", names(meta), sapply(meta, "[", 1L)), collapse=", "))
		
		# Align
		if(alignEnabled) {
			align.fsa.args <- design$align.fsa
			align.fsa.args$x <- x
			align.fsa.args$title <- sprintf("%s (alignment)", basename(f))
			align.fsa.args$disable <- NULL
			x <- do.call("align.fsa", align.fsa.args)
		}
		
		# Plot
		if(plotEnabled) {
			# Layout
			if(classifyPlotEnabled) { layout(matrix(1:2, nrow=1), widths=c(7,2))
			} else                  { layout(1)
			}
			
			# Plot profile
			plot.fsa.args <- design$plot.fsa
			plot.fsa.args$x <- x
			plot.fsa.args$title <- basename(f)
			plot.fsa.args$disable <- NULL
			do.call("plot.fsa", plot.fsa.args)
			
			if(peaksAnnotated) {
				# Plot peak intervals
				rect(xleft=sapply(design$PEAKS$ranges, "[", 1), xright=sapply(design$PEAKS$ranges, "[", 2), ybottom=-1e6, ytop=1e6, col=design$PEAKS$backgrounds, border=NA)
				text(x=sapply(design$PEAKS$ranges, mean), y=par("usr")[4]+diff(par("usr")[3:4])/50, labels=names(design$PEAKS$ranges), srt=30, adj=c(0, 0), cex=gene.cex, font=2, xpd=NA, col=design$PEAKS$colors)
				box()
			}
		}
		
		if(peaksAnnotated) {
			# Peak detection
			peaks.fsa.args <- design$peaks.fsa
			peaks.fsa.args$x <- x
			peaks.fsa.args$ranges <- design$PEAKS$ranges
			peaks.fsa.args$channels <- design$PEAKS$channels
			peaks.fsa.args$disable <- NULL
			peaks <- do.call("peaks.fsa", peaks.fsa.args)
			
			# Add peaks to raw output
			peakTable <- rbind(peakTable, cbind(file=basename(f), peaks))
			
			# Add peaks to expression matrix
			row <- peaks$normalized
			names(row) <- peaks$gene
			peakMatrix <- rbind(peakMatrix, row)
			rownames(peakMatrix)[ nrow(peakMatrix) ] <- basename(f)
			
			# Plot peak values
			if(plotEnabled) {
				abline(v=peaks$peak.size, col=design$PEAKS$colors)
				box()
			}
			
			if(classifyEnabled) {
				# Classify class
				classify.args <- design$classify
				classify.args$peaks <- peaks
				classify.args$model <- model
				classify.args$disable <- NULL
				pred <- do.call("classify", classify.args)
				
				# Add prediction to output
				predTable <- data.frame(file=basename(f), score=pred$score, p1=pred$p[1], p2=pred$p[2], class=pred$class, row.names=FALSE)
				colnames(predTable)[3:4] <- sprintf("p.%s", names(pred$p))
				sampleTable <- rbind(sampleTable, predTable)
			}
		}
		
		# Progress bar
		if(!is.null(progressBar)) {
			tcltk::tcl(progressBar, "step", 1)
			tcltk::tcl("update")
		}
	}
	
	# Close plotting device
	if(plotEnabled || rescueEnabled) dev.off()
	
	if(alignEnabled) {
		if(peaksAnnotated) {
			# Write normalized matrix (by total RNA quantity)
			write.table(peakMatrix, sub("\\.[^\\.]+$", ".expr.tsv", output), sep="\t", dec=".", row.names=TRUE, col.names=NA, quote=FALSE)
	
			# Write peak table
			write.table(peakTable, sub("\\.[^\\.]+$", ".peaks.tsv", output), sep="\t", dec=".", row.names=FALSE, col.names=TRUE, quote=FALSE)
		
			if(classifyEnabled) {
				# Write sample table
				write.table(sampleTable, sub("\\.[^\\.]+$", ".pred.tsv", output), sep="\t", dec=".", row.names=FALSE, col.names=TRUE, quote=FALSE)
			} else message("Partial output : as classification is disabled, no 'pred.tsv' file have been produced")
		} else message("Partial output : as no peak interval is defined, no 'tsv' file have been produced")
	} else message("Partial output : as alignment is disabled, no 'tsv' file have been produced")
	
	# Done
	message("All done")
}

