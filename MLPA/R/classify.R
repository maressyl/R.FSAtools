# Apply a LPS-like bayesian predictor to peaks produced by peaks.fsa()
# Author : Sylvain Mareschal <maressyl@gmail.com>
classify <- function(
		peaks,
		model,
		plot = TRUE
		)
	{
	# Class checks
	if(!is(model, "fsaModel"))   stop("'model' must be a 'fsaModel' object")
	if(!is(peaks, "data.frame")) stop("'peaks' must be a data.frame")
	
	# Check model genes
	missingGenes <- setdiff(model$geneNames, peaks$gene)
	if(length(missingGenes) > 0)    stop("Genes required by the model not measured : ", paste(missingGenes, collapse=", "))
	if(any(duplicated(peaks$gene))) stop("Duplicated gene names in peak table")
	
	# Cross data
	rownames(peaks) <- peaks$gene
	peaks <- peaks[ model$geneNames ,]
	peaks$T <- model$geneTs
	peaks$M <- model$geneMs
	
	# Compute score
	score <- sum(peaks$T * (peaks$normalized / peaks$M))
	
	# Compute probability
	D1 <- dnorm(score, mean=model$groupMeans[1], sd=model$groupSDs[1])
	D2 <- dnorm(score, mean=model$groupMeans[2], sd=model$groupSDs[2])
	p <- D1 / (D1 + D2)
	p <- c(p, 1 - p[1])
	names(p) <- model$groupNames
	
	# Class
	if(p[1] > model$threshold)        { class <- model$groupNames[1]
	} else if(p[2] > model$threshold) { class <- model$groupNames[2]
	} else                            { class <- as.character(NA)
	}
	
	if(isTRUE(plot)) {
		# Plot values
		xlim <- c(min(score, model$groupMeans - 4*model$groupSDs), max(score, model$groupMeans + 4*model$groupSDs))
		xval <- seq(from=xlim[1], to=xlim[2], length.out=1000)
		yval1 <- dnorm(xval, mean=model$groupMeans[1], sd=model$groupSDs[1])
		yval2 <- dnorm(xval, mean=model$groupMeans[2], sd=model$groupSDs[2])
		ylim <- range(c(yval1, yval2), na.rm=TRUE)
		ylim[2] <- ylim[2] * 1.1
		
		# Plot gaussians
		savePar <- par(mar=c(5,1,5,1))
		on.exit(par(savePar))
		plot(x=xval, y=yval1, type="l", col=model$groupColors[1], xlim=xlim, ylim=ylim, xlab="Score", ylab="", yaxt="n", lwd=3)
		lines(x=xval, y=yval2, col=model$groupColors[2], lwd=3)
		
		# Add current prediction
		abline(v=score)
		mtext(side=3, at=score, text=signif(score, 3))
		
		# Add probabilities
		if(diff(model$groupMeans) > 0) { lowest <- 1; highest <- 2;	
		} else                         { lowest <- 2; highest <- 1;	
		}
		mtext(side=3, line=2.75, adj=1, at=par("usr")[2], font=2, cex=1.5, col=ifelse(p[lowest]  >= model$threshold, model$groupColors[lowest], "darkgrey"),  text=sprintf("p(%s) = %g", model$groupNames[lowest],  signif(p[lowest], 3)))
		mtext(side=3, line=1.25, adj=1, at=par("usr")[2], font=2, cex=1.5, col=ifelse(p[highest] >= model$threshold, model$groupColors[highest], "darkgrey"), text=sprintf("p(%s) = %g", model$groupNames[highest], signif(p[highest], 3)))
		
		# Estimate unclassified regions
		P1 <- yval1 / (yval1 + yval2)
		P2 <- yval2 / (yval1 + yval2)
		grayPoints <- rle(P1 < model$threshold & P2 < model$threshold)
		grayBreaks <- c(1, cumsum(grayPoints$lengths))
		
		# Plot
		for(i in which(grayPoints$values)) {
			# X coordinates
			left <- xval[ grayBreaks[i] ]
			right <- xval[ grayBreaks[i+1] ]
			
			# Plot
			rect(xleft=left, xright=right, ybottom=-100, ytop=100, col="#00000030", border=NA)
		}
	}
	
	return(
		list(
			score = score,
			p = p,
			class = class
		)
	)
}

