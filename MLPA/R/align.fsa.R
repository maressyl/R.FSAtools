# Size estimation based on few markers searched in predefined ranges
# Author : Sylvain Mareschal <maressyl@gmail.com>
align.fsa <- function(
		x,
		channel = "ROX",
		fullLadder = c(50, 60, 90, 100, 120, 150, 160, 180, 190, 200, 220, 240, 260, 280, 290, 300, 320, 340, 360, 380, 400),
		useLadder = c(50, 60, 90, 100, 120),
		outThreshold = 0.15,
		noiseLevel = 10,
		surePeaks = 5,
		trim = c("forward", "backward", "none"),
		maskOffScale = FALSE,
		rMin = 0.999,
		rescue = FALSE,
		ylim = NULL,
		...
	) {
	# Args
	trim <- match.arg(trim)
	if(is.null(useLadder)) useLadder <- fullLadder
	
	# Unrecoverable error
	if(!channel %in% colnames(x)) stop("channel \"", channel, "\" not found (available: ", paste(colnames(x), collapse=", "), ")")
	
	# Postpone errors
	rSquared <- NULL
	status <- try(silent=TRUE, expr={
		# Smoothed channel (offScale masked)
		y <- x[, channel ]
		if(isTRUE(maskOffScale)) y[ attr(x, "offScale") ] <- NA
		y <- ksmooth(x=1:length(y), y=y, kernel="normal", bandwidth=5, n.points=length(y))$y

		# Local maxima (above noise level)
		allPeaks <- which(head(tail(y,-1),-1) >= tail(y,-2) & head(y,-2) < head(tail(y,-1),-1) & head(tail(y,-1),-1) > noiseLevel) + 1L

		# Exclude outliers comparing intensity with sure peaks ('surePeaks' last ones)
		sureIntensity <- median(y[ tail(allPeaks, surePeaks) ])
		if(outThreshold < 1) { threshold <- outThreshold * sureIntensity
		} else               { threshold <- outThreshold
		}
		truePeaks <- allPeaks[ abs(y[allPeaks] - sureIntensity) < threshold ]
		
		if(length(truePeaks) > length(fullLadder)) {
			# Too many peaks retained (artefacts)
			if(trim == "forward")         {
				trueSizes <- fullLadder
				truePeaks <- head(truePeaks, length(fullLadder))
				warning("Detected more peaks than described in full size ladder, discarding last ones", call.=FALSE)
			} else if(trim == "backward") {
				trueSizes <- fullLadder
				truePeaks <- tail(truePeaks, length(fullLadder))
				warning("Detected more peaks than described in full size ladder, discarding first ones", call.=FALSE)
			} else {
				stop("Detected more peaks than described in full size ladder", call.=FALSE)
			}
		} else if(length(truePeaks) < length(fullLadder)) {
			# Not enough peaks retained (premature run ending)
			if(length(truePeaks) < length(useLadder)) stop("Detected less peaks than described in use ladder (", length(truePeaks), ")", call.=FALSE)
			if(trim == "forward")         { trueSizes <- head(fullLadder, length(truePeaks))
			} else if(trim == "backward") { trueSizes <- tail(fullLadder, length(truePeaks))
			} else                        { stop("Detected less peaks than described in full size ladder", call.=FALSE)
			}
		} else {
			# Fine
			trueSizes <- fullLadder
		}
		
		# Restrict modelization on interest zone
		usePeaks <- truePeaks[ trueSizes %in% useLadder ]
		useSizes <- trueSizes[ trueSizes %in% useLadder ]
	
		# Linear transformation of indexes in bp
		model <- lm(useSizes~usePeaks)
	
		# Check alignment
		rSquared <- summary(model)$adj.r.squared
		if(is.na(rSquared)) stop("Unable to compute R² (not enough points for a linear model ?)", call.=FALSE)
		if(rSquared < rMin) warning(sprintf("Poor alignment model (rSquared = %.5f)", rSquared), call.=FALSE)
	
		# Store model
		attr(x, "ladderModel") <- coefficients(model)
		attr(x, "ladderExact") <- usePeaks
		names(attr(x, "ladderExact")) <- useSizes
	})

	# Alignment rescue data
	if(isTRUE(rescue)) {
		# Replace raw intensities with smoothed ones (for rescue plot only)
		object <- x
		object[,channel] <- y
		
		# Plot profile
		if(is.null(ylim)) ylim <- c(min(y, sureIntensity-threshold, na.rm=TRUE), max(y, sureIntensity+threshold, na.rm=TRUE))
		plot(object, units=ifelse(is(status, "try-error"), "index", "bp"), ladder=!is(status, "try-error"), channels=channel, chanColors="#000000", ylim=ylim, ...)
		
		# All detected peaks
		xcor <- attr(object, "ladderModel")[2] * allPeaks + attr(object, "ladderModel")[1]
		axis(side=3, at=xcor, labels=FALSE, lwd.ticks=3, col.ticks="#BB3333")
		
		# Intensity filter
		xlim <- par()$usr[1:2]
		segments(x0=xlim[1], x1=xlim[2], y0=sureIntensity, y1=sureIntensity, lty="dotted")
		rect(xleft=xlim[1], xright=xlim[2], ybottom=sureIntensity-threshold, ytop=sureIntensity+threshold, col="#00000033", border=NA)
		
		# Peaks retained as size markers
		xcor <- attr(object, "ladderModel")[2] * truePeaks + attr(object, "ladderModel")[1]
		axis(side=3, at=xcor, labels=FALSE, lwd.ticks=3, col.ticks="#33BB33")
		
		# R squared
		if(!is.null(rSquared)) text(x=par("usr")[2], y=par("usr")[4], labels=substitute(R^2==value, list(value=round(rSquared, 6))), adj=c(1.25, 1.5))
	}
	
	# Call postponed errors
	if(is(status, "try-error")) stop(conditionMessage(attr(status, "condition")), call.=FALSE)
	
	return(x)
}
