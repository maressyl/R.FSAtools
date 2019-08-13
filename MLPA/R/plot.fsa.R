# Plot method for "fsa" S3 class
# Author : Sylvain Mareschal <maressyl@gmail.com>
plot.fsa <- function(
		x,
		units = NA,
		channels = NA,
		chanColors = NA,
		ladder = TRUE,
		offScaleCol = "#FF0000",
		offScalePch = "+",
		offScaleCex = 0.4,
		bg = "white",
		fg = "black",
		title = "",
		xlab = NA,
		ylab = "Intensity",
		xlim = NA,
		ylim = NA,
		xaxt = "s",
		yaxt = "s",
		bty = "o",
		xaxp = NA,
		nticks = 5L,
		all.bp = TRUE,
		...
	) {
	# Defaults
	if(identical(channels, NA))   channels <- colnames(x)
	if(identical(chanColors, NA)) chanColors <- attr(x, "colors")
	if(identical(units, NA)) {
		if(is.null(attr(x, "ladderModel"))) { units <- "index"
		} else                              { units <- "bp"
		}
	}
	if(identical(xlab, NA))       xlab <- units
	
	# Checks
	if(length(channels) == 0) stop("No channel selected for plotting")
	
	# 'channels' must be a character vector
	if(is.numeric(channels)) channels <- colnames(x)[ channels ]
	
	# Missing channel(s)
	missing <- setdiff(channels, colnames(x))
	if(length(missing) > 0) stop("channel(s) ", paste(missing, collapse=", "), " not found (available: ", paste(colnames(x), collapse=", "), ")")
	
	# 'chanColors' must be a named character vector
	if(is.null(names(chanColors))) {
		# Unamed, supposed to be parallel with channels
		if(length(channels) != length(chanColors)) stop("'chanColors' must be named, or have the same length as 'channels'")
		names(chanColors) <- channels
	} else {
		# Already named, check if all are described
		if(any(! channels %in% names(chanColors))) stop("Some 'channels' elements can't be found in 'chanColors' names")
	}
	
	# Empty plot
	if(nrow(x) == 0 && (identical(xlim, NA) || identical(ylim, NA))) stop("Cannot plot an empty object without 'xlim' and 'ylim'")
	
	# X scale
	units <- match.arg(units, choices=c("index", "bp"))
	if(units == "index") {
		if(nrow(x) > 0) { xcor <- 1:nrow(x)
		} else          { xcor <- integer(0)
		}
		if(identical(xlim, NA)) xlim <- c(1, nrow(x))
		if(identical(xaxp, NA)) xaxp <- NULL
	} else if(units == "bp") {
		if(is.null(attr(x, "ladderModel"))) stop("Can't plot an unaligned object in 'bp' units")
		if(identical(xlim, NA)) xlim <- attr(x, "ladderModel")[2] * c(1, nrow(x)) + attr(x, "ladderModel")[1]
		if(identical(xaxp, NA)) xaxp <- c((xlim[1] %/% nticks)*nticks, (xlim[2] %/% nticks)*nticks, diff(xlim %/% nticks))
		if(nrow(x) > 0) { xcor <- attr(x, "ladderModel")[2] * 1:nrow(x) + attr(x, "ladderModel")[1]
		} else          { xcor <- integer(0)
		}
	}
	
	# ylim
	if(identical(ylim, NA)) {
		# Values to consider
		tmp <- as.vector(x[ xcor >= xlim[1] & xcor <= xlim[2] , channels ])
		tmp <- tmp[ !is.na(tmp) ]
		
		# Max in any
		if(length(tmp) > 0) { ylim <- c(0, max(tmp))
		} else              { ylim <- c(-1, 1)
		}
	}
	
	# Graphical parameters
	savePar <- par(bg=bg, fg=fg, col=fg, col.axis=fg, col.lab=fg, col.main=fg, col.sub=fg)
	on.exit(par(savePar))
	
	# Mask off-scale values
	offScaleValues <- x[ attr(x, "offScale") , , drop=FALSE ]
	x[ attr(x, "offScale") ,] <- NA
	
	# Plot channels
	first <- TRUE
	for(h in channels) {
		if(!first) par(new=TRUE)
		plot(
			x = xcor, y = x[,h],
			xlim = xlim, ylim = ylim,
			xlab = ifelse(first, xlab, ""),
			ylab = ifelse(first, ylab, ""),
			xaxt = ifelse(first, xaxt, "n"),
			yaxt = ifelse(first, yaxt, "n"),
			bty = ifelse(first, bty, "n"),
			col = chanColors[h],
			type="l", xaxp = xaxp,
			...
		)
		if(first) first <- FALSE
	}
	
	# Plot bp axis
	if(units == "bp" && isTRUE(all.bp)) axis(side=1, at=xlim[1]:xlim[2], labels=FALSE)
	
	# Off scale points
	if(length(attr(x, "offScale")) > 0) {
		# Coordinates
		xoff <- xcor[ attr(x, "offScale") ]
		yoff <- offScaleValues
		
		if(length(xoff) > 0) {
			for(h in channels) {
				# Color
				if(is.na(offScaleCol)) { koff <- chanColors[h]
				} else                 { koff <- offScaleCol
				}
				
				# Points
				points(x=xoff, y=yoff[,h], col=koff, type="p", pch=offScalePch, cex=offScaleCex)
			}
		}
	}
	
	# Exact ladder
	if(isTRUE(ladder)) {
		if(!is.null(attr(x, "ladderExact"))) {
			# Add ladder sizes
			at <- switch(units,
				"index" = attr(x, "ladderExact"),
				"bp" = attr(x, "ladderModel")[2] * attr(x, "ladderExact") + attr(x, "ladderModel")[1]
			)
			segments(x0=at, y0=par("usr")[4], x1=at, y1=par("usr")[3]-diff(par("usr")[3:4])/8, lty="dotted", xpd=NA)
			mtext(side=1, at=at, text=names(attr(x, "ladderExact")), line=2)
		} else warning("Can't add ladder without alignment ('ladderExact' attribute)")
	}
	
	# Legend
	inset <- c(par("din")[2]/1000, par("din")[1]/1000)
	legend("topleft", inset=inset, legend=colnames(x[, channels, drop=FALSE]), col=chanColors[channels], lty="solid", bg="#FFFFFF")
	
	# Title
	title(main=title, adj=0)
}

