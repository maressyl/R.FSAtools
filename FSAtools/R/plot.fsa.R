# Plot method for "fsa" S3 class
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
		title.adj = 0,
		title.line = NA,
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
		peaks.alpha = 48L,
		peaks.srt = 30,
		peaks.adj = c(0, 0),
		peaks.cex = 1.3,
		peaks.font = 2,
		legend.x = "topleft",
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
	
	# Background
	plot(
		x = NA, y = NA,
		xlim = xlim, ylim = ylim,
		xlab = xlab,
		ylab = ylab,
		xaxt = xaxt,
		yaxt = yaxt,
		bty = bty,
		xaxp = xaxp,
		...
	)
	
	# Peaks
	peaks <- attr(x, "peaks")
	if(!is.null(peaks)) {
		# Ignore invisible peaks
		peaks <- peaks[ !is.na(peaks$color) ,]
		
		# Transparent version of peak colors
		peaks$background <- sprintf(
			"#%s",
			apply(
				rbind(
					col2rgb(peaks$color),
					peaks.alpha
				),
				2,
				function(x) {
					paste(sprintf("%02x", x), collapse="")
				}
			)
		)
		
		# Full height rectangles
		rect(
			xleft = peaks$size.min,
			xright = peaks$size.max,
			ybottom = -1e6,
			ytop = 1e6,
			col = peaks$background,
			border = NA
		)
		
		if(all(c("N0", "N1", "N2") %in% colnames(peaks))) {
			# Allele rectangles
			rect(
				xleft = peaks$size.min,
				xright = peaks$size.max,
				ybottom = c(peaks$N0, peaks$N0) * peaks$height / peaks$normalized,
				ytop = c(peaks$N1, peaks$N2) * peaks$height / peaks$normalized,
				col = peaks$background,
				border = NA
			)
		}
		
		# Names
		text(
			x = (peaks$size.min + peaks$size.max) / 2,
			y = par("usr")[4] + diff(par("usr")[3:4]) / 50,
			labels = rownames(peaks),
			srt = peaks.srt,
			adj = peaks.adj,
			cex = peaks.cex,
			font = peaks.font,
			xpd = NA,
			col = peaks$color
		)
		
		# Genotyping N1-based
		if("present" %in% colnames(peaks)) {
			text(
				x = (peaks$size.min + peaks$size.max) / 2,
				y = par("usr")[4],
				adj = c(0.5, 1),
				col = chanColors[ peaks$channel ],
				labels = ifelse(peaks$present, "+", "-"),
				font = 2
			)
		}
	}
	
	# Plot channels
	for(h in channels) {
		points(
			x = xcor,
			y = x[,h],
			col = chanColors[h],
			type = "l"
		)
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
	legend(legend.x, inset=inset, legend=colnames(x[, channels, drop=FALSE]), col=chanColors[channels], lty="solid", bg="#FFFFFF")
	
	# Title
	title(main=title, adj=title.adj, line=title.line)
	
	invisible(TRUE)
}

