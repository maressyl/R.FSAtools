plot.fusions <- function(sampleName, forward, reverse, forwardFile, reverseFile, identified, design, pal=NULL) {
	# Top length
	n <- nrow(identified$top)
	
	# Default palette
	if(is.null(pal)) {
		pal <- c(
			"unileft" = "red",
			"left.seq" = "orange",
			"right.seq" = "darkgreen",
			"uniright" = "blue",
			"M13" = "darkgrey"
		)
	}
	
	# Layout
	layout(matrix(1:(n*2+2), ncol=1), heights=c(5, rep(1, n), 5, rep(1, n)))
	
	# Iterable objects
	objects <- list(
		forward = forward,
		reverse = reverse
	)
	
	for(way in c("forward", "reverse")) {
		# Shortcuts
		seq <- attr(objects[[way]], "seq")
		peaks <- attr(objects[[way]], "peaks")

		# Plot the profile
		par(mar=c(1,4,1,10), cex=1)
		plot.sanger(objects[[way]], xaxt="n", xlab="", ladder=FALSE, lwd=2, las=2, bty="n")
		mtext(side=1, line=0.1, col="black", at=peaks, text=seq)
		mtext(side=4, las=1, line=1, at=0, text="Align", font=2, adj=0.5)
		mtext(side=4, las=1, line=5, at=0, text="Phred", font=2, adj=0.5)
		mtext(side=4, las=1, line=8, at=0, text="Cover", font=2, adj=0.5)
		
		# Meta-data
		at <- par("usr")[1] + 0.1 * diff(par("usr")[1:2])
		mtext(side=3, line=-1, adj=0, cex=1, at=at, text=sprintf("%s (%s)", sampleName, way), font=2)
		mtext(side=3, line=-2, adj=0, cex=0.7, at=at, text=ifelse(way == "forward", forwardFile, reverseFile))
		mtext(side=3, line=-2.8, adj=0, cex=0.7, at=at, text=sprintf("Run started on %s", attr(objects[[way]], "runMetaData")$runStart))
		mtext(side=3, line=-3.6, adj=0, cex=0.7, at=at, text=sprintf("Machine : %s", attr(objects[[way]], "runMetaData")$machine))
		mtext(side=3, line=-4.4, adj=0, cex=0.7, at=at, text=sprintf("Module : %s %s", attr(objects[[way]], "runMetaData")$runModule.name, attr(objects[[way]], "runMetaData")$runModule.version))
		mtext(side=3, line=-5.2, adj=0, cex=0.7, at=at, text=sprintf("Protocole : %s %s", attr(objects[[way]], "runMetaData")$runProtocole.name, attr(objects[[way]], "runMetaData")$runProtocole.version))
		
		for(match in identified$top$i) {
			# Background
			par(mar=c(0,4,0,10))
			plot(x=NA, y=NA, xlim=par("usr")[1:2], ylim=c(0,4), xlab="", ylab="", xaxt="n", yaxt="n", xaxs="i", yaxs="i", bty="n")
			
			# Peak annotation
			ann <- annotate.alignment(object=objects[[way]], align=identified$alignments[[way]][match], way=way, design_row=design[match,])
			
			# Simplify insertions
			i <- 1L
			while(i <= nrow(ann)) {
				if(is.na(ann[i,"x"])) {
					# NA range
					start <- i
					while(is.na(ann[i,"x"]) & i < nrow(ann)) i <- i +1L
					end <- i - 1L
					
					# Simplify
					ann[start,"x"] <- (ann[start-1,"x"] + ann[end+1,"x"]) / 2
					ann[start:end,"alignment"] <- NA
					ann[start,"alignment"] <- if(start == end) { "+" } else { sprintf("+\n%i", end - start + 1L) }
				} else {
					# Not a NA, pass
					i <- i + 1L
				}
			}
			ann <- ann[ !is.na(ann$x) ,]
			
			# Clean reference sequence
			text(x=ann$x, y=3, col=pal[ ann$group ], labels=ann$alignment, cex=0.8)
			
			for(element in unique(ann$group)) {
				# Underline the component
				x0 <- min(ann[ ann$group == element, "x" ])
				x1 <- max(ann[ ann$group == element, "x" ])
				segments(x0=x0, x1=x1, y0=2, y1=2, col=pal[element], xpd=NA)
				
				# Add a meaningful title
				title <- element
				if(title == "left.seq")  title <- sprintf("%s (%g%%)", identified$top[ identified$top$i == match , "left.name" ],  round(100*identified$top[ identified$top$i == match , sprintf("coverage.left.%s", way) ]))
				if(title == "right.seq") title <- sprintf("%s (%g%%)", identified$top[ identified$top$i == match , "right.name" ], round(100*identified$top[ identified$top$i == match , sprintf("coverage.right.%s", way) ]))
				
				# Print element title
				text(x=(x0+x1)/2, y=1, labels=title, col=pal[element], xpd=NA)
			}

			# Alignment score
			class <- "OK"
			if(identified$top[ identified$top$i == match, "score.sum" ] < 300) class <- "warn"
			if(identified$top[ identified$top$i == match, "score.sum" ] < 150) class <- "fail"
			mtext(
				side=4, las=1, line=1, adj=0.5,
				col = switch(class, OK="darkgreen", warn="orange", fail="red"),
				text = sprintf(
					"%g / %g",
					round(identified$top[ identified$top$i == match, sprintf("score.%s", way) ]),
					round(identified$top[ identified$top$i == match, "score.sum" ])
				)
			)
			
			# Phred score
			class <- "OK"
			if(identified$top[ identified$top$i == match, sprintf("phred.%s", way) ] < 35) class <- "warn"
			if(identified$top[ identified$top$i == match, sprintf("phred.%s", way) ] < 10) class <- "fail"
			mtext(
				side=4, las=1, line=5, adj=0.5,
				col = switch(class, OK="darkgreen", warn="orange", fail="red"),
				text = round(identified$top[ identified$top$i == match, sprintf("phred.%s", way) ], 1)
			)
			
			# Minimum cover
			left.cover <- max(
				identified$top[ identified$top$i == match, "coverage.left.forward" ],
				identified$top[ identified$top$i == match, "coverage.left.reverse" ]
			)
			right.cover <- max(
				identified$top[ identified$top$i == match, "coverage.right.forward" ],
				identified$top[ identified$top$i == match, "coverage.right.reverse" ]
			)
			class <- "OK"
			symbol <- c("\U2714", "\U2714")
			if(left.cover < 0.8)  { symbol[1] <- "!"; class <- "warn" }
			if(right.cover < 0.8) { symbol[2] <- "!"; class <- "warn" }
			if(left.cover < 0.2)  { symbol[1] <- "\U2718"; class <- "fail" }
			if(right.cover < 0.2) { symbol[2] <- "\U2718"; class <- "fail" }
			mtext(
				side=4, las=1, line=8, adj=0.5,
				col = switch(class, OK="darkgreen", warn="orange", fail="red"),
				text = paste(symbol, collapse="")
			)
		}
	}
}

