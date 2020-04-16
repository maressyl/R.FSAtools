# Plot method for "fsaModel" S3 class
# Author : Sylvain Mareschal <maressyl@gmail.com>
plot.fsaModel <- function(
		x,
		xlab = "Score",
		lwd = 3,
		...
	) {
	# Plot values
	xlim <- c(min(x$groupMeans - 4*x$groupSDs), max(x$groupMeans + 4*x$groupSDs))
	xval <- seq(from=xlim[1], to=xlim[2], length.out=1000)
	yval1 <- dnorm(xval, mean=x$groupMeans[1], sd=x$groupSDs[1])
	yval2 <- dnorm(xval, mean=x$groupMeans[2], sd=x$groupSDs[2])
	ylim <- range(c(yval1, yval2), na.rm=TRUE)
	ylim[2] <- ylim[2] * 1.1
	
	# Plot gaussians
	plot(x=xval, y=yval1, type="l", col=x$groupColors[1], xlim=xlim, ylim=ylim, xlab=xlab, ylab="", yaxt="n", lwd=lwd, ...)
	lines(x=xval, y=yval2, col=x$groupColors[2], lwd=3)
	
	# Estimate unclassified regions
	P1 <- yval1 / (yval1 + yval2)
	P2 <- yval2 / (yval1 + yval2)
	grayPoints <- rle(P1 < x$threshold & P2 < x$threshold)
	grayBreaks <- c(1, cumsum(grayPoints$lengths))
	
	# Plot
	for(i in which(grayPoints$values)) {
		# X coordinates
		left <- xval[ grayBreaks[i] ]
		right <- xval[ grayBreaks[i+1] ]
		
		# Plot
		rect(xleft=left, xright=right, ybottom=-100, ytop=100, col="#00000030", border=NA)
	}
	
	# Legend
	legend(x="topleft", lty="solid", col=x$groupColors, legend=x$groupNames, bty="n", lwd=lwd)
}

