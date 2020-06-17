# Plot method for "fsa" S3 class
# Author : Sylvain Mareschal <maressyl@gmail.com>
addPeaks <- function(
		ranges,
		colors,
		backgrounds,
		border = NA,
		srt = 30,
		adj = c(0, 0),
		cex = 1.3,
		font = 2,
		...
	) {

	# Rectangles
	rect(
		xleft = sapply(ranges, "[", 1),
		xright = sapply(ranges, "[", 2),
		ybottom = -1e6,
		ytop = 1e6,
		col = backgrounds,
		border = border
	)
	
	# Names
	text(
		x = sapply(ranges, mean),
		y = par("usr")[4] + diff(par("usr")[3:4]) / 50,
		labels = names(ranges),
		srt = srt,
		adj = adj,
		cex = cex,
		font = font,
		xpd = NA,
		col = colors
	)
	
	# Retrace box
	box()
}

