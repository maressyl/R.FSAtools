# Plot method for "sanger" S3 class
# Author : Sylvain Mareschal <maressyl@gmail.com>
plot.sanger <- function(
		x,
		xaxt = "n",
		xlab = "",
		ladder = FALSE,
		lwd = 2,
		xlim = NA,
		...
	) {
	# Auto-xlim
	if(identical(xlim, NA)) {
		if(nrow(x) == 0) stop("Cannot plot an empty object without 'xlim'")
		xlim <- c(0, max(attr(x, "peaks")))
	}
	
	# Profile
	plot.fsa(x=x, xaxt=xaxt, xlab=xlab, ladder=ladder, lwd=lwd, xlim=xlim, ...)
	
	# Calls	
	mtext(side=1, line=0.1, col="black", at=peaks <- attr(x, "peaks"), text=attr(x, "seq"))
}

