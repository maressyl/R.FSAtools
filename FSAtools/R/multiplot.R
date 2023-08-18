# Wrapper for layout with only atomic arguments
multiplot <- function(
			nrow,
			ncol,
			widths = rep.int(1, ncol),
			heights = rep.int(1, nrow),
			indexes = 1:(nrow*ncol),
			byrow = FALSE,
			respect = FALSE,
			cex = 1
		) 
	{
	
	# Transfer to layout
	layout(
		mat = matrix(indexes, nrow=nrow, ncol=ncol, byrow=byrow),
		widths = widths,
		heights = heights,
		respect = respect
	)
	if(!is.na(cex)) par(cex=cex)
	
	invisible(TRUE)
}

