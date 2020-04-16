# Imports a .ab1 / .fsa file from Applied Biosystems, using the provided converter
# Author : Sylvain Mareschal <maressyl@gmail.com>
identify.fusions <- function(forward, reverse, design, top=5, extra=NULL) {
	# Check arguments
	if(!is.data.frame(design)) stop("'design' must be a data.frame")
	
	# Aggregate objects
	objects <- list(
		forward = forward,
		reverse = reverse
	)
	
	# Pairwise alignment using quality scores
	align <- list()
	for(way in c("forward", "reverse")) {
		if(!"sanger" %in% class(objects[[way]])) stop("'", way, "' must be a 'sanger' object")
		align[[way]] <- Biostrings::pairwiseAlignment(
			pattern = design[[ sprintf("seq_%s", way) ]],
			patternQuality = Biostrings::PhredQuality(65L),
			subject = paste(attr(objects[[way]], "seq"), collapse=""),
			subjectQuality = Biostrings::PhredQuality(attr(objects[[way]], "phred")),
			type = "local"
		)
	}

	# Best candidate(s) in both ways
	score <- Biostrings::score(align$forward) + Biostrings::score(align$reverse)
	selection <- head(order(score, decreasing=TRUE), top)

	# Aggregate annotation
	out <- design[ selection , c("left.name", "left.symbol", "left.GRCh38", "left.GRCh38_band", "right.name", "right.symbol", "right.GRCh38", "right.GRCh38_band", extra) ]
	out$score.forward <- Biostrings::score(align$forward[selection])
	out$score.reverse <- Biostrings::score(align$reverse[selection])
	out$score.sum <- score[selection]
	out$i <- selection
	
	# Aligned part of the reference (forward)
	sequenced.forward.start <- Biostrings::start(Biostrings::pattern(align$forward[selection]))
	sequenced.forward.end   <- Biostrings::end(Biostrings::pattern(align$forward[selection]))
	
	# Aligned part of the reference (reverse)
	sequenced.reverse.start <- Biostrings::start(Biostrings::pattern(align$reverse[selection]))
	sequenced.reverse.end   <- Biostrings::end(Biostrings::pattern(align$reverse[selection]))
	
	# Specific part of the reference (forward)
	specific.forward.left.start  <- nchar(design[selection,"left.unileft"])
	specific.forward.left.end    <- specific.forward.left.start + nchar(design[selection,"left.seq"])
	specific.forward.right.start <- specific.forward.left.end + 1L
	specific.forward.right.end   <- specific.forward.right.start + nchar(design[selection,"right.seq"]) - 1L
	
	# Specific part of the reference (reverse)
	specific.reverse.right.start <- nchar(design[selection,"right.uniright"])
	specific.reverse.right.end   <- specific.reverse.right.start + nchar(design[selection,"right.seq"])
	specific.reverse.left.start  <- specific.reverse.right.end + 1L
	specific.reverse.left.end    <- specific.reverse.left.start + nchar(design[selection,"left.seq"]) - 1L
	
	# Proportion sequenced
	out$coverage.left.forward  <- pmax(0, (pmin(sequenced.forward.end, specific.forward.left.end)  - pmax(sequenced.forward.start, specific.forward.left.start)  + 1L) / (specific.forward.left.end  - specific.forward.left.start  + 1L))
	out$coverage.left.reverse  <- pmax(0, (pmin(sequenced.reverse.end, specific.reverse.left.end)  - pmax(sequenced.reverse.start, specific.reverse.left.start)  + 1L) / (specific.reverse.left.end  - specific.reverse.left.start  + 1L))
	out$coverage.right.forward <- pmax(0, (pmin(sequenced.forward.end, specific.forward.right.end) - pmax(sequenced.forward.start, specific.forward.right.start) + 1L) / (specific.forward.right.end - specific.forward.right.start + 1L))
	out$coverage.right.reverse <- pmax(0, (pmin(sequenced.reverse.end, specific.reverse.right.end) - pmax(sequenced.reverse.start, specific.reverse.right.start) + 1L) / (specific.reverse.right.end - specific.reverse.right.start + 1L))
	
	# Average Phred over the aligned region (forward)
	for(way in c("forward", "reverse")) {
		out[[ sprintf("phred.%s", way) ]] <- as.double(NA)
		for(i in 1:top) {
			aligned.range <- Biostrings::start(Biostrings::subject(align[[way]][selection][i])) : Biostrings::end(Biostrings::subject(align[[way]][selection][i]))
			phred <- as.integer(align[[way]][selection][i]@subject@unaligned@quality)
			out[ i , sprintf("phred.%s", way) ] <- mean(phred[ aligned.range ])
		}
	}
	
	# Flags
	negative <- cbind(
		"WT NPM1"        = out$left.name == "FwtNPM1" & out$right.name == "RmutNPM1",
		"Phred < 10"     = out$phred.forward < 10 | out$phred.reverse < 10,
		"Score < 150"    = out$score.sum < 150,
		"Coverage < 20%" = (out$coverage.left.forward < 0.2 & out$coverage.left.reverse < 0.2) | (out$coverage.right.forward < 0.2 & out$coverage.right.reverse < 0.2)
	)
	warning <- cbind(
		"Phred < 35"     = out$phred.forward < 35 | out$phred.reverse < 35,
		"Score < 300"    = out$score.sum < 300,
		"Coverage < 80%" = (out$coverage.left.forward < 0.8 & out$coverage.left.reverse < 0.8) | (out$coverage.right.forward < 0.8 & out$coverage.right.reverse < 0.8)
	)
	out$negative <- apply(negative, 1, function(x) { paste(colnames(negative)[x], collapse=", ") })
	out$warning  <- apply(warning,  1, function(x) { paste(colnames(warning)[x],  collapse=", ") })
	
	return(list(top=out, alignments=align))
}

