# Imports a .ab1 / .fsa file from Applied Biosystems, using the provided converter
# Author : Sylvain Mareschal <maressyl@gmail.com>
annotate.alignment <- function(object, align, way, design_row) {
	# Check arguments
	if(!is(align, "PairwiseAlignmentsSingleSubject")) stop("'align' must be a Biostring's 'PairwiseAlignmentsSingleSubject' object")
	
	# Peak positions in the electrophoregram
	peaks <- attr(object, "peaks")
	
	# Aligned sequences
	aln.pattern <- strsplit(as.character(Biostrings::pattern(align)), split="")[[1]]
	aln.subject <- strsplit(as.character(Biostrings::subject(align)), split="")[[1]]

	# Ranges occupied by the reference pattern and the sequenced subject (local alignment)
	i.pattern <- Biostrings::start(Biostrings::pattern(align)) : Biostrings::end(Biostrings::pattern(align))
	i.subject <- Biostrings::start(Biostrings::subject(align)) : Biostrings::end(Biostrings::subject(align))

	# Deletion (in the subject / sequence obtained)
	if(any(aln.subject == "-")) {
		i.subject <- i.subject[ cumsum(aln.subject != "-") ]
		i.subject[ aln.subject == "-" ] <- NA
	}

	# Insertion (in the subject / sequence obtained)
	if(any(aln.pattern == "-")) {
		i.pattern <- i.pattern[ cumsum(aln.pattern != "-") ]
		i.pattern[ aln.pattern == "-" ] <- NA
	}

	# Split reference sequence into components
	if(way == "forward") {
		k <- c(
			rep("unileft", nchar(design_row$left.unileft)),
			rep("left.seq", design_row$left.seq.size),
			rep("right.seq", design_row$right.seq.size),
			rep("uniright", nchar(design_row$right.uniright)),
			rep("M13", nchar(design_row$right.M13.right))
		)
	} else {
		k <- c(
			rep("uniright", nchar(design_row$right.uniright)),
			rep("right.seq", design_row$right.seq.size),
			rep("left.seq", design_row$left.seq.size),
			rep("unileft", nchar(design_row$left.unileft)),
			rep("M13", nchar(design_row$left.M13.left))
		)
	}

	# Trim to local alignment window
	k <- k[i.pattern]
	k[ is.na(k) ] <- ""
	
	# Annotate each peak
	ann <- data.frame(
		x = peaks[i.subject],
		group = k,
		alignment = ifelse(aln.subject == "-", sprintf("+%s", aln.pattern), ifelse(aln.pattern == aln.subject, "|", aln.pattern)),
		stringsAsFactors = FALSE
	)
	
	return(ann)
}

