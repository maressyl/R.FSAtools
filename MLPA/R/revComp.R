# Revserse complement a vector of character strings (ATCG)
# Author : Sylvain Mareschal <maressyl@gmail.com>
revComp <- function(seq) {
	out <- character(length(seq))
	seq <- strsplit(seq, split="")
	dic <- c("A"="T", "C"="G", "G"="C", "T"="A")
	for(i in 1:length(seq)) out[i] <- paste(dic[ rev(seq[[i]]) ], collapse="")
	return(out)
}
