# MLPA gene fusion TCL-TK worker
# Author : Sylvain Mareschal <maressyl@gmail.com>
fusions.process <- function(
		input,
		design,
		sheet = NA,
		output = ".",
		cores = NA,
		...
	) {
	
	requireNamespace("Biostrings")
	
	# Version
	message(date())
	message("MLPA ", as.character(packageVersion("MLPA")), sep="")
	message("Biostrings ", as.character(packageVersion("Biostrings")), "\n", sep="")
	
	# Check input
	if(!file.exists(input) || !file.info(input)$isdir) stop("Select an existing 'input' directory")
	
	# Auto-pick sample sheet
	if(is.na(sheet)) {
		 sheet <- dir(input, pattern="\\.csv$", full.names=TRUE)
		 if(length(sheet) != 1L) stop("'sheet' is NA but 'input' doesn't contain a single CSV file")
	}
	
	# Check sample sheet
	if(sheet == "" || !grepl("\\.(csv)$", sheet)) stop("Select a sample sheet with a \".csv\" extension.")
	if(!file.exists(sheet))                       stop("Select an existing sample sheet.")
	
	# Parse sample sheet
	samples <- read.csv(sheet, stringsAsFactors=FALSE)
	if(!setequal(colnames(samples), c("ID", "way", "file")))     stop("Columns in 'samples.csv' must be 'ID', 'way' and 'file'")
	if(nrow(samples) %% 2L != 0L | nrow(samples) == 0L)          stop("'samples.csv' must have a non-zero and even row count")
	if(any(!samples$way %in% c("forward", "reverse")))           stop("Column 'way' in 'samples.csv' must be 'forward' or 'reverse'")
	if(any(!file.exists(sprintf("%s/%s", input, samples$file)))) stop("Column 'file' in 'samples.csv' must refer to existing files")
	
	# Create output directory if necessary
	if(!file.exists(output)) dir.create(output)
	
	# Automatic core detection
	if(is.na(cores)) cores <- parallel::detectCores()
	
	if(cores == 1L) {
		# Un-parallelized processing
		out <- list()
		for(sampleName in unique(samples$ID)) {
			message("- ", sampleName)
			out[[sampleName]] <- fusions.process.core(
				sampleName = sampleName,
				samples = samples,
				design = design,
				input = input,
				output = output,
				...
			)
		}
	} else {
		# Parallelized processing
		cluster <- parallel::makeCluster(spec=cores)
		out <- parallel::clusterMap(
			cl = cluster,
			fun = fusions.process.handle,
			sampleName = unique(samples$ID),
			MoreArgs = list(
				samples = samples,
				design = design,
				input = input,
				output = output,
				...
			),
			RECYCLE = FALSE,
			SIMPLIFY = FALSE,
			USE.NAMES = FALSE,
			.scheduling = "dynamic"
		)
	}
	
	# Aggregate output
	dialogs <- do.call(rbind, lapply(out, "[[", "dialogs"))
	tab <- do.call(rbind, lapply(out, "[[", "output"))
	
	# Export output table
	write.csv2(dialogs, file=sprintf("%s/Dialogs.csv", output), row.names=FALSE)
	write.csv2(tab, file=sprintf("%s/Top_fusions.csv", output), row.names=FALSE)
	
	invisible(out)
}

