# Interactive tcltk::tcl-TK directory choosing
tk.folder = function(
		title = "Choose a directory",
		mustexist = TRUE,
		mandatory = TRUE
		)
	{
	# Dialog
	suppressWarnings(folder <- tcltk::tclvalue(tcltk::tkchooseDirectory(title=title, mustexist=mustexist)))
	
	# No folder
	if(mandatory && folder == "") stop("No directory selected")
	
	return(folder)
}

