# MLPA Gene Expression Profiling TCL-TK interface
# Author : Sylvain Mareschal <maressyl@gmail.com>
GEP.interface <- function() {
	
	# Version name
	verName <- paste("MLPA", as.character(packageVersion("MLPA")), sep=" ")
	
	
	## FUNCTIONS ##
	
	inputBrowse <- function() {
		tcltk::tclvalue(inputValue) <- tk.folder(
			title = "Input directory",
			mustexist = TRUE,
			mandatory = TRUE
		)
		
		# Guess output
		if(tcltk::tclvalue(outputValue) == "") tcltk::tclvalue(outputValue) <- sprintf("%s.pdf", tcltk::tclvalue(inputValue))
		
		# Guess design
		if(tcltk::tclvalue(designValue) == "") {
			conf <- dir(tcltk::tclvalue(inputValue), recursive=TRUE, pattern=".+\\.conf$", full.names=TRUE)
			if(length(conf) == 1) tcltk::tclvalue(designValue) <- conf
		}
	}
	
	designBrowse <- function() {
		tcltk::tclvalue(designValue) <- tk.file(
			title = "Design file",
			typeNames = "Configuration file",
			typeExt = "*.conf",
			multiple = FALSE,
			mandatory = TRUE,
			type = "open",
			parent = topLevel
		)
	}

	outputBrowse <- function() {
		tcltk::tclvalue(outputValue) <- tk.file(
			title = "Output file names",
			typeNames = "Output log file",
			typeExt = "*.log",
			multiple = FALSE,
			mandatory = TRUE,
			type = "open",
			parent = topLevel
		)
	}
	
	process <- function() {
		# Watch cursor
		tcltk::tkconfigure(topLevel, cursor="watch")
		
		# Error-catching processing
		
		handle(
			expr = {
				# Log file
				logFile <- sub("\\.[^\\.]+$", ".log", tcltk::tclvalue(outputValue))
				warnCount <- 0L
				cat("", file=logFile)
				
				# CLI call
				GEP.process(
					input = tcltk::tclvalue(inputValue),
					design = tcltk::tclvalue(designValue),
					output = tcltk::tclvalue(outputValue),
					progressBar = progressBar
				)
				
				# Done
				if(warnCount == 0) { tcltk::tkmessageBox(parent=topLevel, parent=topLevel, icon="info", type="ok", title="Done", message="Processing achieved without warning.")
				} else             { tcltk::tkmessageBox(parent=topLevel, parent=topLevel, icon="warning", type="ok", title="Done", message=sprintf("Processing achieved with %i warning(s), please refer to the log file.", warnCount))
				}
			},
			messageHandler = function(m) {
				mes <- conditionMessage(m)
				if(grepl("^Processing", mes) || grepl("^Partial output", mes) || grepl("^All done", mes)) {
				         cat(sprintf("\n%-10s%s", "", mes), file=logFile, append=TRUE)
				} else { cat(sprintf("%-10s%s", "", mes), file=logFile, append=TRUE)
				}
			},
			warningHandler = function(w) {
				cat(sprintf("%-10s%s\n", "[WARNING]", conditionMessage(w)), file=logFile, append=TRUE)
				warnCount <<- warnCount + 1L
			},
			errorHandler = function(e) {
				cat(sprintf("%-10s%s\n", "[ERROR]", conditionMessage(e)), file=logFile, append=TRUE)
				tcltk::tkmessageBox(parent=topLevel, icon="error", type="ok", title="R error", message=sprintf("An error occured, check the log file for details.\n(%s)", conditionMessage(e)))
				dev.off()
			}
		)
		
		# Normal cursor
		tcltk::tkconfigure(topLevel, cursor="arrow")
	}
	
	replaceEntry = function(widget) {
		tcltk::tkfocus("-force", widget)
		tcltk::tcl(widget, "delete", "0", "end")
		tcltk::tcl(widget, "icursor", "end")
	}
	
	
	## INTERFACE ##
	
	# Top level
	topLevel <- tcltk::tktoplevel()
	tcltk::tktitle(topLevel) <- "RT-MLPA processor"
	
	# Horizontal resizing
	tcltk::tkgrid.columnconfigure(topLevel, 1, weight=1)
		
		# Preprocessing frame
		allFrame <- tcltk::ttklabelframe(parent=topLevel, relief="groove", borderwidth=2, text=verName)
		tcltk::tkgrid.columnconfigure(allFrame, 2, weight=1)
		
			# Input directory
			inputValue <- tcltk::tclVar("")
			inputButton <- tcltk::tkbutton(parent=allFrame, text="Input directory", command=inputBrowse, width=20)
			inputEntry <- tcltk::tkentry(parent=allFrame, textvariable=inputValue, width=70)
			tcltk::tkgrid(inputButton, column=1, row=1, padx=c(8,2), pady=c(8,2))
			tcltk::tkgrid(inputEntry, column=2, row=1, padx=c(2,8), pady=c(8,2), sticky="ew")
			
			# Design file
			designValue <- tcltk::tclVar("")
			designButton <- tcltk::tkbutton(parent=allFrame, text="Design file", command=designBrowse, width=20)
			designEntry <- tcltk::tkentry(parent=allFrame, textvariable=designValue, width=70)
			tcltk::tkgrid(designButton, column=1, row=2, padx=c(8,2), pady=2)
			tcltk::tkgrid(designEntry, column=2, row=2, padx=c(2,8), pady=2, sticky="ew")
			
			# Output files
			outputValue <- tcltk::tclVar("")
			outputButton <- tcltk::tkbutton(parent=allFrame, text="Output files", command=outputBrowse, width=20)
			outputEntry <- tcltk::tkentry(parent=allFrame, textvariable=outputValue, width=70)
			tcltk::tkgrid(outputButton, column=1, row=3, padx=c(8,2), pady=2)
			tcltk::tkgrid(outputEntry, column=2, row=3, padx=c(2,8), pady=2, sticky="ew")
			
			# Progression
			progressButton <- tcltk::tkbutton(parent=allFrame, text="Process", command=process, width=20)
			progressBar <- tcltk::ttkprogressbar(parent=allFrame, value=0, maximum=101, length=300)
			tcltk::tkgrid(progressButton, column=1, row=6, padx=c(8,2), pady=c(2,8))
			tcltk::tkgrid(progressBar, column=2, row=6, padx=c(2,8), pady=c(2,8), sticky="ew")
		
		tcltk::tkgrid(allFrame, column=1, row=1, padx=5, pady=5, sticky="nsew")
		
		# Help frame
		helpFrame <- tcltk::ttklabelframe(parent=topLevel, relief="groove", borderwidth=2, text="How to", width=100)
			
			# Help message
			helpMessage <- tcltk::tkmessage(parent=helpFrame, width=650, text="- Click buttons to select files and directories, or type their paths directly.\n- The \"Input directory\" must contain \".fsa\" files to process (sub-directories will be explored too).\n- A design file (.conf) stored with the data files (.fsa) will be detected automatically.")
			tcltk::tkgrid(helpMessage, column=1, row=1, padx=8, pady=8)
		
		tcltk::tkgrid(helpFrame, column=1, row=2, padx=5, pady=5, sticky="nsew")
		
	# Freeze size
	tcltk::tkwm.minsize(topLevel, 720, 260)
	tcltk::tkwm.resizable(topLevel, 1, 0)
	
	# Entry events
	tcltk::tkbind(inputEntry,  "<ButtonPress-3>", function(){ replaceEntry(inputEntry) })
	tcltk::tkbind(designEntry, "<ButtonPress-3>", function(){ replaceEntry(designEntry) })
	tcltk::tkbind(outputEntry, "<ButtonPress-3>", function(){ replaceEntry(outputEntry) })
	
	# Wait for closing
	tcltk::tkwait.window(topLevel)
}

