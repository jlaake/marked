print.marked.version <- function()
{ library(help=marked)$info[[1]] -> version
	version <- version[pmatch("Version",version)]
	um <- strsplit(version," ")[[1]]
	version <- um[nchar(um)>0][2]
	hello <- paste("This is marked ",version,"\n",sep="")
	packageStartupMessage(hello)
}
.onAttach <- function(...) { 
	print.marked.version()
}