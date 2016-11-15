#' ADMB setup
#' 
#' Sets up executable for the tpl file by looking for exe in package directory or
#' compiles tpl file in local directory (clean=F) of from package directory.
#' If admb cannot be found will attempt to run prepare_admb function (if exists) to
#' establish connections for compilation.
#' 
#' @param tpl character string for admb template file
#' @param compile if TRUE forces re-compilation of tpl file
#' @param clean if TRUE, deletes the tpl and executable files for amdb in local directory
#' @param re uses admb-re if TRUE for random effects
#' @param safe can be used to set safe mode for admb
#' @export 
setup_admb=function(tpl,compile=FALSE,clean=TRUE,safe=TRUE,re=FALSE)
{
if(R.Version()$os=="mingw32")
	ext=".exe"
else
	ext=""
sdir=system.file(package="marked")
# cleanup any leftover files
clean_admb(tpl)
# if argument clean is TRUE, delete exe and TPL files as well
if(clean)
{
	if(file.exists(paste(tpl,".tpl",sep=""))) unlink(paste(tpl,".tpl",sep=""))
	if(file.exists(paste(tpl,ext,sep=""))) unlink(paste(tpl,ext,sep=""))
}
# if tpl is not available, copy from the package directory
if(!file.exists(paste(tpl,".tpl",sep="")))
{
	file.copy(file.path(sdir,paste(tpl,".tpl",sep="")),file.path(getwd(),paste(tpl,".tpl",sep="")),overwrite=TRUE)
	file.copy( file.path(sdir,list.files(sdir,"*.cpp")),file.path(getwd()),overwrite=TRUE)
} 
# if exe available in the package or workspace directory use it
have.exe=TRUE
if(!file.exists(file.path(getwd(),paste(tpl,ext,sep=""))))
{
	if(file.exists(file.path(sdir,paste(tpl,ext,sep=""))) &!compile)
	{
		file.copy(file.path(sdir,paste(tpl,ext,sep="")),file.path(getwd(),paste(tpl,ext,sep="")),overwrite=TRUE)
	}else
	# if there is no exe in either place then check to make sure ADMB is installed and linked
	{
		have.exe=FALSE
	}
}
if(!have.exe | compile)
{
	# no exe or compile set TRUE; see if admb can be found; this is not a complete test but should catch the novice user who has
	# not setup admb at all
	if(R.Version()$os=="mingw32" & Sys.which(paste("tpl2cpp",ext,sep=""))=="")
	    stop("admb not found; setup links to admb and c++ compiler with environment variables or put in path")
	else
	{
		compile_admb(tpl,re=re,safe=safe,verbose=T)
	}
}
invisible()
}
#' TMB setup
#' 
#' Sets up executable for the .cpp file (tpl) by looking for exe in package directory or
#' compiles cpp file in local directory (clean=FALSE) of from package directory.
#' 
#' @param tpl character string for admb template file
#' @param clean if TRUE, deletes the tpl (.cpp) and executable files in local directory and copies from package directory
#' @export 
setup_tmb=function(tpl,clean=FALSE)
{
	sdir=system.file(package="marked")
#   if argument clean is TRUE, delete dll and cpp files as well
	if(clean)
	{
		cat("\n deleting old TMB program and rebuilding")
		if(file.exists(paste(tpl,".cpp",sep=""))) unlink(paste(tpl,".cpp",sep=""))
		if(is.loaded(dynlib(tpl)))dyn.unload(dynlib(tpl))
		if(file.exists(paste(tpl,".dll",sep=""))) unlink(paste(tpl,".dll",sep=""))
		file.copy(file.path(sdir,paste(tpl,".cpp",sep="")),file.path(getwd(),paste(tpl,".cpp",sep="")),overwrite=TRUE)
		cat("\n compiling and linking TMB program\n")
		compile(paste(tpl,".cpp",sep=""))               # Compile the C++ file
		dyn.load(dynlib(tpl))                           # Dynamically link the C++ code
	} else
	{
# if cpp is not available, copy from the package directory
		if(!file.exists(paste(tpl,".cpp",sep="")))
		{
			file.copy(file.path(sdir,paste(tpl,".cpp",sep="")),file.path(getwd(),paste(tpl,".cpp",sep="")),overwrite=TRUE)
			if(is.loaded(dynlib(tpl)))dyn.unload(dynlib(tpl))
			if(file.exists(paste(tpl,".dll",sep=""))) unlink(paste(tpl,".dll",sep=""))
			cat("\n compiling and linking TMB program\n")
			compile(paste(tpl,".cpp",sep=""))               # Compile the C++ file
			dyn.load(dynlib(tpl))          # Dynamically link the C++ code
		} else
		{
			if(file.exists(paste(tpl,".dll",sep="")))
			{
				if(!is.loaded(dynlib(tpl))) {
					dyn.load(dynlib(tpl))
				}
			} else
			{
				if(file.exists(paste(tpl,".o",sep=""))) unlink(paste(tpl,".o",sep=""))
				cat("\n compiling and linking TMB program\n")
				compile(paste(tpl,".cpp",sep=""))               # Compile the C++ file
				if(is.loaded(dynlib(tpl)))dyn.unload(dynlib(tpl))
				dyn.load(dynlib(tpl))          # Dynamically link the C++ code
			}	
		}
	  }
	  invisible()
}
