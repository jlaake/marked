\name{What's new}
\alias{What's new}
\title{Summary of changes by version}
\description{A good place to look for new changes.  Often I'll add changes here but don't always get to it in the documentation 
for awhile.  They are ordered from newest to oldest.}
\details{
Version 1.0.3 ( May 2011)
\itemize{
\item Code modularization to make package easier to extend.
\item If \code{time.interval} is a field in the Phi design data it
will be used as Phi^time.interval.  This allows variation in time intervals
across animals as with cohorts of pups branded at different times each year.

}
Version 1.0.2 ( March 2011)
\itemize{
\item Fixed issue with capture history accumulator when chunk_size was larger than needed. Added 1 to pieces.
\item Added check on accumulator to report error if sum of frequencies does not match original number of records.
\item Added model convergence check and reporting of model convergence message, if non-null. 
\item Beta estimate names were lost once sparse matrices were implemented and this was fixed.
\item Added function \code{\link{fix.parameters}} to create matrix needed for fixing real parameters
\item Made change to accumulation code to correct error introduced in last version with fixed parameters.
\item Added option to include \code{remove.intercept=TRUE} in \code{model.parameters} list for each parameter to force removal of intercept.
\item Added \code{refit} argument to \code{crm},\code{cjs} and \code{js} to control number of refittings if model doesn't converge.
\item Added function \code{create.links} which works out which real parameters can use a sin link.  Code was added to \code{\link{cjs}} and 
       \code{\link{cjs.lnl}} to use the sin link where appropriate.  This is commented out at present!  Not sure it is working correctly.
\item Uses optimx function for optimization which allows more methods and multiple methods to be selected.
\item Made fixes to \code{\link{js}} to accommodate accumulation.
}    
Version 1.0.1 (22 March 2011)
\itemize{
\item Added use of sparse matrices for design matrices which sped up code and reduced memory consumption. The argument chunk_size was
added to \code{\link{crm}}, \code{\link{create.dm}},\code{\link{cjs}}, and \code{\link{js}} to control amount of memory usage.
\item Added run timing and various print statements to track progress of model.  If debug=FALSE, includes function evaluation counter (every 100) 
and neg lnl which remains on same line unless used in Rterm. 
\item Note that js model needs further testing after these changes.
}    
Version 1.0.0 (Initial posting 2011)
\itemize{
\item Extracted crm and accompanying code from RMark and created base package which was posted on bitbucket.
}
}
\author{Jeff Laake}

