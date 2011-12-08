

#' Dipper capture-recapture data
#' 
#' A capture-recapture data set on European dippers from France that
#' accompanies MARK as an example analysis using the CJS and POPAN models.  The
#' dipper data set was orginally described as an example by Lebreton et al
#' (1992).
#' 
#' This is a data set that accompanies program MARK as an example for CJS and
#' POPAN analyses.  The data can be stratified using sex as a grouping
#' variable.  The functions \code{run.dipper}, \code{run.dipper.alternate},
#' \code{run.dipper.popan} defined below in the examples mimic the models used
#' in the dbf file that accompanies MARK. Note that the models used in the MARK
#' example use PIM coding with the sin link function which is often better at
#' identifying the number of estimable parameters.  The approach used in the R
#' code uses design matrices and cannot use the sin link and is less capable at
#' counting parameters.  These differences are illustrated by comparing the
#' results of \code{run.dipper} and \code{run.dipper.alternate} which fit the
#' same set of "CJS" models.  The latter fits the models with constraints on
#' some parameters to achieve identifiability and the former does not. Although
#' it does not influence the selection of the best model it does infleunce
#' parameter counts and AIC ordering of some of the less competitive models. In
#' using design matrices it is best to constrain parameters that are confounded
#' (e.g., last occasion parameters in Phi(t)p(t) CJS model) when possible to
#' achieve more reliable counts of the number of estimable parameters.
#' 
#' Note that the covariate "sex" defined in dipper has values "Male" and
#' "Female".  It cannot be used directly in a formula for MARK without using it
#' do define groups because MARK.EXE will be unable to read in a covariate with
#' non-numeric values.  By using \code{groups="sex"} in the call the
#' \code{\link{process.data}} a factor "sex" field is created that can be used
#' in the formula.  Alternatively, a new covariate could be defined in the data
#' with say values 0 for Female and 1 for Male and this could be used without
#' defining groups because it is numeric.  This can be done easily by
#' translating the values of the coded variables to a numeric variable.  Factor
#' variables are numbered 1..k for k levels in alphabetic order.  Since Female
#' < Male in alphabetic order then it is level 1 and Male is level 2.  So the
#' following will create a numeric sex covariate.
#' 
#' \preformatted{ dipper$numeric.sex=as.numeric(dipper$sex)-1 }
#' 
#' @name dipper
#' @docType data
#' @format A data frame with 294 observations on the following 2 variables.
#' \describe{ \item{ch}{a character vector containing the encounter history of
#' each bird} \item{sex}{the sex of the bird: a factor with levels
#' \code{Female} \code{Male}} }
#' @source Lebreton, J.-D., K. P. Burnham, J. Clobert, and D. R. Anderson.
#' 1992. Modeling survival and testing biological hypotheses using marked
#' animals: case studies and recent advances. Ecol. Monogr. 62:67-118.
#' @keywords datasets
NULL





#' Summary of changes by version
#' 
#' A good place to look for new changes.  Often I'll add changes here but don't
#' always get to it in the documentation for awhile.  They are ordered from
#' newest to oldest.
#' 
#' Version 1.0.3 ( May 2011) \itemize{ \item Code modularization to make
#' package easier to extend. \item If \code{time.interval} is a field in the
#' Phi design data it will be used as Phi^time.interval.  This allows variation
#' in time intervals across animals as with cohorts of pups branded at
#' different times each year.
#' 
#' } Version 1.0.2 ( March 2011) \itemize{ \item Fixed issue with capture
#' history accumulator when chunk_size was larger than needed. Added 1 to
#' pieces. \item Added check on accumulator to report error if sum of
#' frequencies does not match original number of records. \item Added model
#' convergence check and reporting of model convergence message, if non-null.
#' \item Beta estimate names were lost once sparse matrices were implemented
#' and this was fixed. \item Added function \code{\link{fix.parameters}} to
#' create matrix needed for fixing real parameters \item Made change to
#' accumulation code to correct error introduced in last version with fixed
#' parameters. \item Added option to include \code{remove.intercept=TRUE} in
#' \code{model.parameters} list for each parameter to force removal of
#' intercept. \item Added \code{refit} argument to \code{crm},\code{cjs} and
#' \code{js} to control number of refittings if model doesn't converge. \item
#' Added function \code{create.links} which works out which real parameters can
#' use a sin link.  Code was added to \code{\link{cjs}} and
#' \code{\link{cjs.lnl}} to use the sin link where appropriate.  This is
#' commented out at present!  Not sure it is working correctly. \item Uses
#' optimx function for optimization which allows more methods and multiple
#' methods to be selected. \item Made fixes to \code{\link{js}} to accommodate
#' accumulation. } Version 1.0.1 (22 March 2011) \itemize{ \item Added use of
#' sparse matrices for design matrices which sped up code and reduced memory
#' consumption. The argument chunk_size was added to \code{\link{crm}},
#' \code{\link{create.dm}},\code{\link{cjs}}, and \code{\link{js}} to control
#' amount of memory usage. \item Added run timing and various print statements
#' to track progress of model.  If debug=FALSE, includes function evaluation
#' counter (every 100) and neg lnl which remains on same line unless used in
#' Rterm.  \item Note that js model needs further testing after these changes.
#' } Version 1.0.0 (Initial posting 2011) \itemize{ \item Extracted crm and
#' accompanying code from RMark and created base package which was posted on
#' bitbucket. }
#' 
#' @name Whatsnew
#' @author Jeff Laake
NULL



