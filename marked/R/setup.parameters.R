#' Setup parameter structure specific to model (internal use)
#' 
#' Defines list of parameters used in the specified type of model
#' (\code{model}) and adds default values for each parameter to the list of
#' user specified values (eg formula, link etc). \code{parameters} is a list of lists.
#' 
#' @param model type of model (e.g., "cjs")
#' @param parameters list of model parameter specifications
#' @param nocc number of occasions (value only specified if needed)
#' @param check if TRUE only the vector of parameter names is returned
#' \code{par.list}
#' @param number.of.groups number of groups defined for data
#' @return The return value depends on the argument \code{check}. If it is TRUE
#' then the return value is a vector of the names of the parameters used in the
#' specified type of model. For example, if \code{model="cjs"} then the return
#' value is \code{c("Phi","p")}.  This is used by the function
#' \code{\link{valid.parameters}} to make sure that parameter specifications
#' are valid for the model (i.e., specifying recovery rate r for "CJS" would
#' give an error).  If the function is called with the default of
#' \code{check=FALSE}, the function returns a list of parameter specifications
#' which is a modification of the argument \code{parameters} which adds
#' parameters not specified and default values for all types of parameters that
#' were not specified. The list length and names of the list elements depends
#' on the type of model. Each element of the list is itself a list with varying
#' numbers of elements which depend on the type of parameter although some
#' elements are the same for all parameters.  Below the return value list is
#' shown generically with parameters named p1,...,pk.  \tabular{ll}{ \code{p1}
#' \tab List of specifications for parameter 1 \cr \code{p2} \tab List of
#' specifications for parameter 2 \cr . \tab \cr .  \tab \cr .  \tab \cr
#' \code{pk} \tab List of specifications for parameter k \cr }
#' 
#' The elements for each parameter list all include: \tabular{ll}{ \code{begin}
#' \tab 0 or 1; beginning time for the first parameter relative to first
#' occasion \cr \code{num} \tab 0 or -1; number of parameters relative to
#' number of occassions \cr \code{type} \tab type of PIM structure; either
#' "Triang" or "Square" \cr \code{formula} \tab formula for parameter model
#' (e.g., \code{~time}) \cr \code{link} \tab link function for parameter (e.g.,
#' \code{"logit"}) \cr }
#' 
#' and may include: \tabular{ll}{ \code{share} \tab only valid for p in closed
#' capture models; if TRUE p and c models shared \cr \code{mix} \tab only valid
#' for closed capture heterogeneity models; if TRUE mixtures are used \cr
#' \code{rows} \tab only valid for closed capture heterogeneity models \cr
#' \code{fixed} \tab fixed values specified by user and not used modified in
#' this function \cr }
#' @author Jeff Laake
#' @seealso \code{\link{setup.model}},\code{\link{valid.parameters}}
#' @keywords utility
"setup.parameters" <-
function(model,parameters=list(),nocc=NULL,check=FALSE,number.of.groups=1)
# ----------------------------------------------------------------------------------------
#  setup.parameters  - fills in value for begin and num for each parameter type depending
#                      on the type of c-r model. num defines number of parameters relative to
#                      number of occasions.  begin defines the first occasion which is relevant
#                      to the parameter
#
#  Arguments:
#    model      - type of model ("CJS", "Burnham" etc)
#    parameters - list of model parameter specifications
#    nocc       - number of occasions (value only specified if needed)
#    check      - default is FALSE; if TRUE it only returns list of parameter names
#
#  Value:
#    parameters - updated list of model parameter specifications with new fields added
#
#
#  10 Jan 06; added pim.type to list for triang structured parameters
# ----------------------------------------------------------------------------------------
{
#
#  Create valid parameter list depending on model.
#
   if(model=="probitCJS" | model=="cjs")par.list=c("Phi","p")
   if(model=="js") par.list=c("Phi","p","pent","N")
#
#  If this is just a parameter check, return par.list
#
   if(check)return(par.list)
#
#  For each parameter create an empty list if none specified in input
#   
   if(length(parameters)>0)
   {
      for (i in 1:length(par.list))
      {
         if(is.na(names(parameters[par.list[i]])))
            parameters[[par.list[i]]]=list()
      }
   }
#
#  Next depending on model type, assign non-specified default values
#
   if(model=="probitCJS" | model=="cjs")
   {
      parameters$Phi$num=-1
      parameters$Phi$begin=0
      if(is.null(parameters$Phi$default))parameters$Phi$default=1
      parameters$Phi$type="Triang"
      if(is.null(parameters$Phi$pim.type))parameters$Phi$pim.type="all"
      if(is.null(parameters$Phi$link))parameters$Phi$link="logit"
      parameters$p$num=-1
      parameters$p$begin=1
      parameters$p$type="Triang"
      if(is.null(parameters$p$default))parameters$p$default=0
      if(is.null(parameters$p$pim.type))parameters$p$pim.type="all"
      if(is.null(parameters$p$link))parameters$p$link="logit"
   }
   else
   if(model=="js")
   {
      parameters$Phi$begin=0
      parameters$Phi$num=-1
      parameters$Phi$type="Square"
      if(is.null(parameters$Phi$default))parameters$Phi$default=1
      if(is.null(parameters$Phi$link))parameters$Phi$link="logit"
      parameters$p$num=0
      parameters$p$begin=0
      parameters$p$type="Square"
      if(is.null(parameters$p$default))parameters$p$default=0
      if(is.null(parameters$p$link))parameters$p$link="logit"
      parameters$pent$num=-1
      parameters$pent$begin=1
      parameters$pent$type="Square"
      if(is.null(parameters$pent$default))parameters$pent$default=0
      if(number.of.groups>1)
         if(is.null(parameters$pent$formula))parameters$pent$formula=~group
      else
         if(is.null(parameters$pent$formula))parameters$pent$formula=~1
      if(is.null(parameters$pent$link))parameters$pent$link="mlogit"
      parameters$N$num=-(nocc-1)
      parameters$N$begin=0
      parameters$N$type="Square"
      parameters$N$leave.unused=TRUE
      if(number.of.groups>1)
         if(is.null(parameters$N$formula))parameters$N$formula=~group
      else
         if(is.null(parameters$N$formula))parameters$N$formula=~1
      if(is.null(parameters$N$link))parameters$N$link="log"
   } 
   return(parameters)
}
