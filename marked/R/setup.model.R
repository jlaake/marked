#' Defines model specific parameters (internal use)
#' 
#' Compares \code{model}, the name of the type of model (eg "CJS") to the list
#' of acceptable models to determine if it is supported and then creates some
#' global fields specific to that type of model that are used to modify the
#' operation of the code.
#' 
#' In general, the structure of the different types of models (e.g.,
#' "CJS","Recovery",...etc) are very similar with some minor exceptions.  This
#' function is not intended to be called directly by the user but it is
#' documented to enable other models to be added.  This function is called by
#' other functions to validate and setup model specific parameters.  For
#' example, for live/dead models, the length of the capture history is twice
#' the number of capture occasions and the number of time intervals equals the
#' number of capture occasions because the final interval is included with dead
#' recoveries.  Whereas, for recapture models, the length of the capture
#' history is the number of capture occasions and the number of time intervals
#' is 1 less than the number of occasions.  This function validates that the
#' model is valid and sets up some parameters specific to the model that are
#' used in the code.
#' 
#' @param model name of model type (must be in vector \code{valid.models})
#' @param nocc length of capture history string
#' @param mixtures number of mixtures
#' @return model.list - a list with following elements \item{etype}{encounter
#' type string for MARK input; typically same as model} \item{nocc}{number of
#' capture occasions} \item{num}{number of time intervals relative to number of
#' occasions (0 or -1)} \item{mixtures}{number of mixtures if any}
#' \item{derived}{logical; TRUE if model produces derived estimates}
#' @author Jeff Laake
#' @seealso \code{\link{setup.parameters}}, \code{\link{valid.parameters}}
#' @keywords utility
"setup.model" <-
function(model,nocc,mixtures=1)
{
#
# setup.model - defines list of acceptable models and creates some global fields for the model
#
# Arguments:
# 
#   model    - name of model (must be in valid.models)
#   nocc     - length of capture history string
#   mixtures - number of mixtures
#
# Value: 
#
#   model.list - a list with following elements
#                  etype - encounter type string; typically same as model name
#                  nocc  - number of capture occasions
#                  num   - number of time intervals relative to number of occasions (0 or -1)
#                  mixtures - number of mixtures if any
#                  derived - TRUE if model produces derived parameters
#
#  10 Jan 06; added known fate
#  11 Jan 06; added multistrata
#           ; robust design models
#   5 Oct 07; added Nest
#   9 Dec 07; added Occupancy and OccupHet and robust versions
#   14 Feb 07; added Jolly
#   12 Sept 08; added cjs and js R models
#   3 Mar 10; added Brownie; note ORDMS added year earlier
#   Aug 10; added CRDMS; MsLiveDead done in between
#
  valid.models=c("CJS","Recovery","Burnham","Barker","POPAN","Pradel","Pradrec","LinkBarker","Pradsen","Pradlambda",
                 "Closed","HetClosed","FullHet","Huggins","HugHet","HugFullHet","Known","Multistrata","Robust",
                 "RDHet","RDFullHet","RDHuggins","RDHHet","RDHFHet","Nest","Occupancy","OccupHet",
                 "RDOccupEG","RDOccupPE","RDOccupPG","RDOccupHetEG","RDOccupHetPE","RDOccupHetPG",
                 "OccupRPoisson","OccupRNegBin","OccupRNPoisson","OccupRNNegBin","MSOccupancy","Jolly",
                 "cjs","js","ORDMS","Brownie","MSLiveDead","CRDMS")
  stype=c(rep("mark",39),"crm","crm","mark","mark","mark","mark")
  num=c(-1,0,0,0,rep(-1,12),0,-1,rep(0,9),rep(-1,6),0,0,0,0,-1,-1,-1,-1,0,0,0,0) # of intervals relative to nocc
  divisor=c(1,2,2,2,rep(1,12),2,rep(1,16),2,2,1,1,1,1,1,1,1,2,2,1) # to compute nocc from length of ch
  default.mixtures=c(rep(1,11),2,2,1,2,2,1,1,1,2,2,1,2,2,1,1,2,1,1,1,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1)
  valid.types=c("Live","Dead","Both","Barker","POPAN","Pradel","Pradrec","LinkBarker","Pradsen","Pradlambda",
                "Closed","HetClosed","FullHet","Huggins","HugHet","HugFullHet","Known","Multistrata","Robust",
                 "RDHet","RDFullHet","RDHuggins","RDHHet","RDHFHet","Nest","Occupancy","OccupHet",
                 "RDOccupEG","RDOccupPE","RDOccupPG","RDOccupHetEG","RDOccupHetPE","RDOccupHetPG",
                 "OccupRPoisson","OccupRNegBin","OccupRNPoisson","OccupRNNegBin","MSOccupancy","Jolly",
                 "cjs","js","ORDMS","Brownie","MSLiveDead","CRDMS")
  etype= match(model,valid.models)
  derived=c(rep(FALSE,13),rep(TRUE,3),TRUE,FALSE,rep(TRUE,7),FALSE,FALSE, rep(TRUE,11),FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,TRUE)
  robust=c(rep(FALSE,18),rep(TRUE,6),rep(FALSE,3),rep(TRUE,6),rep(FALSE,8),TRUE,FALSE,FALSE,TRUE)
  closed=c(rep(FALSE,10),rep(TRUE,6),rep(FALSE,28),TRUE)
  occupancy=c(rep(FALSE,25),rep(TRUE,13),FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE)
  if(is.na(etype))
     stop("Invalid type of model = ",model," Valid types are\n", paste(valid.models,collapse=" "))
  if(mixtures==1) mixtures=default.mixtures[etype]
  return(list(etype=valid.types[etype],nocc=nocc/divisor[etype],num=num[etype],mixtures=mixtures,
               derived=derived[etype],robust=robust[etype],closed=closed[etype],occupancy=occupancy[etype],stype=stype[etype]))
}
