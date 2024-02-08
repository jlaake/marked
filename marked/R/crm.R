#' Capture-recapture model fitting function
#' 
#' Fits user specified models to some types of capture-recapture wholly in R
#' and not with MARK.  A single function that processes data, creates the
#' design data, makes the crm model and runs it
#' 
#' This function is operationally similar to the function \code{mark in RMark}
#' in that is is a shell that calls several other functions to perform the
#' following steps: 1) \code{\link{process.data}} to setup data and parameters
#' and package them into a list (processed data),2)
#' \code{\link{make.design.data}} to create the design data for each parameter
#' in the specified model, 3) \code{\link{create.dm}} to create the design
#' matrices for each parameter based on the formula provided for each
#' parameter, 4) call to the specific function for model fitting. As with \code{mark} the calling
#' arguments for \code{crm} are a compilation of the calling arguments for each
#' of the functions it calls (with some arguments renamed to avoid conflicts).
#' expects to find a value for \code{ddl}.  Likewise, if the data have not been
#' processed, then \code{ddl} should be NULL.  This dual calling structure
#' allows either a single call approach for each model or alternatively the preferred method
#' where the data area processed and the design data (\code{ddl}) created once and
#' then a whole series of models can be analyzed without repeating those steps.
#' 
#' There are some optional arguments that can be used to set initial values and
#' control other aspects of the optimization.  The optimization is done with
#' the R package/function \code{optimx} and the arguments \code{method} and
#' \code{hessian} are described with the help for that function.  In addition,
#' any arguments not matching those in the fitting functions (eg \code{cjs_admb}) are passed to
#' \code{optimx} allowing any of the other parameters to be set.  If you set
#' \code{debug=TRUE}, then at each function evaluation (\code{\link{cjs.lnl}}
#' the current values of the parameters and -2*log-likelihood value are output.
#' 
#' In the current implementation, a logit link is used to constrain the
#' parameters in the unit interval (0,1) except for probability of entry which
#' uses an mlogit and N which uses a log link. For the probitCJS model, a probit link is
#' used for the parameters. These could be generalized to
#' use other link functions. Following the notation of MARK, the parameters in
#' the link space are referred to as \code{beta} and those in the actual
#' parameter space of \code{Phi} and \code{p} as reals.
#' 
#' Initial values can be set in 2 ways.  1) Define a list of named vectors with 
#' the initial beta parameter values (eg logit link) in \code{initial}. 
#' The names of the vectors should be the parameter names in the model. Any unspecified
#' values are set to 0. 2) Specify a previously run model for initial. The code will match the names of the
#' current design matrix to the names in \code{beta} and use the appropriate
#' initial values. Any non-specified values are set to 0.  If no value is specified for initial,
#' all beta are started at a value of 0, except for the CJS model which attempts to use a glm approach to
#' setting starting values. If the glm fails then they are set to 0.
#' 
#' If you have a study with unequal time intervals between capture occasions,
#' then these can be specified with the argument \code{time.intervals}.
#' 
#' The argument \code{accumulate} defaults to \code{TRUE}.  When it is
#' \code{TRUE} it will accumulate common capture histories that also have
#' common design and common fixed values (see below) for the parameters.  This
#' will speed up the analysis because in the calculation of the likelihood
#' (\code{\link{cjs.lnl}} it loops over the unique values. In general the
#' default will be the best unless you have many capture histories and are
#' using many individual covariate(s) in the formula that would make each entry
#' unique. In that case there will be no effect of accumulation but the code
#' will still try to accumulate. In that particular case by setting
#' \code{accumulate=FALSE} you can skip the code run for accumulation.
#' 
#' Most of the arguments controlling the fitted model are contained in lists in
#' the arguments \code{model.parameters} and \code{design.parameters} which are
#' similar to their counterparts in \code{mark inb RMark}. Each is a named list
#' with the names being the parameters in the model (e.g., Phi and p in "cjs"
#' and "Phi","p","pent","N" in "js"). Each named element is also a list
#' containing various values defining the design data and model for the
#' parameter. The elements of \code{model.parameters} can include
#' \code{formula} which is an R formula to create the design matrix for the
#' parameter and \code{fixed} is a matrix of fixed values as described below.
#' The elements of \code{design.parameters} can include \code{time.varying},
#' \code{fields}, \code{time.bins},\code{age.bins}, and \code{cohort.bins}. See
#' \code{\link{create.dmdf}} for a description of the first 2 and
#' \code{\link{create.dm}} for a description of the last 3.
#' 
#' Real parameters can be set to fixed values using \code{fixed=x} where x is a
#' matrix with 3 columns and any number of rows.  The first column specifies
#' the particular animal (capture history) as the row number in the dataframe
#' \code{x}.  The second specifies the capture occasion number for the real
#' parameter to be fixed.  For \code{Phi} and \code{pent} these are 1 to
#' \code{nocc}-1 and for \code{p} they are 2 to \code{nocc} for "cjs" and 1 to
#' \code{nocc} for "js". This difference is due to the parameter labeling by
#' the beginning of the interval for Phi (e.g., survival from occasion 1 to 2)
#' and by the occasion for \code{p}.  For "cjs" \code{p} is not estimated for
#' occasion 1. The third element in the row is the real value in the closed
#' unit interval [0,1] for the fixed parameter.  This approach is completely
#' general allowing you to fix a particular real parameter for a specific
#' animal and occasion but it is a bit kludgy. Alternatively, you can set fixed
#' values by specifying values for a field called fix in the design data for a parameter.
#' If the value of fix is NA the parameter is estimated and if it is not NA then the real
#' parameter is fixed at that value.  If you also specify fixed as decribed above, they will over-ride any 
#' values you have also set with fix in the design data. To set all of the real values for a
#' particular occasion you can use the following example with the dipper data
#' as a template:
#' 
#' \code{model.parameters=list(Phi=list(formula=~1,}
#' \code{fixed=cbind(1:nrow(dipper),rep(2,nrow(dipper)),rep(1,nrow(dipper)))))}
#' 
#' The above sets \code{Phi} to 1 for the interval between occasions 2 and 3
#' for all animals. 
#' 
#' Alternatively, you could do as follows:
#' 
#' data(dipper)
#' dp=process.data(dipper)
#' ddl=make.design.data(dp)
#' ddl$Phi$fix=ifelse(ddl$Phi$time==2,1,NA)
#' 
#' At present there is no modification of the parameter count
#' to address fixing of real parameters except that if by fixing reals, a beta is not needed in the design it will be dropped.
#' For example, if you were to use ~time for Phi with survival fixed to 1 for time 2, then then beta for that time would not
#' be included.
#' 
#' To use ADMB (use.admb=TRUE), you need to install: 1) the R package R2admb, 2) ADMB, and 3) a C++ compiler (I recommend gcc compiler).
#' The following are instructions for installation with Windows. For other operating systems see (\url{http://www.admb-project.org/}) and 
#'  (\url{https://www.admb-project.org/tools/gcc/}). 
#' 
#' Windows Instructions:
#'
#'  1) In R use install.packages function or choose Packages/Install Packages from menu and select R2admb.
#' 
#'  2) Install ADMB 11: \url{https://www.admb-project.org/downloads/}. Put the software in C:/admb to
#'  avoid problems with spaces in directory name and for the function below to work.
#' 
#'  3) Install gcc compiler from: \url{https://www.admb-project.org/tools/gcc/}. Put in c:/MinGW
#' 
#' I use the following function in R to setup R2admb to access ADMB rather than adding to my path so gcc versions
#' with Rtools don't conflict. 
#' 
#' \preformatted{
#' prepare_admb=function()
#' {
#'   Sys.setenv(PATH = paste("c:/admb/bin;c:admb/utilities;c:/MinGW/bin;", 
#'         Sys.getenv("PATH"), sep = ";"))
#'     Sys.setenv(ADMB_HOME = "c:/admb")
#'     invisible()
#' }
#' }
#' To use different locations you'll need to change the values used above
#' 
#' Before running crm with use.admb=T, execute the function prepare_admb().  You could put this function or the code it 
#' contains in your .First or .Rprofile so it runs each time you start R. 
#' 
#' @param data Either the raw data which is a dataframe with at least one
#' column named ch (a character field containing the capture history) or a
#' processed dataframe
#' @param ddl Design data list which contains a list element for each parameter
#' type; if NULL it is created
#' @param begin.time Time of first capture(release) occasion
#' @param model Type of c-r model (eg, "cjs", "js") 
#' @param title Optional title; not used at present
#' @param design.parameters Specification of any grouping variables for design
#' data for each parameter
#' @param model.parameters List of model parameter specifications
#' @param initial Optional list of named vectors of initial values for beta parameters or a previously
#' run model
#' @param groups Vector of names factor variables for creating groups
#' @param time.intervals Intervals of time between the capture occasions
#' @param method optimization method 
#' @param debug if TRUE, shows optimization output
#' @param hessian if TRUE, computes v-c matrix using hessian
#' @param accumulate if TRUE, like capture-histories are accumulated to reduce
#' computation
#' @param chunk_size specifies amount of memory to use in accumulating capture
#' histories and design matrices; amount used is 8*chunk_size/1e6 MB (default
#' 80MB)
#' @param control control string for optimization functions
#' @param refit non-zero entry to refit
#' @param itnmax maximum number of iterations for optimization 
#' @param scale if TRUE, scales design matrix
#' @param run if TRUE, it runs model; otherwise if FALSE can be used to test model build components 
#' @param burnin number of iterations for mcmc burnin; specified default not realistic for actual use
#' @param iter number of iterations after burnin for mcmc (not realistic default)
#' @param use.admb if TRUE uses ADMB for cjs, mscjs or mvms models
#' @param use.tmb if TRUE runs TMB for cjs, mscjs or mvms models
#' @param crossed if TRUE it uses cjs.tpl or cjs_reml.tpl if reml=FALSE or TRUE respectively; if FALSE, then it uses cjsre which can use Gauss-Hermite integration
#' @param reml if TRUE uses restricted maximum likelihood
#' @param compile if TRUE forces re-compilation of tpl file
#' @param extra.args optional character string that is passed to admb if use.admb==TRUE
#' @param strata.labels labels for strata used in capture history; they are converted to numeric in the order listed. Only needed to specify unobserved strata. For any unobserved strata p=0..
#' @param clean if TRUE, deletes the tpl and executable files for amdb if use.admb=T
#' @param save.matrices for HMM models this option controls whether the gamma,dmat and delta matrices are saved in the model object
#' @param simplify if TRUE, design matrix is simplified to unique valus including fixed values
#' @param getreals if TRUE, compute real values and std errors for TMB models; may want to set as FALSE until model selection is complete
#' @param real.ids vector of id values for which real parameters should be output with std error information for TMB models; if NULL all ids used
#' @param check if TRUE values of gamma, dmat and delta are checked to make sure the values are valid with initial parameter values.
#' @param prior if TRUE will expect vectors of prior values in list prior.list; currently only implemented for cjsre_tmb
#' @param prior.list which contains list of prior parameters that will be model dependent
#' @param useHess if TRUE, the TMB hessian function is used for optimization; using hessian is typically slower with many parameters but can result in a better solution
#' @param optimize if TRUE, optimizes to get parameter estimates; set to FALSE to extract estimates of ADREPORTed values only
#' @param unit_scale default TRUE, if FALSE any time scaled parameter (e.g. Phi,S) is scaled when computing real value such that it represents the length of the interval rather than a unit interval
#' @param ... optional arguments passed to js or cjs and optimx
#' @importFrom graphics boxplot par
#' @importFrom stats as.formula binomial coef density
#'             glm.fit median model.frame model.matrix optim
#'              plogis pnorm predict rgamma rmultinom
#'              rnorm sd nlminb
#' @importFrom utils capture.output flush.console
#'             read.delim
#' @importFrom methods as is
#' @import data.table
#' @import bookdown
#' @import kableExtra
#' @import knitr
#' @return crm model object with class=("crm",submodel), eg "CJS".
#' @author Jeff Laake
#' @export crm
#' @import optimx Matrix Rcpp numDeriv
#' @useDynLib marked,  .registration = TRUE
#' @seealso \code{\link{cjs_admb}}, \code{\link{js}},
#' \code{\link{make.design.data}},\code{\link{process.data}}
#' @keywords models
#' @examples
#' {
#' # cormack-jolly-seber model
#' # fit cjs models with crm
#' data(dipper)
#' dipper.proc=process.data(dipper,model="cjs",begin.time=1)
#' dipper.ddl=make.design.data(dipper.proc)
#' mod.Phit.pt=crm(dipper.proc,dipper.ddl,
#'    model.parameters=list(Phi=list(formula=~time),p=list(formula=~time)))
#' mod.Phit.pt
#' mod.Phisex.pdot=crm(dipper.proc,dipper.ddl,
#'    model.parameters=list(Phi=list(formula=~sex),p=list(formula=~1)))
#' mod.Phisex.pdot
#' # demo initial value setting
#' mod.Phisex.ptime=crm(dipper.proc,dipper.ddl,
#'    model.parameters=list(Phi=list(formula=~sex),p=list(formula=~time)),initial=mod.Phit.pt)
#' mod.Phisex.ptime
#' mod.Phisex.ptime=crm(dipper.proc,dipper.ddl,
#'    model.parameters=list(Phi=list(formula=~sex),p=list(formula=~time)),initial=list(Phi=0,p=0))
#' mod.Phisex.ptime
#' 
#' ## if you have RMark installed you can use this code to run the same models 
#' ## by removing the comment symbol
#' #library(RMark)
#' #data(dipper)
#' #mod0=mark(dipper,
#' #model.parameters=list(Phi=list(formula=~time),p=list(formula=~time)),output=FALSE)
#' #summary(mod0,brief=TRUE)
#' #mod1=mark(dipper,
#' #model.parameters=list(Phi=list(formula=~1),p=list(formula=~1)),output=FALSE)
#' #summary(mod1,brief=TRUE)
#' #mod2<-mark(dipper,groups="sex",
#' #model.parameters=list(Phi=list(formula=~sex),p=list(formula=~1)),output=FALSE)
#' #summary(mod2,brief=TRUE)
#' # jolly seber model
#' crm(dipper,model="js",groups="sex",
#'    model.parameters=list(pent=list(formula=~sex),N=list(formula=~sex)),accumulate=FALSE)
#' # examples showing use of unit.scale
#' dipper.proc=process.data(dipper,model="cjs",begin.time=1,time.intervals=c(.1,.2,.3,.4,.5,.6))
#' dipper.ddl=make.design.data(dipper.proc)
#' mod.Phit.p=crm(dipper.proc,dipper.ddl,
#'                model.parameters=list(Phi=list(formula=~time),p=list(formula=~1)),
#'                hessian=TRUE,unit_scale=TRUE)
#' mod.Phit.p
#' mod.Phit.p$results$reals
#' 
#' dipper.proc=process.data(dipper,model="cjs",begin.time=1,time.intervals=c(.1,.2,.3,.4,.5,.6))
#' dipper.ddl=make.design.data(dipper.proc)
#' mod.Phit.p=crm(dipper.proc,dipper.ddl,
#'                model.parameters=list(Phi=list(formula=~time),p=list(formula=~1)),
#'                hessian=TRUE,unit_scale=FALSE)
#' mod.Phit.p
#' mod.Phit.p$results$reals
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' # if you have RMark installed you can use this code to run the same models 
#' # by removing the comment 
#' #data(dipper)
#' #data(mstrata)
#' #mark(dipper,model.parameters=list(p=list(formula=~time)),output=FALSE)$results$beta
#' #mark(mstrata,model="Multistrata",model.parameters=list(p=list(formula=~1),
#' # S=list(formula=~1),Psi=list(formula=~-1+stratum:tostratum)),
#' # output=FALSE)$results$beta
#' #mod=mark(dipper,model="POPAN",groups="sex",
#' #   model.parameters=list(pent=list(formula=~sex),N=list(formula=~sex)))
#' #summary(mod)
#' #CJS example with hmm
#' crm(dipper,model="hmmCJS",model.parameters = list(p = list(formula = ~time)))
#' ##MSCJS example with hmm
#' data(mstrata)
#' ms=process.data(mstrata,model="hmmMSCJS",strata.labels=c("A","B","C"))
#' ms.ddl=make.design.data(ms)
#' ms.ddl$Psi$fix=NA
#' ms.ddl$Psi$fix[ms.ddl$Psi$stratum==ms.ddl$Psi$tostratum]=1
#' crm(ms,ms.ddl,model.parameters=list(Psi=list(formula=~-1+stratum:tostratum)))
#' }
#' }
crm <- function(data,ddl=NULL,begin.time=1,model="CJS",title="",model.parameters=list(),design.parameters=list(),initial=NULL,
 groups = NULL, time.intervals = NULL,debug=FALSE, method=NULL, hessian=FALSE, accumulate=TRUE,chunk_size=1e7, 
 control=list(),refit=1,itnmax=5000,scale=TRUE,run=TRUE,burnin=100,iter=1000,use.admb=FALSE,use.tmb=FALSE,crossed=NULL,reml=FALSE,compile=FALSE,extra.args=NULL,
 strata.labels=NULL,clean=NULL,save.matrices=FALSE,simplify=FALSE,getreals=FALSE,real.ids=NULL,check=FALSE,prior=FALSE,prior.list=NULL,useHess=FALSE,optimize=TRUE,unit_scale=TRUE,...)
{
model=toupper(model)
ptm=proc.time()
if(is.null(crossed))crossed=FALSE
if(crossed)accumulate=FALSE
#
#  If the data haven't been processed (data$data is NULL) do it now with specified or default arguments
# 
if(is.null(data$data))
{
   if(!is.null(ddl))
   {
      warning("Warning: specification of ddl ignored, as data have not been processed")
      ddl=NULL
   }
   if(debug)message("Model: ",model,"\n")
   if(debug)message("Processing data...\n")
   flush.console()
   data.proc=process.data(data,begin.time=begin.time, model=model,mixtures=1, 
	   groups = groups, age.var = NULL, initial.ages = NULL, 
	   time.intervals = time.intervals,nocc=NULL,accumulate=accumulate,strata.labels=strata.labels)
}   
else
{
	data.proc=data
	model=data$model
}
if(is.null(method))
{
  if(substr(model,1,4)%in%c("MVMS","MSLD"))
    method="nlminb"
  else
    method="BFGS"
}
if(model%in%c("MSLD"))
  use.tmb=TRUE
#
# Setup parameter list
#
number.of.groups=1
if(!is.null(data.proc$group.covariates))number.of.groups=nrow(data.proc$group.covariates)
par.list=setup.parameters(data.proc$model,check=TRUE)
#
# Check validity of parameter list; stop if not valid
#
if(!valid.parameters(model,model.parameters)) stop()
parameters=setup.parameters(data.proc$model,model.parameters,data$nocc,number.of.groups=number.of.groups)
parameters=parameters[par.list]
# See if any formula contain random effects and set re
re=FALSE
for (i in 1:length(parameters))
{
	if(is.null(parameters[[i]]$formula)) parameters[[i]]$formula=~1
	mlist=proc.form(parameters[[i]]$formula)
	
	if(!is.null(mlist$re.model))
	{
		re_names=sub("^\\s+", "",sapply(strsplit(names(mlist$re.model),"\\|"),function(x)x[2]))
		if(length(re_names)>1 | !"id" %in% re_names) crossed=TRUE  
		if((length(re_names)> 1 || re_names[1]!="time" ||  use.admb ) & any(data.proc$freq>1)) 
			stop("\n data cannot be accumulated (freq>1) except with temporal random effects only; set accumulate=FALSE\n")	
#		else
#            if(use.tmb & any(data.proc$freq>1))
#			{
#				re_names=re_names[re_names!="time"]
#				if(!all(re_names %in% names(data.proc$data)))
#				{
#					message("\n data cannot be accumulated (freq>1) unless the following fields are in the data\n")
#					message(paste(re_names),"\n")
#					stop()
#				}
#			}
		re=TRUE
	}
	if(parameters[[i]]$nointercept)parameters[[i]]$remove.intercept=TRUE
}
# currently if re, then set use.admb to TRUE unless use.tmb=T
if(re&!use.tmb) {
	use.admb=TRUE
	if(is.null(clean))clean=TRUE
}
if(use.admb | (!use.tmb &toupper(model)%in%c("MSCJS")))
{
	if(!re) crossed=FALSE
	if(is.null(clean))clean=TRUE
}
if((use.tmb|toupper(model)%in%c("MSLD"))&is.null(clean))clean=FALSE
#
# If the design data have not been constructed, do so now
#
external.ddl=FALSE
if(is.null(ddl)) 
{
  if(debug)message("Creating design data...\n")
	flush.console()
	ddl=make.design.data(data.proc,design.parameters)
} else
{
  if(is.character(ddl)&&toupper(ddl)=="EXTERNAL")
  {
    external.ddl=TRUE
    if(file.exists("ddl.rda"))
    {
      load(file="ddl.rda")
      if(!exists("ddl")) stop("\nexternal ddl.rda file must contain object named ddl")
    } else
    {
      stop("\nCannot find external file named ddl.rda")
    }
  }
	for (i in 1:length(parameters))
	{
		if(!is.null(ddl[[i]]$order))
		   if(any(ddl[[i]]$order!=1:nrow(ddl[[i]]))) 
			   stop(paste("Design data for parameter",names(parameters)[i],"is out of order."))
	}
    if(!is.null(design.parameters))
		for(parname in names(design.parameters))
			ddl$design.parameters[[parname]]=c(ddl$design.parameters[[parname]],design.parameters[[parname]])
	design.parameters=ddl$design.parameters
}
ddl=set.fixed(ddl,parameters) #   setup fixed values if old way used
if(substr(model,1,4)=="MVMS"&check)
{
  check_mlogit=function(x){ifelse(sum(x,na.rm=TRUE)==0||(any(is.na(x))&!(any(x[!is.na(x)]==1))),TRUE,FALSE)}
	if(is.null(ddl$pi$fix))
		message("\n Warning: No values provided for fix for pi. Must have a reference cell via formula.")
	else
	{
		bad_pi=sapply(split(ddl$pi$fix,ddl$pi$id),check_mlogit)
		if(any(bad_pi))
		  {
		    message("\n Check values of fix for pi with id:")
		    cat(names(which(bad_pi)))
		  }
	}
	if(is.null(ddl$delta$fix))
	{
		message("\n Warning: No values provided for fix for delta. Must have a reference cell via formula.")
	}else
	{
	  xx=list(id=ddl$delta$id,occ=ddl$delta$occ,stratum=ddl$delta$stratum)
		bad_delta=sapply(split(ddl$delta$fix,xx),check_mlogit)
	  if(any(bad_delta))
	  {
	     message("\n Warning: Check values of fix for delta for the following records with id.occ.stratum.")
	     cat(names(which(bad_delta)))  
	  }
	}
	if(is.null(ddl$Psi$fix))
	{
		message("\n Warning: No values provided for fix for Psi. Must have a reference cell via formula.")
	}else
	{
	  xx=list(id=ddl$Psi$id,occ=ddl$Psi$occ,stratum=ddl$Psi$stratum)
	  bad_Psi=sapply(split(ddl$Psi$fix,xx),check_mlogit)
	  if(any(bad_Psi))
	  {
	    message("\n Warning: Check values of fix for Psi for the following records with id.occ.stratum.")
	    cat(names(which(bad_Psi)))  
	  }
	}
}
if(model=="MSCJS" | model%in%c("MSLD") | (substr(model,1,4)=="MVMS" & (use.admb | use.tmb))) 
{
  fullddl=ddl
  if(debug)message("Simplifying design data\n")
	ddl=simplify_ddl(ddl,parameters) # add indices to ddl and reduce ddl to unique values used
}
else 
  fullddl=NULL
if(simplify)
{
	simplify=FALSE
	message("\nsimplify argument has been disabled")
}
# check to see if all values for a parameter have been fixed.  If so, then set formula to ~0
for (i in 1:length(parameters))
{
	if(!is.null(ddl[[i]]$fix))
	{
		if(!is.null(ddl[[i]]$fix) && all(!is.na(ddl[[i]]$fix)))
		{
			message(paste("All values for",names(parameters)[i],"have been fixed. Setting formula to ~0\n"))
			parameters[[i]]$formula=~0
		} else {
			if(parameters[[i]]$formula==~0)
				stop(paste("Cannot use formula ~0 for",names(parameters)[i],"when some of the parameters must be estimated.\n"))
		}
	} else
	   if(parameters[[i]]$formula==~0)
		   stop(paste("Cannot use formula ~0 for",names(parameters)[i],"when some of the parameters must be estimated.\n"))
}
# Create design matrices for each parameter
if(debug)message("Creating design matrices\n")
dml=create.dml(ddl,model.parameters=parameters,design.parameters=design.parameters,chunk_size=chunk_size,simplify=simplify,use.admb=use.admb)
if(model=="MSCJS"| model%in%c("MSLD")|(substr(model,1,4)=="MVMS" & (use.admb | use.tmb)))
  {
  fulldml=dml
  for(parx in names(dml))
  {
    fulldml[[parx]]$fe=dml[[parx]]$fe[ddl[[paste(parx,".indices",sep="")]],,drop=FALSE]
    parameters[[parx]]$indices=ddl[[paste(parx,".indices",sep="")]]    
  }
} 
# For HMM call set.initial to get ptype and set initial values
if(substr(model,1,3)=="HMM"|(nchar(model)>=4 &substr(model,1,4)=="MVMS"))
	initial.list=set.initial(names(dml),dml,initial)
else
	initial.list=NULL
# if not running, return object with data,ddl,dml etc
if(!run) return(list(model=model,data=data.proc,model.parameters=parameters,design.parameters=design.parameters,ddl=ddl,dml=dml,results=initial.list))
# Depending on method set some values
if("SANN"%in%method)
{
	if(length(method)>1)
		warning("***SANN can only be used by itself; other methods ignored.")
	method="SANN"
    control$maxit=itnmax
}
if("nlminb"%in%method)
{
	control$eval.max=itnmax
	control$iter.max=itnmax
}
# Call estimation function which depends on the model
if(debug)message("Fitting model\n")
if(model=="CJS")
	if(use.tmb)
	{
		runmodel=cjs_tmb(data.proc,ddl,dml,parameters=parameters,initial=initial,method=method,hessian=hessian,debug=debug,accumulate=accumulate,chunk_size=chunk_size,
				refit=refit,control=control,itnmax=itnmax,scale=scale,crossed=crossed,compile=compile,extra.args=extra.args,reml=reml,clean=clean,getreals=getreals,
				prior=prior,prior.list=prior.list,useHess=useHess,...)
	} else
		runmodel=cjs_admb(data.proc,ddl,dml,parameters=parameters,initial=initial,method=method,hessian=hessian,debug=debug,accumulate=accumulate,chunk_size=chunk_size,
		          refit=refit,control=control,itnmax=itnmax,scale=scale,use.admb=use.admb,crossed=crossed,compile=compile,extra.args=extra.args,reml=reml,clean=clean,...)
if(model=="JS")
    runmodel=js(data.proc,ddl,dml,parameters=parameters,initial=initial,method=method,hessian=hessian,debug=debug,accumulate=FALSE,chunk_size=chunk_size,
		          refit=refit,control=control,itnmax=itnmax,scaleDM=scale,...)
if(model=="MSCJS")
	if(use.tmb)
	{
		runmodel=mscjs_tmb(data.proc,ddl,fullddl,dml,parameters=parameters,initial=initial,method=method,hessian=hessian,debug=debug,accumulate=accumulate,chunk_size=chunk_size,
				refit=refit,control=control,itnmax=itnmax,scale=scale,re=re,compile=compile,extra.args=extra.args,clean=clean,getreals=getreals,useHess=useHess,...)
  }else{
	    runmodel=mscjs(data.proc,ddl,dml,parameters=parameters,initial=initial,method=method,hessian=hessian,debug=debug,accumulate=accumulate,chunk_size=chunk_size,
				   refit=refit,control=control,itnmax=itnmax,scale=scale,re=re,compile=compile,extra.args=extra.args,clean=clean,...)
  }
if(model%in%c("MSLD"))
{
  save(fullddl,file="tmp.rda")
  rm(fullddl)
  runmodel=msld_tmb(data.proc,ddl,dml,parameters=parameters,initial=initial,method=method,hessian=hessian,debug=debug,accumulate=accumulate,chunk_size=chunk_size,
                    refit=refit,control=control,itnmax=itnmax,scale=scale,re=re,compile=compile,extra.args=extra.args,
                    clean=clean,getreals=getreals,useHess=useHess,savef=save.matrices,...)
  if(save.matrices){
    runmodel$mat=runmodel$f$report()
    names(runmodel$mat)=c("dmat","gamma") 
    runmodel$mat$delta=matrix(0,nrow=nrow(data.proc$data),ncol=length(data.proc$strata.labels)*2+1)
    for(i in 1:nrow(data.proc$data))
        runmodel$mat$delta[i,data.proc$start[i,1]]=1
    # Adjust dmat to have nocc+1 entries for time to work with global decode
    temp=array(0,dim=c(dim(runmodel$mat$dmat)[1],dim(runmodel$mat$dmat)[2]+1,dim(runmodel$mat$dmat)[3],dim(runmodel$mat$dmat)[4]))
    temp[,2:dim(temp)[2],,]=runmodel$mat$dmat
    for(i in 1:nrow(temp))
      for(j in 1:length(data.proc$strata.labels))
      {
        temp[i,data.proc$start[i,2],j+1,j]=1
        temp[i,data.proc$start[i,2],1,j]=0
      }
    runmodel$mat$dmat=temp
  }
  load("tmp.rda")
}
if(model=="PROBITCJS")
{
	if(is.null(initial))
	{
	    imat=process.ch(data.proc$data$ch,data.proc$data$freq,all=FALSE)
	    runmodel=probitCJS(ddl,dml,parameters=parameters,design.parameters=design.parameters,
		               imat=imat,iter=iter,burnin=burnin)
    }else
	    runmodel=probitCJS(ddl,dml,parameters=parameters,design.parameters=design.parameters,
					   initial=initial,iter=iter,burnin=burnin)	   
}
if(substr(model,1,3)=="HMM"|(nchar(model)>=4 &substr(model,1,4)=="MVMS"))
{
	if(substr(model,1,4)=="MVMS")
		sup=data.proc$fct_sup(list(obslevels=data.proc$ObsLevels))
	else
		sup=NULL
	if(is.null(data.proc$strata.list) | substr(model,1,4)=="MVMS"){
		mx=data.proc$m
	}else{
		mx=list(ns=length(data.proc$strata.list$states),na=length(data.proc$strata.list[[names(data.proc$strata.list)[names(data.proc$strata.list)!="states"]]]))
	}
	if((use.admb | use.tmb) & model=="MVMSCJS") 
	{
		# call HMMlikelihood to check for problems in setup
	  if(check)
	  {
	    message("Checking for problems in design data setup\n")
	    xx=HMMLikelihood(par=unlist(initial.list$par),xx=data.proc$ehmat,mx=mx,
	                     type=initial.list$ptype,T=data.proc$nocc,xstart=data.proc$start,freq=data.proc$freq,
	                     fct_dmat=data.proc$fct_dmat,fct_gamma=data.proc$fct_gamma,fct_delta=data.proc$fct_delta,ddl=fullddl,dml=fulldml,
	                     parameters=parameters,sup=sup,check=TRUE)
	  }
		# call mvmscjs to run admb program
		if(use.admb)
		  runmodel=mvmscjs(data.proc,ddl,fullddl,dml,parameters=parameters,initial=initial,method=method,hessian=hessian,debug=debug,accumulate=accumulate,chunk_size=chunk_size,
				refit=refit,control=control,itnmax=itnmax,scale=scale,re=re,compile=compile,extra.args=extra.args,clean=clean,sup=sup,...)
		else
		  runmodel=mvmscjs_tmb(data.proc,ddl,fullddl,dml,parameters=parameters,initial=initial,method=method,hessian=hessian,debug=debug,accumulate=accumulate,chunk_size=chunk_size,
		                   refit=refit,control=control,itnmax=itnmax,re=FALSE,compile=compile,clean=clean,sup=sup,getreals=getreals,real.ids=real.ids,useHess=useHess,optimize=optimize,...)
		par=coef(runmodel)[,1]
		runmodel$options=c(runmodel$options,list(accumulate=accumulate,initial=initial.list$par,method=method,
				chunk_size=chunk_size,itnmax=itnmax,control=control,use.tmb=use.tmb,use.admb=use.admb))
		
	} else
	{
		xx=HMMLikelihood(par=unlist(initial.list$par),xx=data.proc$ehmat,mx=mx,
				type=initial.list$ptype,T=data.proc$nocc,xstart=data.proc$start,freq=data.proc$freq,
				fct_dmat=data.proc$fct_dmat,fct_gamma=data.proc$fct_gamma,fct_delta=data.proc$fct_delta,ddl=ddl,dml=dml,
				parameters=parameters,sup=sup,check=TRUE)
		runmodel=optimx(unlist(initial.list$par),HMMLikelihood,method=method,debug=debug,hessian=hessian,itnmax=itnmax,xx=data.proc$ehmat,mx=mx,
				type=initial.list$ptype,T=data.proc$nocc,xstart=data.proc$start,freq=data.proc$freq,control=control,
				fct_dmat=data.proc$fct_dmat,fct_gamma=data.proc$fct_gamma,fct_delta=data.proc$fct_delta,ddl=ddl,dml=dml,
				parameters=parameters,sup=sup,check=FALSE)
		par=coef(runmodel, order="value")[1, ]
		runmodel=list(optim.details=as.list(summary(runmodel, order="value",par.select=FALSE)[1, ]))
		if(hessian)runmodel$hessian=attr(runmodel$optim.details,"details")$nhatend
		runmodel$convergence=runmodel$optim.details$convcode
		runmodel$options=list(accumulate=accumulate,initial=initial.list$par,method=method,
				chunk_size=chunk_size,itnmax=itnmax,control=control,use.tmb=use.tmb,use.admb=use.admb)
		runmodel$model_data=list(Phi.dm=dml$Phi$fe,p.dm=dml$p$fe,delta.dm=dml$delta$fe,Psi.dm=dml$Psi$fe,pi.dm=dml$pi$fe)
		#S.fixed,r.fixed,p.fixed,Psi.fixed and time.intervals
	}
		if(save.matrices)
		{
		  if(!is.null(fullddl))
			   runmodel$mat=HMMLikelihood(par=par,type=initial.list$ptype,xx=data.proc$ehmat,mx=mx,T=data.proc$nocc,xstart=data.proc$start,freq=data.proc$freq,
				   	fct_dmat=data.proc$fct_dmat,fct_gamma=data.proc$fct_gamma,fct_delta=data.proc$fct_delta,ddl=fullddl,dml=fulldml,parameters=parameters,return.mat=TRUE,sup=sup)
		  else
		    runmodel$mat=HMMLikelihood(par=par,type=initial.list$ptype,xx=data.proc$ehmat,mx=mx,T=data.proc$nocc,xstart=data.proc$start,freq=data.proc$freq,
		                               fct_dmat=data.proc$fct_dmat,fct_gamma=data.proc$fct_gamma,fct_delta=data.proc$fct_delta,ddl=ddl,dml=dml,parameters=parameters,return.mat=TRUE,sup=sup)
		  if(model=="HMMCJS")
			{
				dimnames(runmodel$mat$gamma)[3:4]=list(c("Alive","Dead"),c("Alive","Dead"))
				dimnames(runmodel$mat$dmat)[3:4]=list(c("Missed","Seen"),c("Alive","Dead"))
			}else
			{
				dimnames(runmodel$mat$gamma)[3:4]=list(c(data.proc$strata.labels,"Dead"),c(data.proc$strata.labels,"Dead"))
				dimnames(runmodel$mat$dmat)[3:4]=list(data.proc$ObsLevels,c(data.proc$strata.labels,"Dead"))
			}
			names(dimnames(runmodel$mat$gamma))=c("Id","Occasion","From_state","To_state")
			names(dimnames(runmodel$mat$dmat))=c("Id","Occasion","Observation","State")
		}
		parlist=split(par,initial.list$ptype)
		par=vector("list",length=length(names(initial.list$par)))
		names(par)=names(initial.list$par)
		for(p in names(parlist))
		{
			par[[p]]=parlist[[p]]
			names(par[[p]])=colnames(dml[[p]]$fe)	
		}
		runmodel$beta=par
		runmodel$par=NULL
		if(is.null(runmodel$neg2lnl)) 
			runmodel$neg2lnl=2*runmodel$optim.details$value
		runmodel$AIC=runmodel$neg2lnl+2*sum(sapply(runmodel$beta,length))
		if(!is.null(runmodel$hessian))
		{
			runmodel$beta.vcv=solvecov(runmodel$hessian)$inv
			colnames(runmodel$beta.vcv)=names(unlist(runmodel$beta))
			rownames(runmodel$beta.vcv)=colnames(runmodel$beta.vcv)
		}
		class(runmodel)=c("crm","mle",model)
}
#
# Return fitted MARK model object or if external, return character string with same class and save file
if(!is.null(runmodel$convergence) && runmodel$convergence!=0&!use.admb)
{
	warning("******Model did not converge******")
	msg=attr(runmodel$optim.details,"details")$message
	if(is.null(msg)) msg="Exceeded maximum number of iterations"
	warning(msg)
}

object=list(model=model,data=data.proc,model.parameters=parameters,design.parameters=design.parameters,results=runmodel)
class(object)=class(runmodel)
if(!use.tmb&!re & !model%in%c("MSCJS","MSLD") & (nchar(model)<4 | (nchar(model)>=4 & substr(model,1,4)!="MVMS")))
   object$results$reals=predict(object,ddl=ddl,unique=TRUE,se=hessian,unit_scale=unit_scale)
#if(use.tmb & (nchar(model)>=4 & substr(model,1,4)=="MVMS") & getreals)
#  object$results$reals=predict(object,ddl=ddl,real.ids=real.ids,se=hessian)
if(file.exists("tmp.rda"))unlink("tmp.rda")
message(paste("\nElapsed time in minutes: ",round((proc.time()[3]-ptm[3])/60,digits=4),"\n"))
return(object)
}
# solvecov code was taken from package fpc: Christian
# Hennig <chrish@@stats.ucl.ac.uk> \url{http://www.homepages.ucl.ac.uk/~ucakche/}
solvecov=function (m, cmax = 1e+10)
# from package fpc
{
	options(show.error.messages = FALSE)
	covinv <- try(solve(m))
	if (!is(covinv,"try-error"))
		coll = FALSE
	else {
		p <- nrow(m)
		cove <- eigen(m, symmetric = TRUE)
		coll <- TRUE
		if (min(cove$values) < 1/cmax) {
			covewi <- diag(p)
			for (i in 1:p) if (cove$values[i] < 1/cmax)
					covewi[i, i] <- cmax
				else covewi[i, i] <- 1/cove$values[i]
		}
		else covewi <- diag(1/cove$values, nrow = length(cove$values))
		covinv <- cove$vectors %*% covewi %*% t(cove$vectors)
	}
	options(show.error.messages = TRUE)
	out <- list(inv = covinv, coll = coll)
	out
}

