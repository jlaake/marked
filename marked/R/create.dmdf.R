#' Creates a dataframe with all the design data for a particular parameter in a
#' crm model
#' 
#' Creates a dataframe with all the design data for a particular parameter in a
#' crm model which currently includes "cjs" or "js".  These design data are
#' fundamentally different than the design data created for \code{mark} models
#' as explained below.
#' 
#' This function is intended to be called from \code{\link{make.design.data}}.
#' It takes the data in \code{x} and creates a dataframe with all of the data
#' needed to fit capture-recapture models (crm) which currently includes "cjs"
#' (Cormack-Jolly-Seber) or "js" (POPAN formulation of the Jolly-Seber model).
#' Before jumping into the details it is useful to have an understanding of the
#' differences between MARK (via the \code{mark in RMark} function) and the
#' package \code{mra} written by Trent McDonald and how they relate to the
#' implementation in \code{\link{cjs}}.  With MARK, animals can be placed in
#' groups and the parameters for the model specified via PIMs (parameter
#' information matrices) link the parameters to the specific animals.  For
#' example, if for a particular group the \code{Phi} PIM is
#' 
#' \preformatted{ 
#'   1 2 3 
#'     4 5 
#'       6 
#' }
#' 
#' Then animals in that group that were first caught/released on occasion 1
#' have the parameters 1,2,3 for the 3 occasions.  Those first caught/released
#' on occasion 2 have parameters 4 and 5 for occasions 2 and 3 and those first
#' caught/released on occasion 3 have parameter 6 for occasion 3.  Another
#' group of animals would have a different set of indices for the same set of
#' parameters so they could be discriminated.  Many people find this rather
#' confusing until they get used to it but even then if you have many different
#' groups and many occasions, then the indexing process is prone to error.
#' Thus, the rationale for RMark which automates the PIM construction and its
#' use is largely transparent to the user. What RMark does is to create a data
#' structure called design data that automatically assigns \code{design data}
#' to the indices in the PIM.  For example, 1 to 6 would be given the data used
#' to create that group and 1 to 3 would be assigned to cohort 1, ... and 1
#' would be assigned to occasion 1 and 2 and 4 would be assigned to occasion 2
#' etc.  It also creates an age field which follows the diagonals and can be
#' initialized with the intial age at the time of first capture which is group
#' specific.  With a formula and these design data a design matrix can be
#' constructed for the model where the row in the design matrix is the
#' parameter index (e.g., the first 6 rows would be for parameters 1 to 6 as
#' shown above).  That would be all good except for individual covariates which
#' are not group-specific.  MARK handles individual covariates by specifying
#' the covariate name (eg "weight") as a string in the design matrix.  Then for
#' each capture history in plugs in the actual covariate values for that animal
#' to complete the design matrix for that animal. For more details see Laake
#' and Rexstad (2008).
#' 
#' From a brief look at package \code{mra} and personal communication with
#' Trent McDonald, I give the following brief and possibly incorrect
#' description of the pacakge \code{mra} at the time of writing (28 Aug 2008).
#' In that package, the whole concept of PIMS is abandoned and instead
#' covariates are constructed for each occasion for each animal.  Thus, each
#' animal is effectively put in its own group and it has a parameter for each
#' occasion.  This novel approach is quite effective at blurring the lines
#' between design data and individual covariate data and it removes the needs
#' for PIMS because each animal (or unique capture history) has a real
#' parameter for each occasion.  The downside of the pacakge \code{mra} is that
#' every covariate is assumed to be time-varying and any factor variables like
#' \code{time} are coded manually as dummy variables for each level rather than
#' using the R facilities for handling factor variables in the formula to
#' create the design matrix.
#' 
#' In the \code{crm,cjs,js} functions in this package I have used basic idea in
#' \code{mra} but I have taken a different approach to model development that
#' allows for time-varying covariates but does not restrict each covariate to
#' be time-varying and factor variables are used as such which removes the need
#' to construct dummy variables; although the latter could still be used. First
#' an example is in order to explain how this all works. Consider the follwing
#' set of capture histories for small capture-recapture data set with 4 capture
#' occasions:
#' 
#' \preformatted{ 1001 0111 0011 }
#' 
#' To relate the current structure to the concept of PIMS I define the
#' following matrix
#' 
#' \preformatted{
#' 
#'  1 2 3 
#'  4 5 6 
#'  7 8 9 
#' }
#' 
#' If you think of these as \code{Phi} parameter indices, then 1 to 3 are
#' survivals for the intervals 1-2,2-3,3-4 for animal 1, and 4-6 and 7-9 are
#' the same for animals 2 and 3.  This matrix would have a row for each animal.
#' Now you'll notice that for animal 2 parameter 4 is not needed and for animal
#' 3, parameters 7 and 8 are not needed because they are prior to their entry
#' in the study.  While that is certainly true there is no harm in having them
#' and the advantage comes in being able to have a complete matrix in R rather
#' than having a triangular one.
#' 
#' So now we are finally getting to the description of what this function does.
#' It constructs a dataframe with a row for each animal-occasion.  Following on
#' with the example above, depending on how the arguments are set the following
#' dataframe could be constructed:
#' 
#' \preformatted{ row time Time cohort Cohort age Age initial.age 
#'                 1    1    0     1     0     0     0    0
#'                 2    2    1     1     0     1     1    0 
#'                 3    3    2     1     0     2     2    0 
#'                 4    1    0     2     1     0     0    0 
#'                 5    2    1     2     1     1     1    0 
#'                 6    3    2     2     1     2     2    0 
#'                 7    1    0     3     2     0     0    0
#'                 8    2    1     3     2     1     1    0 
#'                 9    3    2     3     2     2     2    0 
#' }
#' 
#' The fields starting with a lowercase character (time,cohort,age) are created
#' as factor variables and those with an uppercase are created as numeric
#' variables.  Note: the \code{age} field is bounded below by the minimum
#' \code{initial.age} to avoid creating factor levels with non-existent data
#' for occasions prior to first capture that are not used. For example, an
#' animal first caught on occasion 2 with an \code{initial.age=0} is
#' technically -1 on occasion 1 with a \code{time.interval} of 1. However, that
#' parameter would never be used in the model and we don't want a factor level
#' of -1.
#' 
#' A formula of ~time would create a design matrix with 3 columns (one for each
#' factor level) and ~Time would create one with 2 columns with the first being
#' an intercept and the second with the numeric value of Time.
#' 
#' Now here is the simplicity of it.  The following few expressions in R will
#' convert this dataframe into a matrix of real parameters (assuming
#' \code{beta=c(1,1,1)} that are structured like the square PIM matrix without
#' the use of PIMs.
#' 
#' \preformatted{ 
#' nocc=4
#' x=data.frame(ch=c("1001","0111","0011"),stringsAsFactors=FALSE)
#' beta=c(1,1,1) 
#' x.proc=process.data(x,model="cjs")
#' Phi.dmdf=make.design.data(x.proc)$Phi Phi.dm=create.dm(Phi.dmdf,~time)
#' Phimat=matrix(plogis(Phi.dm%*%beta),nrow=nrow(x),ncol=nocc-1,byrow=TRUE) 
#' }
#' 
#' Note that the order of the columns for \code{Phi.dmdf} differs slightly from
#' what is shown above. Also, \code{plogis} is an R function that computes the
#' inverse-logit. Once you have the matrix of \code{Phi} and \code{p} values
#' the calculation of the likelihood is straightforward using the formulation
#' of Pledger et al. (2003) (see \code{\link{cjs.lnl}}). The values in the
#' design dataframe are not limited to these fields.  The 2 arguments
#' \code{time.varying} and \code{fields} are vectors of character strings which
#' specify the names of the dataframe columns in \code{x} that should be
#' included. For example if \code{x} contained a field \code{sex} with the
#' values "M","F","M" for the 3 records in our example, and the argument
#' \code{fields=c("sex")} was used then a column named \code{sex} would be
#' included in design dataframe with the values
#' "M","M","M","F","F","F","M","M","M".  The value of the column \code{sex} in
#' \code{x} is repeated for each of the occasions for that animal(capture
#' history). Now if the value of the field changes for each occasion then we
#' use the argument \code{time.varying} instead.  To differentiate the values
#' in the dataframe \code{x} the columns are named with an occasion number.
#' For example, if the variable was named \code{cov} and it was to be used for
#' \code{Phi}, then the variables would be named \code{cov1,cov2,cov3} in
#' \code{x}. Let's say that x was structured as follows:
#' 
#' \preformatted{ 
#' ch   cov1 cov2 cov3 
#' 1001   1   0     1 
#' 0111   0   2     1 
#' 0011   0   0     0 
#' }
#' 
#' If you specified the argument \code{time.varying=c("cov")} then in the
#' design dataframe a field named \code{cov} would be created and the values
#' would be \code{1,0,1,0,2,1,0,0,0}. Thus the value is both animal and
#' occasion specific.  Had the covariate been used for \code{p} then they would
#' be named \code{cov2,cov3,cov4} because the covariate is for those occasions
#' for \code{p} whereas for \code{Phi} the covariate is labelled with the
#' occasion that begins the interval.  Note that unlike usage with \code{mark
#' in RMark} the covariate labelling does not change based on the value of
#' \code{time} which changes with the value of the argument \code{begin.time}.
#' Any number of fields can be specified in \code{fields} and
#' \code{time.varying} that are specified in \code{x}.
#' 
#' The input dataframe \code{x} has a few minor requirements on its structure.
#' First, it must contain a field called \code{ch} which contains the
#' capture-history as a string.  Note that in general strings are converted to
#' factor variables by default when they are put into a dataframe but as shown
#' above that can be controlled by the argument \code{stringsAsFactors=FALSE}.
#' The capture history should be composed only of 0 or 1 and they should all be
#' the same length (at present no error checking on this). Although it is not
#' necessary, the dataframe can contain a field named \code{freq} which
#' specifies the frequency of that capture history.  If the value of
#' \code{freq} is negative then these are treated as loss on capture at the
#' final capture occasion (last 1).  If \code{freq} is missing then a value of
#' 1 is assumed for all records.  Another optional field is \code{initial.age}
#' which specifies the age of the animal at the time it was first captured.
#' This field is used to construct the \code{age} and \code{Age} fields in the
#' design dataframe.  The default is to assume \code{initial.age=0} which means
#' the \code{age} is really time since first marked.  Any other fields in
#' \code{x} are user-specified and can be a combination of factor and numeric
#' variables that are either time-invariant or time-varying (and named
#' appropriately).
#' 
#' @param x processed dataframe from function \code{\link{process.data}}
#' @param parameter list with fields defining each values for each parameter;
#' as created by \code{\link{setup.parameters}}
#' @param time.varying vector of field names that are time-varying for this
#' parameter
#' @param fields character vector containing field names for variables in x to
#' be included in design matrix dataframe; if NULL all other than ch are
#' included
#' @return A dataframe with all of the individual covariate data and the
#' standard design data of time, Time, cohort, Cohort, age and Age; where
#' lowercase first letter implies a factor and uppercase is a numeric variable
#' for those variables.
#' @author Jeff Laake
#' @references Laake, J. and E. Rexstad (2007). RMark -- an alternative
#' approach to building linear models in MARK. Program MARK: A Gentle
#' Introduction. E. Cooch and G. C. White.
#' 
#' Pledger, S., K. H. Pollock, et al. (2003). Open capture-recapture models
#' with heterogeneity: I. Cormack-Jolly-Seber model. Biometrics 59(4):786-794.
create.dmdf=function(x,parameter,time.varying=NULL,fields=NULL)
{
   last= -parameter$num
   begin.num=parameter$begin+1
   chp = process.ch(x$data$ch)
   firstseen=chp$first
   lastseen=chp$last
   chmat = chp$chmat
#  requires field in each record called initial.age; if missing set to 0 for a
#  time since marked field
   if(is.null(x$data$initial.age)) x$data$initial.age=0
   nocc=x$nocc 
   time.intervals=x$time.intervals
   begin.time=x$begin.time
#  from begin.time compute time labels which are parameter dependent; time at
#  beginning of interval for Phi and time of occasion for p; these determine cohort
#  labels as well.
   times=begin.time+c(0,cumsum(time.intervals))
   cohort.levels=levels(factor(times))
   ntimes=length(times)-last
   occ=begin.num:(ntimes+begin.num-1)
   factor.times=factor(times[occ])
#  if there are any time varying covariates then construct the data matrix for
#  those covariates. That matrix is appended to other data in x that is non-time varying
   tcv=NULL
   if(!is.null(time.varying))
   {
     for (i in 1:length(time.varying))
     {
       vnames=paste(time.varying[i],times[occ],sep="")     
       if( !all(vnames %in% names(x$data)))
         stop("Missing time varying variable ",paste(vnames[!vnames%in%names(x$data)],collapse=","))
       if(i==1) 
         tcv=data.frame(as.vector(t(as.matrix(x$data[, vnames]))))
       else
         tcv=cbind(tcv,as.vector(t(as.matrix(x$data[, vnames]))))
       x$data=x$data[,!names(x$data)%in%vnames]
     }  
     names(tcv)=time.varying
   }  
#  create a list of dataframes - one for each entry in x.  The fields time, Time, cohort
#  Cohort, age and Age are created automatically.  Any time-invariant data in x are
#  repeated in each row of the data.  There is one row for each modelled occasion.            
   fx=function(.row){
        cohort=times[firstseen[.row]]
        Times=times[begin.num:(ntimes+begin.num-1)]-min(times[begin.num:(ntimes+begin.num-1)])
        if(firstseen[.row]>1)
           ages=times[begin.num:(ntimes+begin.num-1)]+x$data$initial.age[.row]-sum(time.intervals[1:(firstseen[.row]-1)])-times[1]
        else
           ages=times[begin.num:(ntimes+begin.num-1)]+x$data$initial.age[.row]-times[1]
        ages[ages<min.age]=min.age
		Y=chmat[.row,begin.num:(ntimes+begin.num-1)]
		Z = rep(NA,nocc)
		Z[firstseen[.row]:lastseen[.row]]=1
		Z = Z[begin.num:(ntimes+begin.num-1)]
        newdm.df=data.frame(time=factor.times,
                             cohort=rep(factor(cohort,levels=cohort.levels),ntimes),
                             Time=Times,Cohort=rep(cohort-begin.time,ntimes),age=factor(ages),Age=ages,Y=Y,Z=Z)     
   }
   if(begin.num==1)
      min.age=min(x$data$initial.age,na.rm=TRUE )
   else
      min.age=min(x$data$initial.age+c(time.intervals,max(time.intervals))[firstseen],na.rm=TRUE)   
   dm.df=sapply(1:nrow(x$data),fx,simplify=FALSE)  
#  Bind all the data into a single dataframe and attach time-varying covariates if any.
#  return dataframe.
   res=do.call("rbind",dm.df)
   if(is.null(fields))
      res=cbind(res,x$data[rep(1:nrow(x$data),each=ntimes),!names(x$data)%in%c("ch"),drop=FALSE])                             
   else
      res=cbind(res,x$data[rep(1:nrow(x$data),each=ntimes),names(x$data)%in%fields,drop=FALSE])  
   if("group"%in%names(res))
	   levels(res$group)=apply(x$group.covariates,1,paste,collapse="")
   if(!is.null(tcv))  
      res=cbind(res,tcv)
   return(res)
}
