#' Create design dataframes for crm
#' 
#' For each type of parameter in the analysis model (e.g, p, Phi), this
#' function creates design data for model fitting from the original data. 
#' These design data expand the original data record for a single(freq=1)/multiple(freq>1) individuals
#' to many records where each record is for a single occasion. The design data can
#' be specific to the parameter and a list of design data frames are returned with a dataframe
#' for each parameter.
#' 
#' For each record in the design data, default data fields are created based on the model. For example,
#' for CJS, the default fields include cohort, age, time which are factor variables and Cohort, Age, and Time
#' which are numeric versions of the same fields.  The values of these can be altered by values of 
#' begin.time and initial.ages set in the call to process.data. In addition if groups are identifed the
#' grouping factor variables are added to the design data. Any other fields in the data are repeated on each record
#' for the animal(s) (eg weight), unless they are defined as time.varying in which case the fields should be named
#' with the convention xn where x is the base name for the variable and n is the time value (eg, td1999, td2000,...).
#' For time varying fields the variable name in the design data is the base name and the value for the record is
#' the one for that occasion(time).
#' 
#' @param data Processed data list; resulting value from process.data
#' @param parameters Optional list containing a list for each type of parameter
#' (list of lists); each parameter list is named with the parameter name (eg
#' Phi); each parameter list can contain vectors named age.bins,time.bins and
#' cohort.bins, time.varying, fields \tabular{ll}{ \code{age.bins} \tab bins
#' for binning ages\cr \code{time.bins} \tab bins for binning times\cr
#' \code{cohort.bins} \tab bins for binning cohorts\cr \code{time.varying} \tab
#' vector of field names that are time varying for this parameter\cr
#' \code{fields} \tab vector of field names that are to be included that are
#' not time varying for this parameter\cr }
#' @return The function value is a list of data frames. The list contains a
#' data frame for each type of parameter in the model (e.g., Phi and p for
#' CJS).  The names of the list elements are the parameter names (e.g., Phi).
#' The structure of the dataframe depends on the calling arguments and the
#' model & data structure.
#' @author Jeff Laake
#' @export
#' @seealso \code{\link{process.data}},\code{\link{merge_design.covariates}}
#' @keywords utility
#' @examples
#' \donttest{
#' data(dipper)
#' dipper.proc=process.data(dipper)
#' ddl=make.design.data(dipper.proc)
#' }
"make.design.data" <-                                                      
function(data,parameters=list())
{
#
# Check validity of parameter list; stop if not valid
#
  if(!valid.parameters(data$model,parameters)) stop()
#
#  Add following elements based on type of model
#           begin            - index for compute.design.data
#           num              - number of parameters relative to number of occasions
#
  par.list=setup.parameters(data$model,check=TRUE)
  parameters=setup.parameters(data$model,parameters,data$nocc,check=FALSE,
          number.of.groups=dim(data$freq)[2])
  parameters=parameters[par.list]
  model.list=setup.model(data$model,data$nocc,data$mixtures)
#
# Create data for the each parameter in the model with age, year and cohort for each index
# This data matrix (design.data) is used below to create the design matrix from the formulas
# If age,cohort or year bins are given, use those.  Otherwise each is treated as a factor 
# wihtout binning.
full.design.data=vector("list",length=length(parameters))
#
# Create design data for crm models
#  Loop over each parameter in the model
   for(i in 1:length(parameters))
   {
#    Compute design data for parameters other than N
     if(names(parameters)[i]!="N")
	 {
		 if(is.null(parameters[[i]]$time.varying))
			 time.varying=NULL
		 else
			 time.varying=parameters[[i]]$time.varying
		 if(is.null(parameters[[i]]$fields))
			 fields=NULL
		 else
			 fields=parameters[[i]]$fields   
		 full.design.data[[i]]=create.dmdf(data,parameters[[i]],time.varying=time.varying,fields=fields)
	 }else
#    Compute design data for N
	 {
	    if(!is.null(data$group.covariates))
        {
           xlist=split(data$data,data$data[,names(data$group.covariates),drop=FALSE])
           xlist=xlist[as.vector(sapply(xlist,function(x) nrow(x)))>0]
        } else
        {                              
           xlist=data$data
        }
        numvar=sapply(data$data,is.numeric)
        numvar=numvar[names(numvar)!="freq"]
        if(any(numvar))
        {
           numvar=names(data$data[,names(data$data)!="freq"])[numvar]
           if(!is.null(data$group.covariates))
           {
              xmeans=sapply(xlist,function(x) sapply(subset(x,select=numvar),mean))
			  if(length(numvar)==1)
			  {
				  dd=data.frame(group=1:nrow(data$group.covariates))
				  dd[,numvar]=xmeans
			  }else
			  {
				  dd=data.frame(group=1:nrow(data$group.covariates))
				  dd=cbind(dd,t(xmeans))
			  }
              full.design.data[[i]]=merge(cbind(data.frame(group=1:nrow(data$group.covariates),data$group.covariates)),
					                         dd,by.x="group",all.x=TRUE)
#             full.design.data[[i]]=cbind(data$group.covariates,group=factor(1:nrow(data$group.covariates)),xmeans)
              row.names(full.design.data[[i]])=NULL
              names(full.design.data[[i]])=c("group",names(data$group.covariates),numvar[numvar!="group"])          
           }
           else
           {
              xmeans=colMeans(subset(data$data,select=numvar))         
              if(ncol(t(xmeans))==0)
                 full.design.data[[i]]=data.frame(N=1)
              else
              {
                 full.design.data[[i]]=data.frame(t(xmeans))  
                 row.names(full.design.data[[i]])=NULL
                 names(full.design.data[[i]])=numvar          
              }
           }
         }else
         {
           numvar=NULL
           if(!is.null(data$group.covariates))
           {
              full.design.data[[i]]=data$group.covariates
              full.design.data[[i]]$group=factor(row.names(data$group.covariates))
           }
           else
              full.design.data[[i]]=data.frame(N=1)
         }     
      } 
   }
   names(full.design.data)=names(parameters)
   return(full.design.data)
}

        
         
