"make.design.data" <-                                                      
function(data,parameters=list())
{
#------------------------------------------------------------------------------------------------------
# make.design.data -  creates a design dataframe that is used to construct the design matrix for mark
#                     in make.mark.model
#
# Arguments:
#
#    data             - data list after using process.data
#    parameters       - list with an element for each parameter
#                       each element is a list with age.bins, time.bins and cohort.bins
#                          age.bins         - bins for grouping ages
#                          time.bins        - bins for grouping times
#                          cohort.bins      - bins for grouping cohorts
#                          time.varying     - vector of field names that are time varying for this parameter
#                          fields           - vector of field names to be included in design data that are not time varying
#
# Value:
#    full.design.data - list of design data frames for each type of parameter in the model
#
#
# Functions used: setup.parameters, compute.design.data, valid.parameters, setup.model
#
#----------------------------------------------------------------------------------------------------
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
# Create a data matrix for the each parameter in the model with age, year and cohort for each index
# This data matrix (design.data) is used below to create the design matrix from the formulas
# If age,cohort or year bins are given, use those.  Otherwise each is treated as a factor 
# wihtout binning.
#
# 10 Jan 06 ; added pim.type argument in call to compute.design.data
#
full.design.data=vector("list",length=length(parameters))
#
# Create design data for crm models
#
   for(i in 1:length(parameters))
   {
     if(names(parameters)[i]=="N")
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
              full.design.data[[i]]=merge(cbind(data.frame(group=1:nrow(data$group.covariates),data$group.covariates)),
                              t(xmeans),by.x="group",all.x=TRUE)
#             full.design.data[[i]]=cbind(data$group.covariates,group=factor(1:nrow(data$group.covariates)),xmeans)
              row.names(full.design.data[[i]])=NULL
              names(full.design.data[[i]])=c("group",names(data$group.covariates),numvar[numvar!="group"])          
           }
           else
           {
              xmeans=mean(subset(data$data,select=numvar))         
              if(ncol(t(xmeans))==0)
                 full.design.data[[i]]=data.frame(N=1)
              else
              {
                 full.design.data[[i]]=data.frame(xmeans)  
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
     } else
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
     }
   }
   names(full.design.data)=names(parameters)
   return(full.design.data)
}

        
         