create.dmdf=function(x,parameter,time.varying=NULL,fields=NULL)
##############################################################################
# create.dmdf - create design matrix dataframe with nch*(nocc-1) rows
#             where nch is number of capture histories and nocc is number of
#             occasions. 
#
# Arguments:
#
#    x              - processed dataframe
#    parameter      - list with fields defining parameter
#    time.varying   - character vector containing base names for variables in x
#                     to be used to construct time-varying covariates
#    fields         - character vector containing field names for variables in x
#                     to be included in design matrix dataframe; if NULL all other than
#                     ch are included
#
# Value:      - design matrix dataframe for model.matrix application (create.dm)
##############################################################################
{
   last= -parameter$num
   begin.num=parameter$begin+1
   first=process.ch(x$data$ch)$first
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
        cohort=times[first[.row]]
        Times=times[begin.num:(ntimes+begin.num-1)]-min(times[begin.num:(ntimes+begin.num-1)])
        if(first[.row]>1)
           ages=times[begin.num:(ntimes+begin.num-1)]+x$data$initial.age[.row]-sum(time.intervals[1:(first[.row]-1)])-times[1]
        else
           ages=times[begin.num:(ntimes+begin.num-1)]+x$data$initial.age[.row]-times[1]
        ages[ages<min.age]=min.age
        newdm.df=data.frame(time=factor.times,
                             cohort=rep(factor(cohort,levels=cohort.levels),ntimes),
                             Time=Times,Cohort=rep(cohort-begin.time,ntimes),age=factor(ages),Age=ages)     
   }
#   if(parameter=="Phi" | js)
   if(begin.num==1)
      min.age=min(x$data$initial.age,na.rm=TRUE )
   else
      min.age=min(x$data$initial.age+c(time.intervals,max(time.intervals))[first],na.rm=TRUE)   
   dm.df=sapply(1:dim(x$data)[1],fx,simplify=FALSE)  
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
