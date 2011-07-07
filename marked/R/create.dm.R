create.dm=function(x, formula, time.bins=NULL, cohort.bins=NULL, age.bins=NULL, chunk_size=1e7, remove.intercept=NULL)
##############################################################################
# create.dm - create design matrix with nch*(nocc-1) rows
#             where nch is number of capture histories and nocc is number of
#             occasions. 
#
# Arguments:
#
#    x              - design matrix dataframe created by create.dmdf 
#    formula        - formula for parameter
#    time.bins      - bins for times to reduce number of parameters by collapsing
#    cohort.bins    - bins for cohorts to reduce number of parameters by collapsing
#    age.bins       - bins for ages to reduce number of parameters by collapsing
#    chunk_size      - specifies amount of memory to use in creating design matrix
#                        use is 8*chunk_size/1e6 MB (default 80MB)
#
# Value:      - design matrix for fitting
#
##############################################################################
{
   if(!is.null(time.bins))
   {
      factime=factor(cut(as.numeric(levels(x$time)[x$time]),time.bins,include.lowest=TRUE))  
      if(any(is.na(factime)))
         stop(paste("Time bins do not span all values. Min time:", min(as.numeric(levels(x$time))),
           "Max time:", max(as.numeric(levels(x$time)))))
      if(length(levels(factime))==1) stop(paste("Need to specify at least 2 intervals in the time data",
        "Min time:", min(as.numeric(levels(x$time))),"Max time:", max(as.numeric(levels(x$time)))))
      x$time=factime
   }
   if(!is.null(cohort.bins))
   {
      faccohort=factor(cut(as.numeric(levels(x$cohort)[x$cohort]),cohort.bins,include.lowest=TRUE))  
      if(any(is.na(faccohort)))
         stop(paste("Cohort bins do not span all values. Min cohort:", min(as.numeric(levels(x$cohort))),
           "Max cohort:", max(as.numeric(levels(x$cohort)))))
      if(length(levels(faccohort))==1) stop(paste("Need to specify at least 2 intervals in the cohort data",
        "Min cohort:", min(as.numeric(levels(x$cohort))),"Max cohort:", max(as.numeric(levels(x$cohort)))))
      x$cohort=faccohort
   }
   if(!is.null(age.bins))
   {
      facage=factor(cut(as.numeric(levels(x$age)[x$age]),age.bins,include.lowest=TRUE))  
      if(any(is.na(facage)))
         stop(paste("Age bins do not span all values. Min age:", min(as.numeric(levels(x$age))),
           "Max age:", max(as.numeric(levels(x$age)))))
      if(length(levels(facage))==1) stop(paste("Need to specify at least 2 intervals in the age data",
              "Min age:", min(as.numeric(levels(x$age))),"Max age:", max(as.numeric(levels(x$age)))))
      x$age=facage
   }
#  Create design matrix from formula and data; do so based on chunks of data to reduce space requirements
   mm=model.matrix(formula,x[1,,drop=FALSE])
   npar=ncol(mm)
   nrows=nrow(x)
   upper=0
   dm=Matrix(0,nrow=nrows,ncol=npar)
   pieces=floor(npar*nrows/chunk_size+1)
   rows_in_piece=floor(nrows/pieces)
   if(npar*nrows>chunk_size)
   {
      for(i in 1:pieces)
	  {
		  lower=(i-1)*rows_in_piece+1
		  upper=i*rows_in_piece
		  dm[lower:upper,]=as(model.matrix(formula,x[lower:upper,,drop=FALSE]),"sparseMatrix") 
	  }
   }
   if(upper<nrow(x))
	   dm[(upper+1):nrow(x),]=as(model.matrix(formula,x[(upper+1):nrow(x),,drop=FALSE]),"sparseMatrix")    
#  Remove any unused columns; this is slower but uses less memory
   select=vector("logical",length=npar)
	for (i in 1:npar)
	   select[i]=any(dm[,i]!=0)
   if(!is.null(remove.intercept)&&remove.intercept)select[1]=FALSE 
#  Return dm with selected columns
   colnames(dm)=colnames(mm)
   return(dm[,select,drop=FALSE])
}


