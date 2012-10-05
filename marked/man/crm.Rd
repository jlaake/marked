\name{crm}
\alias{crm}
\title{Capture-recapture model fitting function}
\usage{
  crm(data, ddl = NULL, begin.time = 1, model = "CJS",
    title = "", model.parameters = list(),
    design.parameters = list(), initial = NULL,
    groups = NULL, time.intervals = NULL, debug = FALSE,
    method = "BFGS", hessian = FALSE, accumulate = TRUE,
    chunk_size = 1e+07, control = NULL, refit = 1,
    itnmax = 5000, scale = NULL, run = TRUE, burnin = 100,
    iter = 1000, ...)
}
\arguments{
  \item{data}{Either the raw data which is a dataframe with
  at least one column named ch (a character field
  containing the capture history) or a processed dataframe}

  \item{ddl}{Design data list which contains a list element
  for each parameter type; if NULL it is created}

  \item{begin.time}{Time of first capture(release)
  occasion}

  \item{model}{Type of c-r model ("cjs" or "js" at
  present)}

  \item{title}{Optional title; not used at present}

  \item{design.parameters}{Specification of any grouping
  variables for design data for each parameter}

  \item{model.parameters}{List of model parameter
  specifications}

  \item{initial}{Optional vector of initial values for beta
  parameters; if named from previous analysis only relevant
  values are used}

  \item{groups}{Vector of names factor variables for
  creating groups}

  \item{time.intervals}{Intervals of time between the
  capture occasions}

  \item{method}{optimization method for function
  \code{optimx}}

  \item{debug}{if TRUE, shows optimization output}

  \item{hessian}{if TRUE, computes v-c matrix using
  hessian}

  \item{accumulate}{if TRUE, like capture-histories are
  accumulated to reduce computation}

  \item{chunk_size}{specifies amount of memory to use in
  accumulating capture histories and design matrices;
  amount used is 8*chunk_size/1e6 MB (default 80MB)}

  \item{control}{control string for optimization functions}

  \item{refit}{non-zero entry to refit}

  \item{itnmax}{maximum number of iterations for
  optimization}

  \item{scale}{vector of scale values for parameters}

  \item{run}{if TRUE, it runs model; otherwise if FALSE can
  be used to test model build components}

  \item{burnin}{number of iterations for mcmc burnin;
  specified default not realistic for actual use}

  \item{iter}{number of iterations after burnin for mcmc
  (not realistic default)}

  \item{...}{optional arguments passed to js or cjs and
  optimx}
}
\value{
  crm model object with class=("crm",submodel) where
  submodel is either "CJS" or "JS" at present.
}
\description{
  Fits user specified models to some types of
  capture-recapture wholly in R and not with MARK.  A
  single function that processes data, creates the design
  data, makes the crm model and runs it
}
\details{
  This function is operationally similar to the function
  \code{mark in RMark} in that is is a shell that calls
  several other functions to perform the following steps:
  1) \code{\link{process.data}} to setup data and
  parameters and package them into a list (processed
  data),2) \code{\link{make.design.data}} to create the
  design data for each parameter in the specified model, 3)
  \code{\link{create.dm}} to create the design matrices for
  each parameter based on the formula provided for each
  parameter, 4) call to the specific function for model
  fitting (now either \code{\link{cjs}} or
  \code{\link{js}}). As with \code{mark} the calling
  arguments for \code{crm} are a compilation of the calling
  arguments for each of the functions it calls (with some
  arguments renamed to avoid conflicts). If data is a
  processed dataframe (see \code{\link{process.data}}) then
  it expects to find a value for \code{ddl}.  Likewise, if
  the data have not been processed, then \code{ddl} should
  be NULL.  This dual calling structure allows either a
  single call approach for each model or alternatively for
  the data to be processed and the design data (\code{ddl})
  to be created once and then a whole series of models can
  be analyzed without repeating those steps.

  There are some optional arguments that can be used to set
  initial values and control other aspects of the
  optimization.  The optimization is done with the R
  package/function \code{optimx} and the arguments
  \code{method} and \code{hessian} are described with the
  help for that function.  In addition, any arguments not
  matching those for \code{cjs} (the ...) are passed to
  \code{optimx} allowing any of the other parameters to be
  set.  If you set \code{debug=TRUE}, then at each function
  evaluation (\code{\link{cjs.lnl}} the current values of
  the parameters and -2*log-likelihood value are output.

  In the current implementation, a logit link is used to
  constrain the parameters in the unit interval (0,1)
  except for probability of entry which uses an mlogit and
  N which uses a log link.  These could be generalized to
  use other link functions. Following the notation of MARK,
  the parameters in the link space are referred to as
  \code{beta} and those in the actual parameter space of
  \code{Phi} and \code{p} as reals.

  Initial values can be set in 2 ways.  To set a baseline
  intial value for the intercept of \code{Phi} \code{p} set
  those arguments to some real value in the open interval
  (0,1). All non-intercept beta parameters are set to zero.
  Alternatively, you can specify in \code{initial}, a
  vector of initial values for the beta parameters (on the
  logit scale).  This is most easily done by passing the
  results from a previous model run using the result list
  element \code{beta} as described below.  The code will
  match the names of the current design matrix to the names
  in \code{beta} and use the appropriate initial values.
  Any non-specified values are set to 0.  If there are no
  names associated with the \code{initial} vector then they
  are simply used in the specified order. If you do not
  specify initial values it is equivalent to setting
  \code{Phi} and \code{p} to 0.5.

  If you have a study with unequal time intervals between
  capture occasions, then these can be specified with the
  argument \code{time.intervals}.

  The argument \code{accumulate} defaults to \code{TRUE}.
  When it is \code{TRUE} it will accumulate common capture
  histories that also have common design and common fixed
  values (see below) for the parameters.  This will speed
  up the analysis because in the calculation of the
  likelihood (\code{\link{cjs.lnl}} it loops over the
  unique values. In general the default will be the best
  unless you have many capture histories and are using many
  individual covariate(s) in the formula that would make
  each entry unique. In that case there will be no effect
  of accumulation but the code will still try to
  accumulate. In that particular case by setting
  \code{accumulate=FALSE} you can skip the code run for
  accumulation.

  Most of the arguments controlling the fitted model are
  contained in lists in the arguments
  \code{model.parameters} and \code{design.parameters}
  which are similar to their counterparts in \code{mark inb
  RMark}. Each is a named list with the names being the
  parameters in the model (e.g., Phi and p in "cjs" and
  "Phi","p","pent","N" in "js"). Each named element is also
  a list containing various values defining the design data
  and model for the parameter. The elements of
  \code{model.parameters} can include \code{formula} which
  is an R formula to create the design matrix for the
  parameter and \code{fixed} is a matrix of fixed values as
  described below. The elements of \code{design.parameters}
  can include \code{time.varying}, \code{fields},
  \code{time.bins},\code{age.bins}, and \code{cohort.bins}.
  See \code{\link{create.dmdf}} for a description of the
  first 2 and \code{\link{create.dm}} for a description of
  the last 3.

  Real parameters can be set to fixed values using
  \code{fixed=x} where x is a matrix with 3 columns and any
  number of rows.  The first column specifies the
  particular animal (capture history) as the row number in
  the dataframe \code{x}.  The second specifies the capture
  occasion number for the real parameter to be fixed.  For
  \code{Phi} and \code{pent} these are 1 to \code{nocc}-1
  and for \code{p} they are 2 to \code{nocc} for "cjs" and
  1 to \code{nocc} for "js". This difference is due to the
  parameter labeling by the beginning of the interval for
  Phi (e.g., survival from occasion 1 to 2) and by the
  occasion for \code{p}.  For "cjs" \code{p} is not
  estimated for occasion 1. The third element in the row is
  the real value in the closed unit interval [0,1] for the
  fixed parameter.  This approach is completely general
  allowing you to fix a particular real parameter for a
  specific animal and occasion but it is a bit kludgy and
  later some other functionality will be included.  To set
  all of the real values for a particular occasion you can
  use the following example with the dipper data as a
  template:

  \code{model.parameters=list(Phi=list(,fixed=cbind(1:nrow(dipper),rep(2,nrow(dipper)),rep(1,nrow(dipper)))))}

  The above sets \code{Phi} to 1 for the interval between
  occasions 2 and 3 for all animals. At present there is no
  modification of the parameter count to address fixing of
  real parameters.

  Be cautious with this function at present.  It does not
  include many checks to make sure values are in the
  specified range of the data.  Normally this would not be
  a big problem but because \code{\link{cjs.lnl}} calls an
  external FORTRAN subroutine, if it gets a subscript out
  of bounds, it will cause R to terminate.  So make sure to
  save your workspace frequently if you use this function
  in its current implementation.
}
\examples{
# cormack-jolly-seber model
# fit 3 cjs models with crm
data(dipper)
dipper.proc=process.data(dipper,model="cjs",begin.time=1)
dipper.ddl=make.design.data(dipper.proc)
mod.Phit.pt=crm(dipper.proc,dipper.ddl,model.parameters=list(Phi=list(formula=~time),p=list(formula=~time)))
mod.Phit.pt
mod.Phidot.pdot=crm(dipper.proc,dipper.ddl,model.parameters=list(Phi=list(formula=~1),p=list(formula=~1)))
mod.Phidot.pdot
mod.Phisex.pdot=crm(dipper.proc,dipper.ddl,groups="sex",model.parameters=list(Phi=list(formula=~sex),p=list(formula=~1)))
mod.Phisex.pdot
# fit same 3 models with calls to mark
\donttest{
require(RMark)
mod0=mark(dipper,model.parameters=list(Phi=list(formula=~time),p=list(formula=~time)),output=FALSE)
summary(mod0,brief=TRUE)
mod1=mark(dipper,model.parameters=list(Phi=list(formula=~1),p=list(formula=~1)),output=FALSE)
summary(mod1,brief=TRUE)
mod2<-mark(dipper,groups="sex",model.parameters=list(Phi=list(formula=~sex),p=list(formula=~1)),output=FALSE)
summary(mod2,brief=TRUE)
}
# jolly seber model
crm(dipper,model="js",groups="sex",model.parameters=list(pent=list(formula=~sex),N=list(formula=~sex)),accumulate=FALSE)
\donttest{
mark(dipper,model="POPAN",groups="sex",model.parameters=list(pent=list(formula=~sex),N=list(formula=~sex)))
}
}
\author{
  Jeff Laake
}
\seealso{
  \code{\link{cjs}}, \code{\link{js}},
  \code{\link{make.design.data}},\code{\link{process.data}}
}
\keyword{models}

