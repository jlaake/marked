\name{js.lnl}
\alias{js.lnl}
\title{Likelihood function for Jolly-Seber model using Schwarz-Arnason POPAN formulation}
\description{
For a given set of parameters and data, it computes -2*log Likelihood value but does not include 
data factorials. Factorials for unmarked are not needed but are included in final result
by \code{\link{js}} so the result matches output from MARK for the POPAN model.}
\usage{
js.lnl(par,imat,Phi.dm,p.dm,pent.dm,N.dm,Phi.fixed=NULL,p.fixed=NULL,
                 pent.fixed=NULL,debug=FALSE,time.intervals=NULL,nobstot)
}
\arguments{
  \item{par}{vector of parameter values}
  \item{imat}{list of vectors and matrices constructed by \code{\link{process.ch}} 
              from the capture history data.}
  \item{Phi.dm}{design matrix for Phi constructed by \code{\link{create.dm}} }
  \item{p.dm}{design matrix for p constructed by \code{\link{create.dm}} }
  \item{pent.dm}{design matrix for probability of entry constructed by \code{\link{create.dm}} }
  \item{N.dm}{design matrix for estimates of number of animals not caught from super-population constructed by \code{\link{create.dm}}}
  \item{Phi.fixed}{matrix with 3 columns: ch number(i), occasion number(j), fixed value(f)
                       to fix phi(i,j)=f  }
  \item{p.fixed}{matrix with 3 columns: ch number(i), occasion number(j), fixed value(f)
                       to fix p(i,j)=f  }
  \item{pent.fixed}{matrix with 3 columns: ch number(i), occasion number(j), fixed value(f)
                       to fix pent(i,j)=f  }
  \item{debug}{if TRUE will printout values of \code{par} and function value}
  \item{time.intervals}{intervals of time between occasions if not all 1}
  \item{nobstot}{number of unique caught at least once by group if applicable}
}
\details{
This functions uses \code{\link{cjs.lnl}} and then supplements with the remaining
calculations to compute the likelihood for the POPAN formulation (Arnason and Schwarz 1996) 
of the Jolly-Seber model.
}
\value{
 -2*log likelihood value (excluding data factorials)
}
\references{ 
Schwarz, C. J., and A. N. Arnason. 1996. A general methodology for the 
analysis of capture-recapture experiments in open populations. Biometrics 52:860-873.
}
\author{Jeff Laake}
