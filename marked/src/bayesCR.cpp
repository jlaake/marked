#include <Rcpp.h>
using namespace Rcpp;

/*
 * Sample a discrete uniform variable with unnormalized probability vector 'prob'
 */

int sampleDU(NumericVector ppp){
	double norm = std::accumulate(ppp.begin(), ppp.end(), 0.0);
	ppp = ppp/norm;
	NumericVector cdf = ppp;
	RNGScope scope;
	double U = as<double>(runif(1));
	int out = 1;
	if(U<= cdf[0]) return(out);
	else
	{
		for(int i=1; i<ppp.size(); i++){ 
			cdf(i) += cdf(i-1);
			if(U <= cdf[i]){
  				out = i+1;
  				return(out);
			}
		}
	}
}

/*
 * Indexing for beta.z estimation
 */

RcppExport SEXP makeZtildeIdx(SEXP ID, SEXP zvec){
	IntegerVector id(ID);
	IntegerVector z(zvec);
	int n = z.size();
	LogicalVector out(n);
	out(0) = TRUE;
	for(int i=1; i<n; i++){
    if(id(i) != id(i-1)) out(i) = TRUE;
    else if(z(i-1) == 1) out(i) = TRUE;
    else out(i) = FALSE;
	}
	return wrap(out);
}

RcppExport SEXP makeZCovIdx(SEXP ID, SEXP zvec){
  IntegerVector id(ID);
	IntegerVector z(zvec);
	int n = z.size();
	LogicalVector out(n);
	if(id(0) == id(1)) out(0) = TRUE;
  else out(0) = FALSE;
  out(n-1) = FALSE;
	for(int i=1; i<n-1; i++){
    if(id(i) != id(i-1) & id(i)==id(i+1)) out(i) = TRUE;
    else if(z(i) == 1 & id(i)==id(i+1)) out(i) = TRUE;
    else out(i) = FALSE;
	}
  
	return wrap(out);
}


/*
 * Sample the 'alive' or 'dead' state of an individual
 */

RcppExport SEXP sampleZ(SEXP ID, SEXP PVec, SEXP PhiVec, SEXP yVec){
	
// Setting up the variables
IntegerVector indiv(ID);
NumericVector p(PVec);
NumericVector phi(PhiVec);
IntegerVector y(yVec);
int n = y.size();
NumericMatrix pi(n,2);
NumericVector prob(2);
IntegerVector out(n);

// Starting forward loop
 if(y(0)==1) pi(0,1) = 1;
 else {
   pi(0,0) = 1-phi(0);
   pi(0,1) = phi(0)*(1-p(0));
 }
for(int i=1; i<n; i++)
{
  if(y(i)==1) pi(i,y(i)) = 1;
  else if(indiv(i) != indiv(i-1)) {
    pi(i,0) = 1-phi(i);
    pi(i,1) = phi(i)*(1-p(i));
  }
  else {
    pi(i,0) = pi(i-1,0) + pi(i-1,1)*(1-phi(i-1));
    pi(i,1) = pi(i-1,1)*phi(i-1)*(1-p(i));   
  }
}

// Starting backwrd pass with sampling
prob = pi(n-1,_);  
out(n-1) = sampleDU(pi(n-1,_))-1;
for(int i=n-2; i>=0; i--)
{
  prob = pi(i,_);
  if(indiv(i+1)==indiv(i)){
    prob(0) *= 1-out(i+1);
    prob(1) *= out(i+1)*phi(i+1)*(1-p(i+1))  + (1-out(i+1))*(1-phi(i+1));
  }
  out(i) = sampleDU(prob)-1;
}
return wrap(out);
    
}
