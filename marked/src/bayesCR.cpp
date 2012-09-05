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
    pi(0,y(0)) = 1;
    for(int i=1; i<n; i++)
    {
      if(y(i)==1) 
        pi(i,y(i)) = 1;
      else
      {
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
        prob(1) *= out(i+1)*phi(i)*(1-p(i+1))  + (1-out(i+1))*(1-phi(i));
      }
      out(i) = sampleDU(prob)-1;
    }
    return wrap(out);
    
}
