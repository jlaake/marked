marked - mark-recapture analysis
=======================================

We developed this R package for analysis with marked animals in 
contrast to the R package unmarked (Fiske and Chandler 2011) that 
focuses on analysis with unmarked animals. The original impetus 
for the package was to implement the CJS model using the hierarchical 
likelihood construction described by Pledger et al. (2003) and to 
improve on execution times with RMark/MARK (White and Burnham 1999;Laake 2013)
for analysis of our own large data sets with many time-varying 
individual (animal-specific) covariates. Subsequently, we implemented 
the Jolly-Seber model with the Schwarz and Arnason (1996) POPAN 
structure where the hierarchical likelihood construction idea extended 
to the entry of animals into the population. We also added a 
Bayesian Markov Chain Monte Carlo (MCMC) implementation of the 
CJS model based on the approach used by Albert and Chib (1993) for 
analyzing binary data with a probit regression model. 

This package requires installation of ADMB if you set use.admb=TRUE. 
See [readme.txt](https://github.com/jlaake/marked/blob/master/marked/inst/README.txt) for installation instructions.

Download [Windows package binary](https://docs.google.com/folder/d/0B77g1ScdUwVeOVJNUVVGS0YtWE0/edit). From link, browse to marked and then click on
the version of the package you want. You should see a listing of the package contents as files.  Select File/Download. 
To install in R, from the R menu, use Packages\Install from Local Zip file and browse to location of downloaded zip. 

Download [package source files](https://github.com/jlaake/marked/archive/master.zip)

The following are references from above:

Albert, J. and Chib, S. (1993). Bayesian-analysis of binary and polychotomous
reponse data. Journal of the American Statistical Association, 88(422):669:679.

Fiske, I. J. and Chandler, R. B. (2011). unmarked : An R Package for tting
hierarchical models of wildlife occurrence and abundance. 43(10):1:23.

Laake, J.L. (2013) RMark : An R Interface for Analysis of Capture-Recapture Data
with MARK. AFSC Processed Rep 2013-01, 25p. Alaska Fish. Sci. Cent., NOAA,
Natl. Mar. Fish. Serv., 7600 Sand Point Way NE, Seattle WA 98115.

Pledger, S., Pollock, K. H., and Norris, J. L. (2003). Open capture-recapture models
with heterogeneity: I. Cormack-Jolly-Seber model. Biometrics, 59(4):786:794.

Schwarz, C. J. and Arnason, A. N. (1996). A general methodology for the analysis
of capture-recapture experiments in open populations. Biometrics, 52(3):860:873.
16

White, G. C. and Burnham, K. P. (1999). Program MARK: survival estimation from
populations of marked animals. Bird Study, 46:120:139.
