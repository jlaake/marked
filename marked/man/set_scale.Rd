% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scale.r
\name{set_scale}
\alias{set_scale}
\alias{scale_dm}
\alias{scale_par}
\alias{unscale_par}
\title{Scaling functions}
\usage{
set_scale(pars,model_data,scale)
 
        scale_dm(model_data,scale)

         scale_par(par,scale)

         unscale_par(par,scale)
}
\arguments{
\item{pars}{character vector of parameter names}

\item{model_data}{list of data/design objects}

\item{scale}{list or vector of parameter scales}

\item{par}{list of parameter vectors or vector of parameter values}
}
\value{
List of scale values for set_scale, model.data with scaled design matrices for scale_dm,
vector of scaled parameter values for scale_par, and list of unscaled parameter vectors for unscale_par
}
\description{
Set scale, scale dm and scale/unscale parameters
}
\author{
Jeff Laake
}
