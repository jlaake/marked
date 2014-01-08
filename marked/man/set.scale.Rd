\name{set.scale}
\alias{scale.dm}
\alias{scale.par}
\alias{set.scale}
\alias{unscale.par}
\title{Scaling functions}
\usage{
set.scale(pars,model_data,scale)

	        scale.dm(model_data,scale)

         scale.par(par,scale)

         unscale.par(par,scale)
}
\arguments{
  \item{pars}{character vector of parameter names}

  \item{par}{list of parameter vectors or vector of
  parameter values}

  \item{scale}{list or vector of parameter scales}

  \item{model_data}{list of data/design objects}
}
\value{
List of scale values for set.scale, model.data with scaled
design matrices for scale.dm, vector of scaled parameter
values for scale.par, and list of unscaled parameter
vectors for unscale.par
}
\description{
Set scale, scale dm and scale/unscale parameters
}
\author{
Jeff Laake <jeff.laake@noaa.gov>
}

