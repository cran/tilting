\name{lse.beta}
\alias{lse.beta}

\title{
Compute the least squares estimate on a given index set
}
\description{
The function returns an estimate of the coefficient vector for a linear regression problem by setting the coefficients corresponding to a given index set to be the least squares estimate and the rest to be equal to zero.
}
\usage{
lse.beta(X, y, active = NULL)
}

\arguments{
  \item{X}{
design matrix.
}
  \item{y}{
response vector.
}
  \item{active}{
the index set on which the least squares estimate is computed. 
}
}

\value{
An estimate of the coefficient vector is returned as above. 
If active=NULL, a vector of zeros is returned.
}

\author{
Haeran Cho
}
