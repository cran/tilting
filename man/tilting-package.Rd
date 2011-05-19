\name{tilting-package}
\alias{tilting-package}
\docType{package}
\title{
Variable selection via Tilted Correlation Screening algorithm
}

\description{
The package implements the tilted correlation screening algorithm proposed in Cho and Fryzlewicz (2011) for variable selection in high-dimensional linear regression. As a way of measuring the contribution of each variable to the response, the algorithm uses the "tilted correlation", which takes into account the high correlations among the variables in a data-driven way. 
}

\details{
\tabular{ll}{
Package: \tab tilting\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2011-05-19\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}

The main function of the package is \code{\link{tilting}}.

}
\author{Haeran Cho
Maintainer: Haeran Cho <h.cho1@lse.ac.uk>
}

\references{
H. Cho and P. Fryzlewicz (2011) "High-dimensional variable selection via tilting". To appear in Journal of the Royal Statistical Society: Series B.
}


\examples{

X<-matrix(rnorm(100*100), 100, 100) # 100-by-100 design matrix
y<-apply(X[,1:5], 1, sum)+rnorm(100) # first five variables are significant

tilt<-tilting(X, y, op=2)
tilt$active.hat # returns the finally selected variables

}
