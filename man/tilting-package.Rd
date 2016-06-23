\name{tilting-package}
\alias{tilting-package}
\docType{package}
\title{
Variable Selection via Tilted Correlation Screening Algorithm
}

\description{
Implements an algorithm for variable selection in high-dimensional linear regression using the "tilted correlation", a new way of measuring the contribution of each variable to the response which takes into account high correlations among the variables in a data-driven way. 
}

\details{
\tabular{ll}{
Package: \tab tilting\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2016-06-20\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}

The main function of the package is \code{\link{tilting}}.

}
\author{Haeran Cho
Maintainer: Haeran Cho <haeran.cho@bristol.ac.ukk>
}

\references{
H. Cho and P. Fryzlewicz (2012) High-dimensional variable selection via tilting, Journal of the Royal Statistical Society Series B, 74: 593-622.
}


\examples{

X <- matrix(rnorm(100*100), 100, 100) # 100-by-100 design matrix
y <- apply(X[,1:5], 1, sum)+rnorm(100) # first five variables are significant

tilt <- tilting(X, y, op=2)
tilt$active.hat # returns the finally selected variables

}
