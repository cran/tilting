\name{tilting}
\alias{tilting}

\title{
Variable selection via Tilted Correlation Screening algorithm
}

\description{
Given a design matrix and a response vector, the aim is to select the variables which have significant contributions to the response variable. The function implements the tilted correlation screening algorithm, which selects a threshold for the sample correlation matrix, computes an adaptive measure for the contribution of each variable to the response variable based on the thus-thresholded correlation matrix, and chooses one variable at each iteration. After iteratively applying the algorithm, an "active" set of variables is returned and the extended BIC is adopted to select the final model from the active set.
}

\usage{
tilting(X, y, thr.step = NULL, thr.rep = 1, max.size = NULL, max.count = NULL, 
op = 2, bic.gamma = 1, eps = 1e-10)
}

\arguments{
  \item{X}{
design matrix with the rows representing the observations and the columns representing the variables.
}
  \item{y}{
response vector.
}
  \item{thr.step}{
step size used for the threshold selection. When thr.step=NULL, it is chosen automatically.
}
  \item{thr.rep}{
the number of times for which the threshold selection procedure is repeated. As described in H. Cho and P. Fryzlewicz (2011), the procedure involves generating a random matrix of the same size as X and thr.rep controls how many times such a procedure is repeated. 
}
  \item{max.size}{
the maximum number of the variables conditional on which the contribution of each variable to the response is measured (when max.size=NULL, it is set to be half the number of observations).
}
  \item{max.count}{
the maximum number of iterations for the tilted correlation screening algorithm. It determines the number of variables to be included in the active set.
}
  \item{op}{
when op=1, rescaling 1 is used to compute the tilted correlation. If op=2, rescaling 2 is used. See H. Cho and P. Fryzlewicz (2011) for the details of the rescaling methods. 
}
  \item{bic.gamma}{
parameter which determines the amount of penalty imposed on the model complexity in the extended BIC.
}
  \item{eps}{
effective zero.
}
}

\value{
\item{active}{active set which contains the variables returned by the tilted correlation screening algorithm in the order of selection.}
\item{thr.seq}{a sequence of thresholds selected over the iterations.}
\item{bic.seq}{extended BIC computed over the iterations.}
\item{active.hat}{finally chosen variables from the active set using the extended BIC.}
}

\references{
H. Cho and P. Fryzlewicz (2011) "High-dimensional variable selection via tilting". To appear in Journal of the Royal Statistical Society: Series B.

J. Chen and Z. Chen (2008) "Extended Bayesian information criteria for model selection with large model spaces," Biometrika, 95, 759-771.
}

\author{
Haeran Cho
}

\examples{

X<-matrix(rnorm(100*100), 100, 100) # 100-by-100 design matrix
y<-apply(X[,1:5], 1, sum)+rnorm(100) # first five variables are significant

tilt<-tilting(X, y, op=2)
tilt$active.hat # returns the finally selected variables

}
