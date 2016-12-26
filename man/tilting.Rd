\name{tilting}
\alias{tilting}

\title{
Variable selection via Tilted Correlation Screening algorithm
}

\description{
Given a design matrix and a response vector, the function selects a threshold for the sample correlation matrix, computes an adaptive measure for the contribution of each variable to the response variable based on the thus-thresholded sample correlation matrix, and chooses a variable at each iteration. Once variables are selected in the "active" set, the extended BIC is used for the final model selection.
}

\usage{
tilting(X, y, thr.step = NULL, thr.rep = 1, max.size = NULL, max.count = NULL,
op = 2, bic.gamma = 1, eps = 1e-10)
}

\arguments{
  \item{X}{
design matrix.
}
  \item{y}{
response vector.
}
  \item{thr.step}{
a step size used for threshold selection. When thr.step==NULL, it is chosen automatically.
}
  \item{thr.rep}{
the number of times for which the threshold selection procedure is repeated.
}
  \item{max.size}{
the maximum number of the variables conditional on which the contribution of each variable to the response is measured (when max.size==NULL, it is set to be half the number of observations).
}
  \item{max.count}{
the maximum number of iterations.
}
  \item{op}{
when op==1, rescaling 1 is used to compute the tilted correlation. If op==2, rescaling 2 is used.
}
  \item{bic.gamma}{
a parameter used to compute the extended BIC.
}
  \item{eps}{
an effective zero.
}
}

\value{
\item{active}{active set containing the variables selected over the iterations.}
\item{thr.seq}{a sequence of thresholds selected over the iterations.}
\item{bic.seq}{extended BIC computed over the iterations.}
\item{active.hat}{finally chosen variables using the extended BIC.}
}

\references{
H. Cho and P. Fryzlewicz (2012) High-dimensional variable selection via tilting, Journal of the Royal Statistical Society Series B, 74: 593-622.
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
