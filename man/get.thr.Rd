\name{get.thr}
\alias{get.thr}

\title{
Select a threshold for sample correlation matrix
}

\description{
The function selects a threshold for sample correlation matrix.
}
\usage{
get.thr(C, n, p, max.num = 1, alpha = NULL, step = NULL)
}

\arguments{
  \item{C}{
sample correlation matrix of a design matrix.
}
  \item{n}{
the number of observations of the design matrix.
}
  \item{p}{
the number of variables of the design matrix.
}
  \item{max.num}{
the number of times for which the threshold selection procedure is repeated. Usually max.num==1 is used.
}
  \item{alpha}{
The level at which the false discovery rate is controlled. When alpha==NULL, it is set to be 1/sqrt(p).
}
  \item{step}{
the size of a step taken when screening the p(p-1)/2 off-diagonal elements of C.
}
}

\value{
\item{thr}{selected threshold.}
\item{thr.seq}{when max.num>1, the sequence of thresholds selected at each iteration.}
}

\references{
H. Cho and P. Fryzlewicz (2012) High-dimensional variable selection via tilting, Journal of the Royal Statistical Society Series B, 74: 593-622.
}

\author{
Haeran Cho
}