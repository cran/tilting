\name{thresh}
\alias{thresh}

\title{
Hard-threshold a matrix
}
\description{
For a given matrix and a threshold, the function performs element-wise hard-thresholding based on the absolute value of each element.
}

\usage{
thresh(C, alph, eps = 1e-10)
}

\arguments{
  \item{C}{
a matrix on which the hard-thresholding is performed.
}
  \item{alph}{
threshold.
}
  \item{eps}{
effective zero.
}
}
\value{
Returns the matrix C after hard-thresholding.
}

\author{
Haeran Cho
}