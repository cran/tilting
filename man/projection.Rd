\name{projection}
\alias{projection}

\title{
Compute the projection matrix onto a given set of variables
}
\description{
The function computes the projection matrix onto a set of columns of a given matrix. 
}
\usage{
projection(X, active = NULL)
}

\arguments{
  \item{X}{
a matrix containing the columns onto which the projection matrix is computed.
}
  \item{active}{
an index set of the columns of X.
}
}

\value{
Returns the projection matrix onto the columns of "X" whose indices are included in "active". When active==NULL, a null matrix is returned.
}

\author{
Haeran Cho
}

