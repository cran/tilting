\name{select.model}
\alias{select.model}

\title{
Select the final model
}
\description{
The function returns the final model as a subset of the active set chosen by the Tilted Correlation Screening algorithm, for which the extended BIC is minimised.
}
\usage{
select.model(bic.seq, active)
}
\arguments{
  \item{bic.seq}{
the sequence of extended BIC at each iteration.
}
  \item{active}{
the index set of selected variables by Tilted Correlation Screening algorithm.
}
}


\value{
The index set of finally selected variables is returned.
}

\author{
Haeran Cho
}

\references{
H. Cho and P. Fryzlewicz (2011) "High-dimensional variable selection via tilting". To appear in Journal of the Royal Statistical Society: Series B.

J. Chen and Z. Chen (2008) "Extended Bayesian information criteria for model selection with large model spaces," Biometrika, 95, 759-771.
}
