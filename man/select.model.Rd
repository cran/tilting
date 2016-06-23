\name{select.model}
\alias{select.model}

\title{
Select the final model
}
\description{
The function returns the final model as a subset of the active set chosen by Tilted Correlation Screening algorithm, for which the extended BIC is minimised.
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