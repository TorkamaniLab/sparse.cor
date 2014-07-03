\name{sparse.pcor}
\alias{sparse.pcor}
\title{
sparse.pcor
}
\description{
Calculates a sparse partial correlation matrix using shrinkage estimation
}

\usage{
  sparse.pcor(x, threshold=0.1, n.pieces=10)
}

\arguments{
  \item{x}{matrix or data frame of data.}
  \item{threshold}{threshold for retaining correlation}
  \item{n.pieces}{sparse partial correlation matrix is calculated in n.pieces^2 chunks}
}

\value{
  A sparse matrix of partial orrelations
}

\details{
}

\author{
Eric Kramer <ericransomkramer@gmail.com>
}
\references{
~~ Literature or other references for background information ~~
}

\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
  ## longley data set is highly correlated
  
  ## generates sparse partial correlation matrix
  
  sigma = sparse.pcor(longley)
  
  ## sample to find distribution and then threshold
  
  s = sample.pcor(longley, n=100)
  threshold = quantile(s, probs=c(0.75))
  sigma = sparse.pcor(longley, threshold=threshold)

  ## use smaller pieces to reduce memory usage
  
  sigma = sparse.pcor(longley, n.pieces=15)
  
}
