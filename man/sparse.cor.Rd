\name{sparse.cor}
\alias{sparse.cor}
\title{
sparse.cor
}
\description{
Calculates a sparse (pearson) correlation matrix
}

\usage{
  sparse.cor(x, threshold=0.9, n.pieces=10, shrink=F)
}

\arguments{
  \item{x}{matrix or data frame of data.}
  \item{threshold}{threshold for retaining correlation}
  \item{n.pieces}{sparse correlation matrix is calculated in n.pieces^2 chunks}
  \item{shrink}{whether to use the shrinkage estimate of correlation matrix}
}

\value{
  A sparse matrix of correlations
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
  
  ## generates sparse correlation matrix
  
  sigma = sparse.cor(longley)
  
  ## sample to find distribution and then threshold
  
  s = sample.cor(longley, n=100)
  threshold = quantile(s, probs=c(0.75))
  sigma = sparse.cor(longley, threshold=threshold)

  ## use smaller pieces to reduce memory usage
  
  sigma = sparse.cor(longley, n.pieces=15)

}
