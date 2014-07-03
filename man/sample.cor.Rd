\name{sample.cor}
\alias{sample.cor}
\title{
sample.cor
}
\description{
Random samples from a correlation matrix with replacement
}

\usage{
  sample.cor(x, n=1, shrink=F)
}

\arguments{
  \item{x}{matrix or data frame of observations}
  \item{n}{number of samples}
  \item{shrink}{whether to use the shrinkage estimate of correlation matrix}
}

\value{
  Vector of samples from correlation matrix
}

\details{
  This function randomly sampling from high-dimensional correlation matrices without generating the full correlation matrix. 
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
  
  ## ten samples from correlation matrix
  sample.pcor(longley, n=10) 
  
  ## ten samples using shrinkage estimation
  sample.cor(longly, n=10, shrink=T)
}
