\name{sample.pcor}
\alias{sample.pcor}
\title{
sample.pcor
}
\description{
Randomly samples from a partial correlation matrix with replacement
}

\usage{
  sample.pcor(x, n=1)
}

\arguments{
  \item{x}{matrix or data frame of observations}
  \item{n}{number of samples}
}

\value{
  Vector of samples from partial correlation matrix
}

\details{
  This function for the random sampling from high-dimensional partial correlation matrices without generating the full parital correlation matrix. 
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
  
  ## ten samples from partial correlation matrix
  sample.pcor(longley, n=10) 
  
}
