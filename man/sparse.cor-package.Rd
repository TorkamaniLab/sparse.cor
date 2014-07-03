\name{sparse.cor-package}
\alias{sparse.cor-package}
\docType{package}
\title{
sparse.cor
}
\description{
Methods for creating and sampling from sparse correlation and sparse partial correlation matrices. 
}
\details{
\tabular{ll}{
Package: \tab sparse.cor\cr
Type: \tab Package\cr
Version: \tab 0.1\cr
Date: \tab 2014-06-27\cr
License: \tab What license is it under?\cr
}
Many large data sets contain hundreds of thousands of variables. Correlation matrices for these datasets can be used for clustering and infering dependencies between variables. However, these correlation matrices typically require excessive amounts of memory. This package implements methods for generating sparse correlation matrices (i.e. matrices which retain only correlations above a certain threshold) in a memory-efficient manner.

Additionally, we re-implement the partial correlation method from the corpcor package, which allows for the extremely rapid calculation of large partial correlation matrices. For instance, we can calculate a 100,000 x 100,000 partial correlation matrix on a MacBook pro in under 2 hours, while using less than 2 Gb of ram. 
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

}
