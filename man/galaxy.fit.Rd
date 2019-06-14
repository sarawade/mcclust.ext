\name{galaxy.fit}
\docType{data}
\alias{galaxy.fit}
\title{Fitted density values from a Dirichlet process mixture model for the galaxy dataset}
\description{
Fitted density values of a Dirichlet process scale-location mixture model with normal components fitted to the \code{\link{galaxies}} dataset.}
\usage{
data(galaxy.fit)
}
\source{
Roeder, K. (1990) Density estimation with confidence sets exemplified by superclusters and voids in the galaxies, Journal of the American Statistical Association, 85: 617-624.

Wade, S. and Ghahramani, Z. (2015) Bayesian cluster analysis: Point estimation and credible balls. 
Submitted. arXiv:1505.03339.}
\format{
\enumerate{
\item The data points are contained in the first column \code{x} of length 82, see \code{\link{galaxies}} for more information. 
\item The second column \code{fx} of length 82 contains the density estimate from the Dirichlet process mixture at each data point. 
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datasets }
\examples{
data(galaxy.fit)
x=galaxy.fit$x
fx=galaxy.fit$fx
plot(x,fx,xlab="x",ylab="f(x)")
}