\name{galaxy.pred}
\docType{data}
\alias{galaxy.pred}
\title{Predicted densty values from a Dirichlet process mixture model for the galaxy dataset}
\description{
Predicted density values at a grid of new data points from a Dirichlet process scale-location mixture model with normal components fitted to the \code{\link{galaxies}} dataset.}
\usage{
data(galaxy.pred)
}
\source{
Roeder, K. (1990) Density estimation with confidence sets exemplified by superclusters and voids in the galaxies, Journal of the American Statistical Association, 85: 617-624.

Wade, S. and Ghahramani, Z. (2015) Bayesian cluster analysis: Point estimation and credible balls. 
Submitted. arXiv:1505.03339.}
\format{
\enumerate{ 
\item The column \code{x} of length 141 contains a grid of new data points from \code{5} to \code{40} by \code{0.25}.
\item The column \code{fx} of length 141 contains the density estimate from the Dirichlet process mixture at each new data point.
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datasets }
\examples{
data(galaxy.fit)
x=galaxy.fit$x
fx=galaxy.fit$fx
data(galaxy.pred)
xgrid=galaxy.pred$x
fxgrid=galaxy.pred$fx
plot(xgrid,fxgrid,xlab="x",ylab="f(x)",type="l")
points(x,fx)
}