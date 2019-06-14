\name{ex2.draw}
\docType{data}
\alias{ex2.draw}
\title{MCMC samples of Bayesian cluster model for a simulated dataset}
\description{
MCMC samples of clusterings from a Dirichlet process scale-location mixture model with normal components fitted to a simulated dataset, see \code{\link{ex2.data}}.
True clusters are located at (+/- 2, +/- 2) with a standard deviation of 1, 0.5, 1, and 1.5 in the first, second, third, and fourth quadrant respectively.}
\usage{
data(ex2.draw)
}
\source{Wade, S. and Ghahramani, Z. (2015) Bayesian cluster analysis: Point estimation and credible balls. 
Submitted. arXiv:1505.03339.}
\format{
The matrix \code{ex2.draw} has 10,000 rows and 200 columns, with each row representing a MCMC posterior sample of the clustering of the 200 data points contained in \code{\link{ex2.data}}. 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datasets }
\examples{
data(ex2.data)
data(ex2.draw)
x=data.frame(ex2.data[,c(1,2)])
cls.true=ex2.data$cls.true
plot(x[,1],x[,2],xlab="x1",ylab="x2")
k=max(cls.true)
for(l in 2:k){
points(x[cls.true==l,1],x[cls.true==l,2],col=l)}

# Find representative partition of posterior
psm=comp.psm(ex2.draw)
ex2.VI=minVI(psm,ex2.draw,method=("all"),include.greedy=TRUE)
summary(ex2.VI)
plot(ex2.VI,data=x)

# Uncertainty in partition estimate
ex2.cb=credibleball(ex2.VI$cl[1,],ex2.draw)
summary(ex2.cb)
plot(ex2.cb,data=x)
}