\name{ex1.draw}
\docType{data}
\alias{ex1.draw}
\title{MCMC samples of Bayesian cluster model for a simulated dataset}
\description{
MCMC samples of clusterings from a Dirichlet process scale-location mixture model with normal components fitted to a simulated dataset, see \code{\link{ex1.data}}. True clusters are located at (+/- 2, +/- 2) with a standard deviation of 1.}
\usage{
data(ex1.draw)
}
\source{Wade, S. and Ghahramani, Z. (2015) Bayesian cluster analysis: Point estimation and credible balls. 
Submitted. arXiv:1505.03339.}
\format{
The matrix \code{ex1.draw} has 10,000 rows and 200 columns, with each row representing a MCMC posterior sample of the clustering of the 200 data points. 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datasets }
\examples{
data(ex1.data)
data(ex1.draw)
x=data.frame(ex1.data[,c(1,2)])
cls.true=ex1.data$cls.true
plot(x[,1],x[,2],xlab="x1",ylab="x2")
k=max(cls.true)
for(l in 2:k){
points(x[cls.true==l,1],x[cls.true==l,2],col=l)}

# Find representative partition of posterior
psm=comp.psm(ex1.draw)
ex1.VI=minVI(psm,ex1.draw,method=("all"),include.greedy=TRUE)
summary(ex1.VI)
plot(ex1.VI,data=x)

# Uncertainty in partition estimate
ex1.cb=credibleball(ex1.VI$cl[1,],ex1.draw)
summary(ex1.cb)
plot(ex1.cb,data=x)
}