\name{ex2.data}
\docType{data}
\alias{ex2.data}
\title{A simulated dataset from a mixture of normals}
\description{
A simulated dataset from a mixture of four normals. True clusters are located at (+/- 2, +/- 2) with a standard deviation of 1, 0.5, 1, and 1.5 in the first, second, third, and fourth quadrant respectively.}
\usage{
data(ex2.data)
}
\source{Wade, S. and Ghahramani, Z. (2015) Bayesian cluster analysis: Point estimation and credible balls. 
Submitted. arXiv:1505.03339.}
\format{
\enumerate{
\item The data points are contained in the first two columns \code{x1} and \code{x2} of length 200; the 200 data points were simulated from a mixture of four normals at locations (+/- 2, +/- 2) with a standard deviation of 1, 0.5, 1, and 1.5 in the first, second, third, and fourth quadrant respectively. 
\item The third column \code{cls.true} of length 200 contains the true clustering of the 200 data points.  
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ datasets }
\examples{
data(ex2.data)
x=ex2.data[,c(1,2)]
cls.true=ex2.data$cls.true
plot(x[,1],x[,2],xlab="x1",ylab="x2")
k=max(cls.true)
for(l in 2:k){
points(x[cls.true==l,1],x[cls.true==l,2],col=l)}
}