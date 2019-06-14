\name{VI}
\alias{VI}
\alias{VI.lb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute the posterior expected Variation of Information or the modified version from swapping log and expectation
}
\description{
Based on MCMC samples of partitions, computes the posterior expected Variation of Information or the modified Variation of Information which switches the log and expectation.
}
\usage{

VI(cls, cls.draw)

VI.lb(cls, psm)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cls}{a matrix of partitions where the posterior expected (modified) Variation of Information is to be evaluated. Each row corresponds to a clustering of \code{ncol(cls)} data points.
}
  \item{cls.draw}{a matrix of MCMC samples of partitions. Each row corresponds to a clustering of \code{ncol(cls)} data points.
}
\item{psm}{a posterior similarity matrix, which can be obtained from MCMC samples of clusterings through a call to \code{comp.psm}.}
}
\details{
The Variation of Information (VI) between two clusterings is defined as the sum of the entropies minus two times the mutual information. Computation of the posterior expected VI can be expensive, as it requires a Monte Carlo estimate. The modified posterior expected VI, obtained by swapping the log and expectation, is much more computationally efficient as it only depends on the posterior through the posterior similarity matrix. From Jensen's inequality, the problem of finding the optimal partition which minimizing the posterior expected modified VI can be viewed as minimizing a lower bound to the posterior expected VI. 
}
\value{
vector of length \code{nrow(cls)} of the posterior expected (modified) VI.
}
\references{
Meila, M. (2007) Bayesian model based clustering
procedures, \emph{Journal of Multivariate Analysis} \bold{98}, 873--895.

Wade, S. and Ghahramani, Z. (2015) Bayesian cluster analysis: Point estimation and credible balls. 
Submitted. arXiv:1505.03339.
}
\author{
Sara Wade, \email{sara.wade@eng.cam.ac.uk}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{minVI}} which locates the partition that minimizes the posterior expected modified VI.
}
\examples{
data(ex2.data)
x=data.frame(ex2.data[,c(1,2)])
cls.true=ex2.data$cls.true
plot(x[,1],x[,2],xlab="x1",ylab="x2")
k=max(cls.true)
for(l in 2:k){
points(x[cls.true==l,1],x[cls.true==l,2],col=l)}

# Find representative partition of posterior
data(ex2.draw)
psm=comp.psm(ex2.draw)
ex2.VI=minVI(psm,ex2.draw,method=("all"),include.greedy=TRUE)
summary(ex2.VI)

# Compute posterior expected VI for each partition
VI(ex2.VI$cl,ex2.draw)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
