\name{greedy}
\alias{greedy}

\title{
Optimizes the posterior expected loss with the greedy search algorithm
}
\description{
Finds a representative partition of the posterior by minimizing the posterior expected loss with possible loss function of Binder's loss, the Variation of Information, and the modified Variation of Information through a greedy search algorithm.
}
\usage{
greedy(psm, cls.draw = NULL, loss = NULL, start.cl = NULL, maxiter = NULL, L = NULL, suppress.comment = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{psm}{a posterior similarity matrix, which can be obtained from MCMC samples of clusterings through a call to \code{comp.psm}.}
 \item{cls.draw}{a matrix of the MCMC samples of clusterings of the \code{ncol(cls)} data points that have been used to compute \code{psm}. Note: \code{cls.draw} has to be provided if 
  \code{loss="VI"}.}
  \item{loss}{the loss function used. Should be one of \code{"Binder"}, \code{"VI"}, or \code{"VI.lb"}. Defaults to \code{"VI.lb"}.}
  \item{start.cl}{clustering used as starting point. If \code{NULL} \code{start.cl= 1:nrow(psm)} is used.}
  \item{maxiter}{integer, maximum number of iterations. Defaults to \code{2*nrow(psm)}.}
  \item{L}{integer, specifies the number of local partitions considered at each iteration. Defaults to \code{2*nrow(psm)}.}
  \item{suppress.comment}{logical, for \code{method="greedy"}, prints a description of the current state (iteration number, number of clusters, posterior expected loss) at each iteration if set to FALSE. Defaults to TRUE.}
}
\details{
This function is called by \code{\link{minVI}} and \code{\link{minbinder.ext}} to optimize the posterior expected loss via a greedy search algorithm. Possible loss functions include Binder's loss (\code{"Binder"}) and the Variation of Information (\code{"VI"}). As computation of the posterior expected Variation of Information is expensive, a third option (\code{"VI.lb"}) is to minimize a modified Variation of Information by swapping the log and expectation. From Jensen's inequality, this can be viewed as minimizing a lower bound to the posterior expected Variation of Information.

At each iteration of the algorithm, we consider the \code{L} closest ancestors or descendants and move in the direction of minimum posterior expected; the distance is measured by Binder's loss or the Variation of Information, depending on the choice of \code{loss}. We recommend trying different starting locations \code{cl.start} and values of \code{l} that control the amount of local exploration. A description of the algorithm at every iteration is printed if \code{suppress.comment=FALSE}. 
}
\value{
\item{cl}{clustering with minimal value of expected loss.}
  \item{value}{value of posterior expected loss.}
  \item{iter.greedy}{the number of iterations the method needed to converge.}
}
\references{
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
\code{\link{minVI}} or \code{\link{minbinder.ext}} which call \code{\link{greedy}} to find the point estimate that minimizes the posterior expected loss.
}
\examples{
data(ex1.data)
x=ex1.data[,c(1,2)]
cls.true=ex1.data$cls.true
plot(x[,1],x[,2],xlab="x1",ylab="x2")
k=max(cls.true)
for(l in 2:k){
points(x[cls.true==l,1],x[cls.true==l,2],col=l)}

# Find representative partition of posterior
data(ex1.draw)
psm=comp.psm(ex1.draw)
ex1.VI=minVI(psm,method=("greedy"),suppress.comment=FALSE)
summary(ex1.VI)
# Different initlization
ex1.VI.v2=minVI(psm,method=("greedy"),suppress.comment=FALSE,start.cl=ex1.draw[nrow(ex1.draw),])
summary(ex1.VI.v2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ cluster }
\keyword{ optimize }% __ONLY ONE__ keyword per line
