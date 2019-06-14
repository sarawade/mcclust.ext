\name{minVI}
\alias{minVI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Minimize the posterior expected Variation of Information
}
\description{
Finds a representative partition of the posterior by minimizing the lower bound to the posterior
expected Variation of Information from Jensen's Inequality.
}
\usage{
minVI(psm, cls.draw=NULL, method=c("avg","comp","draws","greedy","all"), 
      max.k=NULL, include.greedy=FALSE, start.cl=NULL, maxiter=NULL,
      l=NULL, suppress.comment=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{psm}{a posterior similarity matrix, which can be obtained from MCMC samples of clusterings through a call to \code{comp.psm}.}
 \item{cls.draw}{a matrix of the MCMC samples of clusterings of the \code{ncol(cls)} data points that have been used to compute \code{psm}. Note: \code{cls.draw} has to be provided if 
  \code{method="draw"} or \code{"all"}.}
 \item{method}{the optimization method used. Should be one of \code{"avg"}, \code{"comp"}, \code{"draws"}, \code{"greedy"} or \code{"all"}. Defaults to \code{"avg"}.}
 \item{max.k}{ integer, if \code{method="avg"} or \code{"comp"} the maximum number of clusters up to which the hierarchical clustering is cut.
  Defaults to \code{ceiling(nrow(psm)/4)}. }
  \item{include.greedy}{logical, should method \code{"greedy"} be included when \code{method="all"}? Defaults to FALSE.}
  \item{start.cl}{clustering used as starting point for \code{method="greedy"}. If \code{NULL} \code{start.cl= 1:nrow(psm)} is used.}
  \item{maxiter}{integer, maximum number of iterations for \code{method="greedy"}. Defaults to \code{2*nrow(psm)}.}
  \item{l}{integer, specifies the number of local partitions considered at each iteration for \code{method="greedy"}. Defaults to \code{2*nrow(psm)}.}
  \item{suppress.comment}{logical, for \code{method="greedy"}, prints a description of the current state (iteration number, number of clusters, posterior expected loss) at each iteration if set to FALSE. Defaults to TRUE.}
% \item{loss}{the loss function used. Should be one of \code{"Binder"}, \code{"VI"}, or \code{"VI.lb"}. Defaults to \code{"VI.lb"}. Is set equal to \code{"VI.lb"} when \code{greedy} is called from \code{minVI}.}
}
\details{
The Variation of Information between two clusterings is defined as the sum of the entropies minus two times the mutual information. Computation of the posterior expected Variation
of Information can be expensive, as it requires a Monte Carlo estimate. We consider a modified posterior expected Variation of Information, obtained by swapping the log and expectation, which is much more computationally efficient as it only depends on the posterior through the posterior similarity matrix. From Jensen's inequality, the problem can be viewed as minimizing a lower bound to the posterior expected loss. 

We provide several optimization methods. For \code{method="avg"} and \code{"comp"}, the search is restricted to the clusterings obtained from a hierarchical clustering with average/complete linkage and \code{1-psm} as a distance matrix (the clusterings with number of clusters \code{1:max.k} are considered).\cr
 Method \code{"draws"} restricts the search to the clusterings sampled in the MCMC algorithm. \cr
 Method \code{"greedy"} implements a greedy search algorithm, where at each iteration, we consider the \code{l} closest ancestors or descendants and move in the direction of minimum posterior expected loss with the VI distance. We recommend trying different starting locations \code{cl.start} and values of \code{l} that control the amount of local exploration. Depending on the starting location and \code{l}, the method can take some time to converge, thus it is only included in \code{method="all"} if \code{include.greedy=TRUE}. If \code{method="all"}, the starting location \code{cl.start} defaults to the best clustering found by the other methods.  A description of the algorithm at every iteration is printed if \code{suppress.comment=FALSE}. 
If \code{method="all"} all minimization methods except \code{"greedy"} are applied by default. 
}
\value{
 \item{cl}{clustering with minimal value of expected loss. If \code{method="all"} a matrix containing the clustering with the smallest value
   of the expected loss over all methods in the first row and the clusterings of the individual methods in the next rows.}
  \item{value}{value of posterior expected loss. A vector corresponding to the rows of \code{cl} if \code{method="all"}.}
  \item{method}{the optimization method used.}
  \item{iter.greedy}{if \code{method="greedy"} or \code{method="all"} and \code{include.greedy=T} the number of iterations the method needed to converge.}
}
\references{
Meila, M. (2007) Bayesian model based clustering
procedures, \emph{Journal of Multivariate Analysis} \bold{98}, 873--895.

Wade, S. and Ghahramani, Z. (2015) Bayesian cluster analysis: Point estimation and credible balls. 
Submitted. arXiv:1505.03339
}
\author{
Sara Wade, \email{sara.wade@eng.cam.ac.uk}
}
\seealso{\code{\link{summary.c.estimate}} and \code{\link{plot.c.estimate}} to summarize and plot the resulting output from \code{\link{minVI}} or \code{\link{minbinder.ext}}; \code{\link{VI}} or \code{\link{VI.lb}} for computing the posterior expected Variation of Information or the modified version from swapping the log and expectation;  \code{\link{comp.psm}} for computing posterior similarity matrix; \code{\link{maxpear}}, \code{\link{minbinder.ext}}, and \code{\link{medv}}
for other point estimates of clustering based on posterior; and \code{\link{credibleball}} to compute credible ball characterizing uncertainty around
the point estimate.
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
plot(ex2.VI,data=x)

# Compare with Binder
ex2.B=minbinder.ext(psm,ex2.draw,method=("all"),include.greedy=TRUE)
summary(ex2.B)
plot(ex2.B,data=x)
}
\keyword{ cluster}
\keyword{ optimize}% __ONLY ONE__ keyword per line
