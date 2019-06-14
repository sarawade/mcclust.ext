\name{minbinder.ext}

\alias{minbinder.ext}

%- Also NEED an '\alias' for EACH other topic documented here.

\title{Minimize the posterior expected Binder's loss
}

\description{

Finds a representative partition of the posterior by minimizing the posterior expected Binder's loss.
}

\usage{


minbinder.ext(psm, cls.draw = NULL, method = c("avg", "comp", "draws", "laugreen", "greedy", "all"), 
	      max.k = NULL, include.lg = FALSE, include.greedy = FALSE, start.cl.lg = NULL, 
              start.cl.greedy = NULL, tol = 0.001, maxiter = NULL, l = NULL, suppress.comment = TRUE)

}

%- maybe also 'usage' for other objects documented here.

\arguments{
  
\item{psm}{a posterior similarity matrix, which can be obtained from MCMC samples of clusterings through a call to \code{comp.psm}.}
 \item{cls.draw}{a matrix of the MCMC samples of clusterings of the \code{ncol(cls)} data points that have been used to compute \code{psm}. Note: \code{cls.draw} has to be provided if 
  \code{method="draw"} or \code{"all"}.}
 \item{method}{the optimization method used. Should be one of \code{"avg"}, \code{"comp"}, \code{"draws"}, \code{"laugreen"}, \code{"greedy"} or \code{"all"}. Defaults to \code{"avg"}.}
 \item{max.k}{ integer, if \code{method="avg"} or \code{"comp"} the maximum number of clusters up to which the hierarchical clustering is cut.
  Defaults to \code{ceiling(nrow(psm)/4)}. }
  \item{include.lg}{logical, should method \code{"laugreen"} be included when \code{method="all"}? Defaults to FALSE.}
  \item{include.greedy}{logical, should method \code{"greedy"} be included when \code{method="all"}? Defaults to FALSE.}
  \item{start.cl.lg}{clustering used as starting point for \code{method="laugreen"}. If \code{NULL} \code{start.cl= 1:nrow(psm)} is used.}
  \item{start.cl.greedy}{clustering used as starting point for \code{method="greedy"}. If \code{NULL} \code{start.cl= 1:nrow(psm)} is used.}
    \item{tol}{convergence tolerance for \code{method="laugreen"}.}
  \item{maxiter}{integer, maximum number of iterations for \code{method="greedy"}. Defaults to \code{2*nrow(psm)}.}
  \item{l}{integer, specifies the number of local partitions considered at each iteration for \code{method="greedy"}. Defaults to \code{2*nrow(psm)}.}
  \item{suppress.comment}{logical, for \code{method="greedy"}, prints a description of the current state (iteration number, number of clusters, posterior expected loss) at each iteration if set to FALSE. Defaults to TRUE.}
}

\details{
This functions extends \code{\link{minbinder}} by implementing the greedy search algorithm to minimize the posterior expected Binder's loss. 

Binder's loss counts the number of disagreements in all possible pairs of data points. The value returned is the posterior expected N-invariant Binder's loss, which is defined by multiplying Binder's loss times 2 and dividing by N^2, N representing the sample size, and is so-called because it only depends on the sample size through the proportion of data points in each cluster intersection.

The function \code{\link{minbinder}} is called for optimization methods \code{method="avg"}, \code{"comp"}, \code{method="draws"}, and \code{"laugreen"}. \cr
 Method \code{"greedy"} implements a greedy search algorithm, where at each iteration, we consider the \code{l} closest ancestors or descendants and move in the direction of minimum posterior expected loss with the N-invariant Binder's loss as the distance. We recommend trying different starting locations \code{cl.start} and values of \code{l} that control the amount of local exploration. Depending on the starting location and \code{l}, the method can take some time to converge, thus it is only included in \code{method="all"} if \code{include.greedy=TRUE}. If \code{method="all"}, the starting location \code{cl.start} defaults to the best clustering found by the other methods. A description of the algorithm at every iteration is printed if \code{suppress.comment=FALSE}. 
If \code{method="all"} all minimization methods except \code{"laugreen"} and \code{"greedy"} are applied by default. 

}

\value{
 \item{cl}{clustering with minimal value of expected loss. If \code{method="all"} a matrix containing the clustering with the smallest value
   of the expected loss over all methods in the first row and the clusterings of the individual methods in the next rows.}
  \item{value}{value of posterior expected loss. A vector corresponding to the rows of \code{cl} if \code{method="all"}.}
  \item{method}{the optimization method used.}
  \item{iter.greedy}{if \code{method="greedy"} or \code{method="all"} and \code{include.greedy=T} the number of iterations the method needed to converge.}
    \item{iter.lg}{if \code{method="laugreen"} or \code{method="all"} and \code{include.lg=T} the number of iterations the method needed to converge.}
}

\references{
Binder, D.A. (1978) Bayesian cluster analysis, \emph{Biometrika} \bold{65}, 31--38.

Fritsch, A. and Ickstadt, K. (2009) An improved criterion for clustering based on the
posterior similarity matrix, \emph{Bayesian Analysis}, \bold{4},367--391.

Lau, J.W. and Green, P.J. (2007) Bayesian model based clustering
procedures, \emph{Journal of Computational and Graphical Statistics} \bold{16}, 526--558.

Wade, S. and Ghahramani, Z. (2015) Bayesian cluster analysis: Point estimation and credible balls. Sumbitted. arXiv:1505.03339.
}

\author{
Sara Wade, \email{sara.wade@eng.cam.ac.uk}
}

%\note{
%%  ~~further notes~~
%}


%% ~Make other sections like Warning with \section{Warning }{....} 

\seealso{
\code{\link{summary.c.estimate}} and \code{\link{plot.c.estimate}} to summarize and plot the resulting output from \code{\link{minVI}} or \code{\link{minbinder.ext}}; \code{\link{comp.psm}} for computing posterior similarity matrix; \code{\link{maxpear}}, \code{\link{minVI}}, and \code{\link{medv}}
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
ex2.B=minbinder.ext(psm,ex2.draw,method=("all"),include.greedy=TRUE)
summary(ex2.B)
plot(ex2.B,data=x)

# Compare with VI
ex2.VI=minVI(psm,ex2.draw,method=("all"),include.greedy=TRUE)
summary(ex2.VI)
plot(ex2.VI,data=x)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the

% R documentation directory.

\keyword{cluster}

\keyword{optimize}% __ONLY ONE__ keyword per line
