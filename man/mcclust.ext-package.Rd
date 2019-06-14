\name{mcclust.ext-package}
\alias{mcclust.ext-package}
\alias{mcclust.ext}
\docType{package}
\title{
Point estimation and credible balls for Bayesian cluster analysis
}
\description{
This is an extension of mcclust package. It provides post-processing tools 
for MCMC samples of partitions to summarize the posterior in Bayesian 
clustering models. Functions for point estimation are provided, giving a 
single representative clustering of the posterior. And, to characterize 
uncertainty in the point estimate, credible balls can be computed.
}
\details{
\tabular{ll}{
Package: \tab mcclust.ext\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-03-24\cr
License: \tab GPL (>= 2)\cr
}
Most important functions:\cr

The functions \code{minVI} and \code{minbinder.ext} find a point estimate of the clustering by minimizing the posterior expected Variation of
Information and Binder's loss, respectively. The function \code{minbinder.ext} extends \code{minbinder} by providing a greedy search optimization 
method to find the optimal clustering. The function \code{minVI} provides several optimization methods to find the optimal clustering. For 
computational reasons, the lower bound to the posterior expected Variation of Information from Jensen's inequality is minimized.

The function \code{credibleball} computes a credible ball around the clustering estimate to characterize uncertainty. It returns the upper vertical,
lower vertical, and horizontal bounds to describe the credible ball.

The function \code{plotpsm} produces a heat map of the posterior similarity matrix.
}

\author{
Sara Wade

Maintainer: Sara Wade <sara.wade@eng.cam.ac.uk>
}
\references{

Binder, D.A. (1978) Bayesian cluster analysis, \emph{Biometrika} \bold{65}, 31--38.

Fritsch, A. and Ickstadt, K. (2009) An improved criterion for clustering based on the
posterior similarity matrix, \emph{Bayesian Analysis}, \bold{4},367--391.

Lau, J.W. and Green, P.J. (2007) Comparing clusters--an information based distance
procedures, \emph{Journal of Computational and Graphical Statistics} \bold{16}, 526--558.

Meila, M. (2007) Bayesian model based clustering
procedures, \emph{Journal of Multivariate Analysis} \bold{98}, 873--895.

Wade, S. and Ghahramani, Z. (2015) Bayesian cluster analysis: Point estimation and credible balls. 
Submitted. arXiv:1505.03339.
}
\keyword{ package }
\seealso{
\code{\link[mcclust:mcclust-package]{mcclust}}
}
\examples{
data(galaxy.fit)
x=data.frame(x=galaxy.fit$x)
data(galaxy.pred)
data(galaxy.draw)

# Find representative partition of posterior
# Variation of Information (minimizes lower bound to VI)
psm=comp.psm(galaxy.draw)
galaxy.VI=minVI(psm,galaxy.draw,method=("all"),include.greedy=TRUE)
summary(galaxy.VI)
plot(galaxy.VI,data=x,dx=galaxy.fit$fx,xgrid=galaxy.pred$x,dxgrid=galaxy.pred$fx)
# Compute Variation of Information
VI(galaxy.VI$cl,galaxy.draw)
# Binder
galaxy.B=minbinder.ext(psm,galaxy.draw,method=("all"),include.greedy=TRUE)
summary(galaxy.B)
plot(galaxy.B,data=x,dx=galaxy.fit$fx,xgrid=galaxy.pred$x,dxgrid=galaxy.pred$fx)

# Uncertainty in partition estimate
galaxy.cb=credibleball(galaxy.VI$cl[1,],galaxy.draw)
summary(galaxy.cb)
plot(galaxy.cb,data=x,dx=galaxy.fit$fx,xgrid=galaxy.pred$x,dxgrid=galaxy.pred$fx)

# Compare with uncertainty in heat map of posterior similarity matrix
plotpsm(psm)
}
