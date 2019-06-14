\name{credibleball}
\alias{credibleball}
\alias{summary.credibleball}
\alias{print.summary.credibleball}
\alias{plot.credibleball}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute a Bayesian credible ball around a clustering estimate
}
\description{
Computes a Bayesian credible ball around a clustering estimate to characterize uncertainty in the posterior, i.e. MCMC samples of clusterings.
}
\usage{
credibleball(c.star, cls.draw, c.dist = c("VI","Binder"), alpha = 0.05)

\method{summary}{credibleball}(object, ...)
\method{plot}{credibleball}(x,data=NULL,dx=NULL,xgrid=NULL,dxgrid=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{c.star}{vector, a clustering estimate of the \code{length(c.star)} data points.
}
  \item{cls.draw}{a matrix of the MCMC samples of clusterings of the \code{ncol(cls.draw)} data points.
}
  \item{c.dist}{the distance function on clusterings to use. Should be one of \code{"VI"} or \code{"Binder"}. Defaults to \code{"VI"}. }
  \item{alpha}{a number in the unit interval, specifies the Bayesian confidence level of \code{1-alpha}. Defaults to \code{0.05}.}
  \item{object}{
      an object of class \code{"credibleball"}.
}
  \item{x}{
      an object of class \code{"credibleball"}.}
\item{data}{the dataset contained in a \code{data.frame} with \code{ncol(cls.draw)} rows of data points.}
\item{dx}{ for \code{ncol(x)}=1, the estimated density at the observed data points.}
\item{xgrid}{for \code{ncol(x)}=1, a grid of data points for density estimation.}
\item{dxgrid}{ for \code{ncol(x)}=1, the estimated density at the grid of data points.}
\item{\dots}{other inputs to \code{summary} or \code{plot}.
}
}
\details{
An advantage of Bayesian cluster analysis is that it provides a posterior over the entire partition space, expressing beliefs in the clustering structure given the data. The credible ball summarizes the uncertainty in the posterior around a clustering estimate \code{c.star} and is defined as the smallest ball around \code{c.star} with posterior probability at least \code{1-alpha}. Possible distance metrics on the partition space are the Variation of Information and the N-invariant Binder's loss (Binder's loss times \code{2/length(c.star)^2}). The posterior probability is estimated from MCMC posterior samples of clusterings.

The credible ball is summarized via the upper vertical, lower vertical, and horizontal bounds, defined, respectively, as the partitions in the credible ball with the fewest clusters that are most distant to \code{c.star}, with the most clusters that are most distant to \code{c.star}, and with the greatest distance to \code{c.star}.

In plots, data points are colored according to cluster membership. For \code{nrow(data)=1}, the data points are plotted against the density (which is estimated via a call to \code{\link{density}} if not provided). For \code{nrow(data)=2} the data points are plotted, and for \code{nrow(data)>2}, the data points are plotted in the space spanned by the first two principal components.
}
\value{
  \item{c.star}{vector, clustering estimate of the \code{length(c.star)} data points.
}
	\item{c.horiz}{A matrix of horizontal bounds of the credible ball, i.e. partitions in the credible ball with the greatest distant to \code{c.star}.
}
	\item{c.uppervert}{A matrix of upper vertical bounds of the credible ball, i.e. partitions in the credible ball with the fewest clusters that are most distant to \code{c.star}.
}
	\item{c.lowervert}{A matrix of lower vertical bounds of the credible ball, i.e. partitions in the credible ball with the most clusters that are most distant to \code{c.star}.
}
	\item{dist.horiz}{the distance between \code{c.star} and the horizontal bounds
}
	\item{dist.uppervert}{the distance between \code{c.star} and the upper vertical bounds
}
	\item{dist.lowervert}{the distance between \code{c.star} and the lower vertical bounds
}
}
\references{
Wade, S. and Ghahramani, Z. (2015) Bayesian cluster analysis: Point estimation and credible balls. 
Submitted. arXiv:1505.03339.
}
\author{
Sara Wade, \email{sara.wade@eng.cam.ac.uk}
}
\seealso{\code{\link{minVI}}, \code{\link{minbinder.ext}}, \code{\link{maxpear}}, and \code{\link{medv}}
to obtain a point estimate of clustering based on posterior MCMC samples; and \code{\link{plotpsm}} for a heat map of posterior similarity matrix.
}
\examples{
data(galaxy.fit)
x=data.frame(x=galaxy.fit$x)
data(galaxy.pred)
data(galaxy.draw)

# Find representative partition of posterior
psm=comp.psm(galaxy.draw)
galaxy.VI=minVI(psm,galaxy.draw,method=("all"),include.greedy=TRUE)
summary(galaxy.VI)
plot(galaxy.VI,data=x,dx=galaxy.fit$fx,xgrid=galaxy.pred$x,dxgrid=galaxy.pred$fx)

# Uncertainty in partition estimate
galaxy.cb=credibleball(galaxy.VI$cl[1,],galaxy.draw)
summary(galaxy.cb)
plot(galaxy.cb,data=x,dx=galaxy.fit$fx,xgrid=galaxy.pred$x,dxgrid=galaxy.pred$fx)

# Compare with heat map of posterior similarity matrix
plotpsm(psm)
}
\keyword{ cluster}
