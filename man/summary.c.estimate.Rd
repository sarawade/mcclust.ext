\name{summary.c.estimate}
%\alias{c.estimate}
\alias{summary.c.estimate}
\alias{print.summary.c.estimate}
\alias{plot.c.estimate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarize and plot the estimate of the partition
}
\description{
Summarizes and plots the estimate of the partition in a Bayesian cluster analysis model.
}
\usage{
\method{summary}{c.estimate}(object, ...)
\method{plot}{c.estimate}(x,data=NULL,dx=NULL,xgrid=NULL,dxgrid=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
      an object of class \code{"c.estimate"}, i.e. a clustering estimate .
}
  \item{x}{
      an object of class \code{"c.estimate"}, i.e. a clustering estimate .
}
\item{data}{the dataset contained in a \code{data.frame} with \code{ncol(x$cl)} rows of data points.}
\item{dx}{ for \code{ncol(data)}=1, the estimated density at the observed data points.}
\item{xgrid}{for \code{ncol(data)}=1, a grid of data points for density estimation.}
\item{dxgrid}{ for \code{ncol(data)}=1, the estimated density at the grid of data points.}
\item{\dots}{other inputs to \code{summary} or \code{plot}.
}
}
\details{
Summarizes and plots the clustering estimates returned by the functions \code{minVI} and \code{minbinder.ext}. In plots, data points are colored according to cluster membership. For \code{nrow(x)=1}, the data points are plotted against the density (which is estimated via a call to \code{\link{density}} if not provided). For \code{nrow(x)=2} the data points are plotted, and for \code{nrow(x)>2}, the data points are plotted in the space spanned by the first two principal components.
}
\value{
	\item{method}{
		the optimization method used to obtain the point estimate.
}
	\item{k}{
		(a vector of) the number of clusters in the point estimate. Returns a vector if \code{n.c>0}.
}
	\item{n.c}{
		the number of point estimates in the object.
}
	\item{t}{
		a list of length \code{n.c} of the table(s) with cluster sizes.
}
	\item{value}{
		(a vector of) the posterior expected loss. Returns a vector if \code{n.c>0}.
}
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
\code{\link{minVI}} and \code{\link{minbinder.ext}}
}
\examples{
data(galaxy.draw)
data(galaxy.fit)
data(galaxy.pred)
x=data.frame(x=galaxy.fit$x)

# Find representative partition of posterior
psm=comp.psm(galaxy.draw)
galaxy.VI=minVI(psm,galaxy.draw,method=("all"),include.greedy=TRUE)
summary(galaxy.VI)
plot(galaxy.VI,data=x,dx=galaxy.fit$fx,xgrid=galaxy.pred$x,dxgrid=galaxy.pred$fx)
galaxy.B=minbinder.ext(psm,galaxy.draw,method=("all"),include.greedy=TRUE)
summary(galaxy.B)
plot(galaxy.B,data=x,dx=galaxy.fit$fx,xgrid=galaxy.pred$x,dxgrid=galaxy.pred$fx)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{cluster}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
