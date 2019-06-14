\name{plotpsm}
\alias{plotpsm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot a heat map of the posterior similarity matrix
}
\description{
Produces a heat map of the posterior similarity matrix with data points reordered by hierarchical clustering.
}
\usage{
plotpsm(psm, method = "complete", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{psm}{a posterior similarity matrix, which can be obtained from MCMC samples of clusterings through a call to \code{comp.psm}.}
  \item{method}{the agglomeration method to be used in hierarchical clustering. Defaults to \code{"complete"}. See \code{hclust}.
}
  \item{\dots}{other inputs to \code{image}.
}
}
\details{
Produces a heatmap of the posterior similarity matrix with red representing high posterior probability of one and white representing low posterior probability of zero. Data points are first reordered by hierarchical clustering to increasing legibility.
}
\value{
Produces a heatmap of the posterior similarity matrix.
%% ...
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Sara Wade, \email{sara.wade@eng.cam.ac.uk}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{comp.psm}} for computing posterior similarity matrix; \code{\link{hclust}} for hierarchical clustering; and \code{\link{credibleball}} for an alternative representation of uncertainty in the posterior on clusterings.
}
\examples{
data(ex1.data)
x=ex1.data[,c(1,2)]
cls.true=ex1.data$cls.true
plot(x[,1],x[,2],xlab="x1",ylab="x2")
k=max(cls.true)
for(l in 2:k){
points(x[cls.true==l,1],x[cls.true==l,2],col=l)}

# Heat map to represent posterior uncertainty
data(ex1.draw)
psm=comp.psm(ex1.draw)
plotpsm(psm)

}
\keyword{ cluster}

