\name{plot.cibn}
\alias{plot.cibn}
\title{Graphic for the DAG of a causal independence Bayesian network}
\description{Obtain the graphic for the DAG of a causal independence Bayesian network.}
\usage{
\method{plot}{cibn}(x, full=FALSE, ...)
}
\arguments{
  \item{x}{An object of class \code{cibn}.}
  \item{full}{Logical value indicating whether the full DAG (i.e., augmented with latent causes and auxiliary nodes implied by the CID) should be displayed. If \code{FALSE} (the default), the DAG before the CID is displayed.}
  \item{...}{Further graphical parameters.}
}
%\value{}
%\details{}
%\note{}
%\references{}
\seealso{\link{new.cibn}.}
\examples{
data(bankrisk_code)
bankrisk_bn <- new.cibn(bankrisk_code)
#
plot(bankrisk_bn, attrs=list(edge=list(arrowsize=0.5)))
plot(bankrisk_bn, attrs=list(edge=list(arrowsize=0.5)), full=TRUE)  ## full DAG
}
