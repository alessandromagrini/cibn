\name{as.graphNEL}
\alias{as.graphNEL}
\title{Conversion into \code{graphNEL} class}
\description{Convert an object of class \code{cibn} into an object of class \code{graphNEL}.}
\usage{as.graphNEL(x, full=FALSE)}
\arguments{
  \item{x}{An object of class \code{cibn}.}
  \item{full}{Logical value indicating whether the full DAG (i.e., augmented with latent causes and auxiliary nodes implied by the CID) should be considered. If \code{FALSE} (the default), the DAG before the CID is considered.}
}
\value{An object of class \code{graphNEL}.}
%\details{}
%\references{}
\seealso{\link{new.cibn}.}
\examples{
data(bankrisk_code)
bankrisk_bn <- new.cibn(bankrisk_code)
#
G <- as.graphNEL(bankrisk_bn)
G
}