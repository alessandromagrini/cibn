\name{sample.cibn}
\alias{sample.cibn}
\title{Sampling from a causal independence Bayesian network}
\description{Draw a random sample from a causal independence Bayesian network}
\usage{sample.cibn(x, nsam, seed=NULL)}
\arguments{
  \item{x}{An object of class \code{cibn}.}
  \item{nsam}{The number of sample units.}
  \item{seed}{The seed for the random number generator. If \code{NULL} (the default), it is chosen randomly.}
}
\value{An object of class \code{data.frame}.}
%\details{This function is an interface to function \code{simulate} applied to the object of class \code{grain}.}
%\references{}
\seealso{\link{new.cibn}; \link{query.cibn}.}
\examples{
data(bankrisk_code)
bankrisk_bn <- new.cibn(bankrisk_code)
#
bankrisk_sam <- sample.cibn(bankrisk_bn, nsam=100)
head(bankrisk_sam)
summary(bankrisk_sam)
}