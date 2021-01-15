\name{query.cibn}
\alias{query.cibn}
\title{Inference in a causal independence Bayesian network}
\description{Perform exact inference in a causal independence Bayesian network through the joint tree algorithm (interface to the \code{gRain} package)}
\usage{query.cibn(x, target=NULL, evidence=NULL, type="marginal")}
\arguments{
  \item{x}{An object of class \code{cibn}.}
  \item{target}{The name of the target variable. If \code{NULL} (the default), all the variables in the Bayesian network will be set as targets.}
  \item{evidence}{A named list with each component indicating the evidence on a specific variable in the form of a vector of state names. See the examples below.}
  \item{type}{A character string indicating the type of inference: \code{"marginal"}, \code{"joint"} or \code{"conditional"}. Default is \code{"marginal"}.}  
}
\value{A list with one component for each variable in the Bayesian network, indicating its probability distribution given the evidence.}
\details{This function is an interface to function \code{querygrain} in the \code{gRain} package.}
\references{
S. L. Lauritzen and D. J. Spiegelhalter (1988). Local computations with probabilities on graphical structures and their application to expert systems. \emph{Journal of the Royal Statistical Society, Series B (Methodological)}, 50(2): 157-224. DOI: 10.1023/A:1008935617754.
}
\seealso{\link{new.cibn}; \link{sample.cibn}.}
\examples{
data(bankrisk_code)
bankrisk_bn <- new.cibn(bankrisk_code)
#
getStates(bankrisk_bn)  ## see the sample spaces
#
query.cibn(bankrisk_bn, target="Risk", evidence=list(Age="31_50",Portf="mixed"))
query.cibn(bankrisk_bn, target="Risk", evidence=list(Age="31_50",Portf="money_market"))
query.cibn(bankrisk_bn, target="Risk", evidence=list(Age="31_50",Portf="stock_market"))
}