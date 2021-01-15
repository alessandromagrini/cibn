\name{getCPT}
\alias{getVariables}
\alias{getTypes}
\alias{getDescription}
\alias{getStates}
\alias{getParSets}
\alias{getCPT}
\title{Functionalities for causal independence Bayesian networks}
\description{Obtain variable names, types, description fields, sample spaces, parent sets and CPTs for a causal independence Bayesian network.}
\usage{getVariables(x)
getTypes(x)
getDescription(x)
getStates(x)
getParSets(x, full=FALSE)
getCPT(x, variables=NULL)
}
\arguments{
  \item{x}{An object of class \code{cibn}.}
\item{full}{Only for function \code{getParSets}: logical value indicating whether the full DAG (i.e., augmented with latent causes and auxiliary nodes implied by the CID) should be considered. If \code{FALSE} (the default), the DAG before the CID is considered.}
\item{variables}{Only for function \code{getCPT}: vector of character strings indicating the name of the variables for which the CPT should be computed. If \code{NULL} (the default), the CPT of all the variables in the Bayesian network will be computed.}
}
%\value{}
\details{Function \code{getVariables} returns the variable names, function \code{getTypes} returns the variable types, function \code{getDescription} returns the description fields, function \code{getStates} returns the sample spaces, function \code{getParSets} returns the parent sets, and function \code{getCPT} computes one or more CPTs of interest.}
%\note{}
%\references{}
\seealso{\link{new.cibn}.}
\examples{
data(bankrisk_code)
bankrisk_bn <- new.cibn(bankrisk_code)
#
getVariables(bankrisk_bn)
getTypes(bankrisk_bn)
getStates(bankrisk_bn)
getDescription(bankrisk_bn)
#
getParSets(bankrisk_bn)
getParSets(bankrisk_bn, full=TRUE)  ## parent sets of the full DAG
#
getCPT(bankrisk_bn)  ## CPTs of all variables
getCPT(bankrisk_bn, variables=c("Portf","Life"))
}
