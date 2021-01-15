\name{cibn-package}
\alias{cibn-package}
\docType{package}
\title{Causal Independence Bayesian Networks}
\description{Elicitation, estimation and inference functionalities for Bayesian networks under the causal independence assumption.}
\details{
  \tabular{ll}{
  Package: \tab cibn\cr
  Type: \tab Package\cr
  Version: \tab 0.0\cr
  Date: \tab 2021-01-07\cr
  License: \tab GPL-2
  }
Causal independence Bayesian networks (Magrini, 2021) are Bayesian networks with non-interacting parent variables (causal independence assumption). They allow three exaustive types of variables (graded, double-graded and multi-valued nominal variables) and admit the Causal Independence Decomposition (CID), which increases efficiency of elicitation, estimation and inference. Ccausal interactions can be added upon need.
The main functions of the package are:
\itemize{
  \item{\link{new.cibn}, }{to create a new network based on prior knowledge;}
  \item{\code{update.cibn}, }{to update an existing network based on possibly incomplete data (not still implemented but available soon);}
  \item{\link{query.cibn}, }{to perform exact inference in a network through an interface to the \code{gRain} package;}
  \item{\link{sample.cibn}, } {to draw a random sample from a network.}
  }
Also, method \link{plot} allows to display the DAG of a network exploiting the \code{graphNEL} class.
}
\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
\references{
A. Magrini (2021). Efficient decomposition of Bayesian networks with non-graded variables. To be appeared on \emph{International Journal of Statistics and Probability}, 10(2).
}
