\name{new.cibn}
\alias{new.cibn}
\title{Create a causal independence Bayesian network}
\description{Create a causal independence Bayesian network based on prior knowledge on the DAG and on CID parameters.}
\usage{new.cibn(model.code = NULL, path = NULL, maximal = TRUE)}
\arguments{
  \item{model.code}{The model code. See details below.}  
  \item{path}{The path to the model code in text format. See details below.}  
  \item{maximal}{Logical value indicating whether the maximal CID should be applied. Default is \code{TRUE}.} 
  }
\details{For each variable, one command \code{variable} and one command \code{model} must be specified.
Further details will be added soon, see the examples below.

The name of a variable must begin with a capital letter and cannot include special characters excepting '_'.
The name of a state cannot include special characters excepting '_', and cannot begin with 'LAMBDA' or 'AUX'.

S3 methods \code{print}, \code{summary} and \code{plot} are available for class \code{cibn}.
}
%\note{}
\references{
A. Magrini (2021). Efficient decomposition of Bayesian networks with non-graded variables. To be appeared on \emph{International Journal of Statistics and Probability}, 10(2).
}
\seealso{\link{query.cibn}; \link{sample.cibn}.}
\examples{
## A simple Bayesian network to infer risk attitude of bank customers
#
#   Variables:
#    - 'Age': age in years, double-graded variable
#        with sample space: (18_30, 31_50, 51_);
#    - 'Edu': education level, double-graded variable
#        with sample space: (primary_or_less, secondary, tertiary);
#    - 'Marital': marital status, graded variable
#        with sample space: (single, convivent);
#    - 'Parent': parentship, graded variable
#        with sample space: (no, yes);
#    - 'Risk': risk attitude, double-graded variable
#        with sample space: (low, normal, high);
#    - 'Portf': type of portfolio, double-graded variable
#        with sample space: (money_market, mixed, stock_market);
#    - 'Life': life insurance, multi-valued nominal variable
#        with sample space: (long_term, short_term, none).
#
#   Edges in the DAG:
#    - 'Age' -> 'Marital'
#    - 'Age' -> 'Parent'
#    - 'Age' -> 'Risk'
#    - 'Edu' -> 'Risk'
#    - 'Marital' -> 'Risk'
#    - 'Parent' -> 'Risk'
#    - 'Risk' -> 'Portf'
#    - 'Risk' -> 'Life'
#
#   Causal interactions:
#    - between 'Marital' and 'Parent' in determining 'Risk'
#

# load model code
data(bankrisk_code)

# create the network
bankrisk_bn <- new.cibn(bankrisk_code)
bankrisk_bn <- new.cibn(bankrisk_code, maximal=FALSE)  ## disable maximal CID

# summary
summary(bankrisk_bn)
}
