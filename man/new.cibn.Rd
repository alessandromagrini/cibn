\name{new.cibn}
\alias{new.cibn}
\title{Create a causal independence Bayesian network}
\description{A causal independence Bayesian network is created based on prior knowledge on the DAG and on CID parameters.}
\usage{new.cibn(model.code = NULL, path = NULL, maximal = TRUE)}
\arguments{
  \item{model.code}{The model code. See details below.}  
  \item{path}{The path to the model code in text format. See details below.}  
  \item{maximal}{Logical value indicating whether the maximal CID should be applied. Default is \code{TRUE}.} 
  }
\details{For each variable, one command \code{variable} and one command \code{model} must be specified.
Further details will be added soon, for the moment see the examples below.

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
## A simple Bayesian network for inferring risk attitude of bank customers
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
#        with sample space: (money_market, mixed, stock_market)
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

# model code
bankrisk_code <- '

  variable Age {
    type = DGRAD
    states = (18_30, 31_50, 51_)
    parents = ()
    description = <Age>
    }

  model Age {
    omitted = (1,3,2)
    }

  variable Edu {
    type = DGRAD
    states = (primary_or_less, secondary, tertiary)
    parents = ()
    description = <Education level>
    }

  model Edu {
    omitted = (1,7,5)
    }

  variable Marital {
    type = GRAD
    states = (single, convivent)
    parents = (Age)
    description = <Single or convivent>
    }

  model Marital {
    omitted = (2,3)
    Age:18_30 = (3,1)
    Age:51_ = (1,3)
    }

  variable Parent {
    type = GRAD
    states = (no, yes)
    parents = (Age)
    description = <Parentship>
    }

  model Parent {
    omitted = (3,1)
    Age:18_30 = (4,1)
    Age:51_ = (2,1)
    }

  variable Risk {
    type = DGRAD
    states = (low, normal, high)
    parents = (Age, Edu, Marital_Parent)
    description = <Risk attitude>
    }

  interaction Marital_Parent {
    from = (Marital, Parent)
    description = <Interaction between marital status and parentship>
    }

  model Risk {
    omitted = (1,3,2)
    Age:18_30 = (2,1,5)
    Age:51_ = (5,1,2)
    Edu:primary_or_less = (3,2,1)
    Edu:tertiary = (1,3,4)
    Marital_Parent:convivent+no = (1,1,2)
    Marital_Parent:single+yes = (2,1,1)
    Marital_Parent:convivent+yes = (3,1,0)
    }

  variable Portf {
    type = DGRAD
    states = (money_market, mixed, stock_market)
    parents = (Risk)
    description = <Type of portfolio>
    }

  model Portf {
    omitted = (3,5,2)
    Risk:low = (3,1,0)
    Risk:high = (0,1,3)
    }

  variable Life {
    type = NOM
    states = (long_term, short_term, none)
    parents = (Risk)
    description = <Life insurance>
    }

  model Life {
    omitted = (2,4,1)
    Risk:low = (0,1,2)
    Risk:high = (2,1,0)
    }
  '

# create the network
bankrisk_bn <- new.cibn(bankrisk_code)
bankrisk_bn <- new.cibn(bankrisk_code, maximal=FALSE)  ## disable maximal CID

# summary
summary(bankrisk_bn)
}
