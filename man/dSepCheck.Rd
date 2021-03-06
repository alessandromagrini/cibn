\name{dSepCheck}
\alias{dSepCheck}
\title{Conditional independence check}
\description{Check conditional independence between two variables based on the d-separation criterion.}
\usage{dSepCheck(x, var1, var2, given = NULL)}
\arguments{
  \item{x}{An object of class \code{cibn}.}
  \item{var1}{The name of the first variable.}
  \item{var2}{The name of the second variable.}
  \item{given}{A vector containing the names of conditioning variables. If \code{NULL}, marginal independence is checked.}
}
\value{Logical}
\details{The d-separation is a necessary and sufficient condition for conditional independence. See Pearl (2000), page 16 and following).}
\note{The result is unchanged if arguments \code{var1} and \code{var2} are switched.

Dependence is a necessary but not sufficient condition for causation: see the discussion in Pearl (2000).}
\references{J. Pearl (2000). Causality: models, reasoning, and inference. Cambridge University Press. Cambridge, UK. ISBN: 978-0-521-89560-6}
\seealso{\link{new.cibn}.}
\examples{
data(bankrisk_code)
bankrisk_bn <- new.cibn(bankrisk_code)
#
dSepCheck(bankrisk_bn, var1="Age", var2="Edu")
dSepCheck(bankrisk_bn, var1="Portf", var2="Edu", given="Risk")
dSepCheck(bankrisk_bn, var1="Portf", var2="Edu", given=c("Risk","Life"))
}