\name{lmekin.object}
\alias{lmekin.object}
\title{lmekin object}
\description{
This class of object is returned by the \code{lmekin} function to 
represent a fittd mixed effect linear model.
Objects of this class currently have methods for \code{print} and
\code{residuals}.
}
\value{ A list with the folling components:
  \item{coefficients}{a list with components \code{fixed} and
    \code{random}; the first will be NULL for a model with no fixed
    effects.  The \code{random} component is itself a list, with an
    element for each random effect.}
  \item{var}{the variance matrix of the fixed effects}
  \item{vcoef}{the parameters of the variance matrix of the random
    effects.}
  \item{residuals}{vector of residuals from the fit}
  \item{method}{either "ML" or "REML"}
  \item{loglik}{the log-likelihood for the fitted model}
  \item{sigma}{the estimated residual error}
  \item{n}{number of observations used}
  \item{call}{a copy of the call}
  \item{na.action}{this will be present if any observations were removed
    due to missing values}
}
\author{Terry Therneau}
\seealso{\code{\link{lmekin}}, \code{\link{coxmeFull}},
  \code{\link{coxmeMlist}}}
\keyword{models}

