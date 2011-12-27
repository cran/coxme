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
  \item{vvar}{variance matrix for vcoef}
  \item{residuals}{vector of residuals from the fit}
  \item{method}{either "ML" or "REML"}
  \item{loglik}{the log-likelihood for the fitted model}
  \item{sigma}{the estimated residual error}
  \item{n}{number of observations used}
  \item{call}{a copy of the call}
  \item{na.action}{this will be present if any observations were removed
    due to missing values}
}
\details{
  The random effects are modeled internally as \eqn{N(0, \sigma^2
    A)}{N(0, s^2A)} where \eqn{\sigma^2}{s^2} is the residual variance.
  The parameters of \code{vcoef} refer to \eqn{A}.  THe printed results
  for the model have been multiplied by \eqn{\sigma}{s}, removing the
  factorization.

  The variance matrix \code{vvar} is for the transformed parameters.
  For the details of the transformation see the individual variance
  functions. The \code{coxmeFull} routine, for instance, uses a log
  transform for variances.  Use this with caution.}
\author{Terry Therneau}
\seealso{\code{\link{lmekin}}, \code{\link{coxmeFull}},
  \code{\link{coxmeMlist}}}
\keyword{models}

