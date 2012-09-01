\name{fixef.lmekin}
\alias{fixef.lmekin}
\alias{ranef.lmekin}
\alias{vcov.lmekin}
\alias{VarCorr.lmekin}
\alias{logLik.lmekin}
\title{Extraction functions for Lmekin}
\description{Extract the fixed effects, random effects, variance of the
  fixed effects, or variance of the random effects from a linear mixed
  effects model fit with lmekin.}
\usage{
\method{fixef}{lmekin}(object, \dots)
\method{ranef}{lmekin}(object, \dots)
\method{vcov}{lmekin}(object, \dots)
\method{VarCorr}{lmekin}(x, \dots)
\method{logLik}{lmekin}(object, \dots)
}
\arguments{
  \item{object}{an object inheriting from class \code{lmekin} 
     representing the result of a mixed effects model.
   }
   \item{x}{an object inheriting from class \code{lmekin} 
     representing the result of a mixed effects model.
   }
 \item{\dots}{some methods for this generic require additional
   arguments.  None are used in this method.} 
}
\value{the fixed effects are a vector and vcov returns their
  variance/covariance matrix.  The random effects are a list with one
  element for each random effect.  The ranef component contains the
  coefficients and VarCorr the estimated variance/covariance matrix.
  The logLik method returns the loglikelihood along with its
  degrees of freedom.}

\details{
  For the random effects model \eqn{y = X\beta + Zb + \epsilon}{y = X beta
    + Zb + e}, let \eqn{\sigma^2}{s^2} be the variance of the error term
  \eqn{\epsilon}{e}.
  Let \eqn{A= \sigma^2 P}{A = sP} be the variance of the random effects
  \eqn{b}{b}.  There is a computational advantage to solving the problem
  in terms of \eqn{P}{P} instead of \eqn{A}{A}, and that is what is
  stored in the returned lmekin object.
  The \code{VarCorr} function returns elements of \eqn{P}{P}; the print
  and summary functions report values of \eqn{A}{A}.
  Pinhiero and Bates call \eqn{P}{P} the precision factor.
  }
\author{Terry Therneau}
\seealso{\code{\link{lmekin}}, \code{\link{random.effects}},
  \code{\link{fixed.effects}}, \code{link{vcov}}, \code{\link{VarCorr}}
}
\references{
  J Pinheiro and D Bates, Mixed-effects models in S and
  S-Plus. Springer, 2000.
  }
  \examples{
data(ergoStool, package="nlme")  # use a data set from nlme
efit <-  lmekin(effort ~ Type + (1|Subject), ergoStool)
ranef(efit)
}
\keyword{models}

