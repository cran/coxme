\name{fixef.lmekin}
\alias{fixef.lmekin}
\alias{ranef.lmekin}
\alias{vcov.lmekin}
\alias{VarCorr.lmekin}
\title{Extraction functions for Lmekin}
\description{Extract the fixed effects, random effects, variance of the
  fixed effects, or variance of the random effects from a linear mixed
  effects model fit with lmekin.}
\usage{
\method{fixef}{lmekin}(object, \dots)
\method{ranef}{lmekin}(object, \dots)
\method{vcov}{lmekin}(object, \dots)
\method{VarCorr}{lmekin}(x, \dots)
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
\value{the fixed effects are a vector and the variance of the fixed
  effects is a matrix.  The random effects will be a list with one
  element for each random effects terms, as will be their variance.}

\author{Terry Therneau}
\seealso{\code{\link{lmekin}}, \code{\link{random.effects}},
  \code{\link{fixed.effects}}
}
\examples{
efit <-  lmekin(effort ~ Type + (1|Subject), ergoStool)
ranef(efit)
}
\keyword{models}

