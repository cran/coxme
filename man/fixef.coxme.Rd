\name{fixef.coxme}
\alias{fixef.coxme}
\alias{ranef.coxme}
\alias{vcov.coxme}
\alias{VarCorr.coxme}
\title{Extraction functions for Coxme}
\description{Extract the fixed effects, randome effects, variance of the
  fixed effects, or variance of the random effects from a coxme model.}
\usage{
\method{fixef}{coxme}(object, \dots)
\method{ranef}{coxme}(object, \dots)
\method{vcov}{coxme}(object, \dots)
\method{VarCorr}{coxme}(x, \dots)
}
\arguments{
  \item{object}{an object inheriting from class \code{coxme} 
     representing the result of a mixed effects Cox model.
   }
   \item{x}{an object inheriting from class \code{coxme} 
     representing the result of a mixed effects Cox model.
   }
 \item{\dots}{some methods for this generic require additional
   arguments.  None are used in this method.} 
}
\value{the fixed effects are a vector and the variance of the fixed
  effects is a matrix.  The random effects will be a list with one
  element for each random effects terms, as will be their variance.}

\author{Terry Therneau}
\seealso{\code{\link{coxme}}, \code{\link{random.effects}},
  \code{\link{fixed.effects}}
}
\examples{
rat1 <- coxme(Surv(time, status) ~ rx + (1|litter), rats)
fixed.effects(rat1)
vcov(rat1)
random.effects(rat1)[[1]]  #one value for each of the 50 litters
VarCorr(rat1)
}
\keyword{models}
\keyword{survival}
