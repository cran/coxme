\name{fixef.coxme}
\alias{fixef.coxme}
\alias{fixed.effects.coxme}
\title{Extract coxme Fixed Effects}
\description{Extract the fixed effects portion from a coxme model.}
\usage{
\method{fixef}{coxme}(object, \dots)
}
\alias{fixed.effects.coxme}
\arguments{
 \item{object}{an object inheriting from class \code{coxme}, representing
   the result of a mixed effects Cox model.
 }
 \item{\dots}{some methods for this generic require additional
   arguments.  None are used in this method.} 
}
\value{the vector of coefficients for the fixed effects}

\author{Terry Therneau}
\seealso{\code{\link{coxme}}, \code{\link{random.effects.coxme}}}

\examples{
rat1 <- coxme(Surv(time, status) ~ rx + (1|litter), rats)
fixed.effects(rat1)
}
\keyword{models}
\keyword{survival}
