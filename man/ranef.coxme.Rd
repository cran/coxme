\name{ranef.coxme}
\alias{ranef.coxme}
\title{Extract coxme Random Effects}
\usage{
\method{ranef}{coxme}(object, \dots)
}
\arguments{
 \item{object}{an object inheriting from class \code{coxme}, representing
   a fitted Cox mixed-effects survival model.}
 \item{\dots}{some methods for this generic require additional
    arguments.  None are used in this method.} 
}
\description{
  The estimated random effects are a list with one element for
  each random term in the model.
  The list elements themselves may be a vector, matrix, or list.
}
\value{
  a list as described above.
}
\author{Terry Therneau}
\seealso{
  \code{\link{coxme}}, \code{\link{fixef.coxme}},
  \code{\link{random.effects}}
}

\examples{
rat1 <- coxme(Surv(time, status) ~ rx + (1|litter), rats)
ranef(rat1)
}
\keyword{models}
\keyword{survival}
