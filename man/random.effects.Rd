\name{random.effects}
\Rdversion{1.1}
\alias{random.effects}
\alias{ranef}
\title{Extract Random Effects}
\description{
     This function is generic; method functions can be written to
     handle specific classes of objects. Classes which already have
     methods for this function include 'coxme', 'lmList' and 'lme'.
}
\usage{
random.effects(object, ...)
ranef(object, ...)
}
\arguments{
  \item{object}{any fitted model object from which random effects estimates
          can be extracted.}
  \item{\dots}{some methods for this generic function require additional
          arguments.}
}
\value{will depend on the method function used; see the appropriate
     documentation.
}
\references{
Pinheiro, J.C., and Bates, D.M. (2000) "Mixed-Effects Models in S
     and S-PLUS", Springer, esp. pp. 100, 461.
}
\author{Terry Therneau}
\seealso{\code{\link{ranef.coxme}}}
\examples{
 ## see the methods function documentation
}
\keyword{ models }
\keyword{ survival }
