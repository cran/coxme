\name{logLik.coxme}
\alias{logLik.coxme}
\title{The logLik method for coxme objects}

\description{logLik is most commonly used for a model fitted by maximum
  likelihood, and some uses, e.g. by AIC.
  This method allows generic functions to easily extract the
  log-likelhood of a coxme model.}

\usage{
\method{logLik}{coxme}(object, type = c("penalized", "integrated"), ...)
}
\arguments{
  \item{object}{a fitted coxme model}

  \item{type}{which of the two types of partial likelihood to extract}
  \item{\dots}{used by other methods}
}
\details{
  The likelihood for a mixed effects Cox model can be viewed in two ways:
  the ordinarly partial likelihood, where the random effects act only as
  a penalty or constraint, or a partial likelihood where the random
  effect has been integrated out.  Both are valid.
}
\value{Returns an object of class logLik.}
\seealso{\code{\link{logLik}}}
\keyword{ survival }
