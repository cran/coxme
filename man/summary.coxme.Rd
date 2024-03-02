\name{summary.coxme}
\alias{summary.coxme}
\alias{print.summary.coxme}
\title{Summary method for coxme objects
}
\description{
Create a short summary of a coxme object
}
\usage{
\method{summary}{coxme}(object, ...)
\method{print}{summary.coxme}(x, digits = max(1L, getOption("digits") - 3L), 
    signif.stars = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{the result of a coxme fit}
  \item{x}{the result of a coxme fit}
  \item{summary.coxme}{a result of summary.coxme}
  \item{digits}{default number of digits to use}
  \item{signif.stars}{display the (idiotic) asterisks next to
    'significant' p values}
  \item{\dots}{other options to the internal \code{printCoefmat} call}
}
\details{
The summary function is a useful way to get the matrix of coefficients,
for other use.
}
\value{a \code{summary.coxme} object with components coefficients (a
  matrix with value, standard error, z, and p), random (a data frame
  containing sd of the random effects), and selected elements from the
  fit object
}
\seealso{\code{\link{fixef.coxme}}, \code{\link{ranef.coxme}},
    \code{\link{VarCorr}}
}
\keyword{survival}
