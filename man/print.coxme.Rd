\name{print.coxme}
\Rdversion{1.1}
\alias{print.coxme}
\alias{summary.coxme}
\title{
Print method for a coxme fit.
}
\description{
Print out the result of a coxme fit.
}
\usage{
\method{print}{coxme}(x, rcoef=FALSE, digits = options()$digits, ...)
\method{summary}{coxme}(object, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
    an object of class \code{coxme}, from the fit of a mixed effects
    survival model.
  }
  \item{rcoef}{print the random (penalized) coefficients, as well as
    the fixed ones.}
  \item{digits}{
    number of significant digits to print}
  \item{object}{
    an object of class \code{coxme}, from the fit of a mixed effects
    survival model.}
  \item{\dots}{
    optional arguments}
}
\note{
The summary function is currently identical to the print function.
The default in R is for print() to be the short form printout and
summary() the long form.
}
\author{
Terry Therneau
}
\seealso{\code{\link{coxme}}}
\keyword{survival}
