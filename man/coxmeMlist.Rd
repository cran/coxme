\name{coxmeMlist}
\Rdversion{1.1}
\alias{coxmeMlist}
\title{Coxme variance function}
\description{
  This variance function accepts a list of matrices, which define a
  correlation structure for a coxme fit.
}
\usage{
coxmeMlist(varlist, rescale = FALSE, pdcheck = TRUE, positive = TRUE)
}
\arguments{
  \item{varlist}{a list containing one or more matrix or bdsmatrix objects.}
  \item{rescale}{
    if TRUE, each input matrix is rescaled to have a diagonal of 1.
    (Kinship matrices for instance are often generated with a diagonal
    of .5 and would be multiplied by 2).
}
  \item{pdcheck}{
    check each matrix to ensure that it is positive definite
}
  \item{positive}{
    constrain coefficients to be positive.  This may also be a vector of
    the same length as \code{varlist}
}
}
\details{
  If two matrices \eqn{A}{A} and \eqn{B}{B} were given, this fits
  the variance structure
  \eqn{V= \sigma_1^2 A + \sigma_2^2 B}{V = s1 A + s2 B},
  where the variances \eqn{\sigma_1^2}{s1} and \eqn{\sigma_2^2}{s2} are
  parameters that will be optimized by \code{coxme},
  treating \eqn{A}{A} and \eqn{B}{B} as fixed.
}
\value{
a coxme variance family object, used by \code{coxme} in the fitting process.
}
\author{Terry Therneau}
\seealso{\code{\link{coxme}}}
\keyword{ survival }
