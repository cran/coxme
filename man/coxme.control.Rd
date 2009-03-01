\name{coxme.control}
\Rdversion{1.1}
\alias{coxme.control}
\title{
Auxillary parameters for controlling coxme fits.
}
\description{
Auxillary function which packages the optional parameters of a
\code{coxme} fit as a single list.
}
\usage{
coxme.control(eps = 1e-08, toler.chol = .Machine$double.eps^0.75,
iter.max = 20, inner.iter = 5, sparse.calc = NULL,
optpar = list(method = "BFGS", control=list(reltol = 1e-5)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eps}{
    convergence criteria for the partial likelihood
}
  \item{toler.chol}{
    tolerance for the underlying Cholesky decomposition.  This is used
    to detect singularity (redundant variables).
}
  \item{iter.max}{
    maximum number of iterations for the final fit
}
  \item{inner.iter}{
    number of iterations for the `inner loop' fits, i.e. when the
    partial likelihood is the objective function of \code{optim}
}
  \item{sparse.calc}{
    choice of method 1 or 2 for a particular portion of the
    calculation.  This can have an effect on run time for problems with
    thousands of random effects.
}
  \item{optpar}{
    parameters passed forward to the \code{optim} routine.
}
}
\details{
 The main flow of \code{coxme} is to use the \code{optim} routine to
 find the best values for the variance parameters.  For any given trial
 value of the variance parameters, an inner loop maximizes the partial
 likelihood to select the regression coefficients beta (fixed) and b
 (random). Within this loop cholesky decomposition is used.  It is
 critical that the convergence criteria of inner loops be less than
 outer ones, thus toler.chol < eps < reltol.
}
\value{
a list of control parameters}
\author{
Terry Therneau
}
\seealso{
\code{\link{coxme}}
}
\keyword{survival}
