\name{lmekin.control}
\alias{lmekin.control}
\title{
Auxillary parameters for controlling lmekin fits.
}
\description{
Auxillary function which packages the optional parameters of a
\code{lmekin} fit as a single list.
}
\usage{
lmekin.control(
optpar = list(method = "BFGS", control=list(reltol = 1e-8)),
varinit=c(.02, .1, .8, 1.5)^2, corinit = c(0, .3)) 
}
\arguments{
  \item{optpar}{
    parameters passed forward to the \code{optim} routine.
  }
  \item{varinit}{the default grid of starting values for variances, used if no
    \code{vinit} argument is supplied in the \code{lmekin} call.}
  \item{corinit}{the default grid of starting values for correlations.}
}
\details{
 The main flow of \code{lmekin} is to use the \code{optim} routine to
 find the best values for the variance parameters.  For any given trial
 value of the variance parameters, a subsidiary computation maximizes
 the likelihood to select the regression coefficients beta (fixed) and b
 (random).

 If no starting values are supplied for the variances of the random
 effects then a grid search is performed to select initial values for
 the main iteration loop.
 The variances and correlations are all scaled by
 \eqn{\sigma^2}{sigma^2},
 making these starting estimates scale free, e.g., replacing y by 10*y in a
 data set will change \eqn{\sigma}{sigma} but not the internal
 representation of any other variance parameters.
 Because we use the log(variance) as our iteration scale the 0--.001
 portion of the
 variance scale is stretched out giving a log-likelihood surface that is almost
 flat; a Newton-Raphson iteration starting at log(.2) may have log(.0001) as its
 next guess and get stuck there, never finding a true maximum that lies in the
 range of .01 to .05.
 Corrleation paramters seem to need fewer starting points.
}
\value{
a list of control parameters}
\author{
Terry Therneau
}
\seealso{
\code{\link{lmekin}}
}
\keyword{survival}
