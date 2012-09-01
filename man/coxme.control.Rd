\name{coxme.control}
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
iter.max = 20, inner.iter = Quote(max(4, fit0$iter+1)),
sparse.calc = NULL,
optpar = list(method = "BFGS", control=list(reltol = 1e-5)),
refine.df=4, refine.detail=FALSE, refine.method="control",
sparse=c(50, .02), 
varinit=c(.02, .1, .4, .8)^2, corinit = c(0, .3)) 
}
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
    partial likelihood is the objective function of \code{optim}.
    The default is to use one more iteration than the baseline coxph
    model \code{fit0}.  The baseline model contains only the fixed
    effects, and is as part of the setup by the main program.
    The minimum value of 4 applies most often to the case where there
    are no fixed effects.
}
  \item{sparse.calc}{
    choice of method 1 or 2 for a particular portion of the
    calculation.  This can have an effect on run time for problems with
    thousands of random effects.
}
  \item{optpar}{
    parameters passed forward to the \code{optim} routine.
  }
  \item{refine.df}{the degrees of freedom for the t-distribution
    used to draw random samples
    for the refine.n option}
  \item{refine.detail}{this option is mostly for debugging.  If TRUE
    then an extra component \code{refine.detail} will be present in
    the output which contains intermediate variables from the
    iterative refinement calculation.}
  \item{refine.method}{method by which the control calculations are
    done. This is a current research/development question, the option
    will likely disappear at some future date, and users should ignore
    it.}
    \item{sparse}{rule for deciding sparsity of a random effect, see
      details below.}  
  \item{varinit}{the default set of starting values for variances, used if no
    \code{vinit} argument is supplied in the \code{coxme} call.}
  \item{corinit}{the default set of starting values for correlations.}
}
\details{
 The main flow of \code{coxme} is to use the \code{optim} routine to
 find the best values for the variance parameters.  For any given trial
 value of the variance parameters, an inner loop maximizes the partial
 likelihood to select the regression coefficients beta (fixed) and b
 (random). Within this loop cholesky decomposition is used.  It is
 critical that the convergence criteria of inner loops be less than
 outer ones, thus toler.chol < eps < reltol.

 If no starting values are supplied for the variances of the random
 effects then a grid search is performed to select initial values for
 the main iteration loop.
 The default values given here are based on experience but without any
 formal arguments for their optimality.
 We have found that the estimated standard
 deviation of a random effect is often between .1 and .3, corresponding to
 exp(.1)= 1.1 to exp(.3)= 1.35 fold ``average'' relative risks associated with
 group membership.  This is  biologically reasonable for a latent trait.
 Other common solutions ane a small random effect corresponding to only 1--5\%
 change in the hazard or likelihood that is maximized at the boundary
 value of 0 variance.
 Variances greater than 2 are very unusual.
 Because we use the log(variance) as our iteration scale the 0--.001 portion of the
 variance scale is stretched out giving a log-likelihood surface that is almost
 flat; a Newton-Raphson iteration starting at log(.2) may have log(.0001) as its
 next guess and get stuck there, never finding a true maximum that lies in the
 range of .01 to .05.
 Corrleation paramters seem to need fewer starting points.

 The sparse option controls a sparse approximation in the code.
 Assume we have a mixed effects model with a random intercept per group,
 and there are 1000 groups.
 In a Cox model (unlike a linear mixed effects model) the resulting
 second derivative matrix used during the solution will be 1000 by 1000
 with no zeros, and
 fitting the model can consume a large amount of both time and memory.
 However, it is almost sparse, in that elements off the diagonal are
 very small and can often be ignored.  Computation with a sparse
 matrix approximation will be many times faster.
 Luckily, as the number of groups increases the accuracy of the
 approximation also increases.
 If \code{sparse=c(50, .03)} this states that sparse approximation will
 be employed for any grouping variable with 50 or more levels, and off
 diagonal elements that relate any two levels both of which represent
 .03 of less of the total sample will be treated as zero.
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
