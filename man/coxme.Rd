\name{coxme}
\alias{coxme}
\title{
Fit a mixed effects Cox model
}
\description{
  Fit a Cox model containing mixed (random and fixed) effects.  Assume
  a Gaussian distribution for the random effects.
}
\usage{
coxme(formula, data, weights, subset, na.action, init, control,
ties = c("efron", "breslow"),
varlist, vfixed, vinit, x = FALSE, y = TRUE,
refine.n = 0, random, fixed, variance, ...)
}
\arguments{
  \item{formula}{
    a two-sided formula with a survival object as the left hand side of a
    \code{~} operator and the fixed and random effects on the right.
}
  \item{data}{
    an optional data frame containing the variables named in the \code{formula}.
}
  \item{subset, weights, na.action}{
    further model specifications arguments as in \code{lm}; see there
    for details.
}
  \item{init}{
    optional initial values for the fixed effects.
}
  \item{control}{
   optional list of control options.  See \code{coxme.control} for details. 
}
  \item{ties}{
    method for handling exact ties in the survival time.
}
  \item{varlist}{
    the variance family to be used for each random term.  If there are
    multiple terms it will be a list of variance functions.
    The default is \code{coxmeFull}.
    Alternatively it can be a list of matrices, in which case the
    \code{coxmeMlist} function is used.
}
  \item{vfixed}{
    optional named list or vector used to fix the value of one or more
    of the variance terms at a constant.
}
  \item{vinit}{
    optional named list or vector giving suggested starting values for
    the variance.
}
  \item{x}{
    if TRUE the X matrix (fixed effects) is included in the output object
}
  \item{y}{
    if TRUE the y variable (survival time) is included in the output object
}
  \item{refine.n}{
    number of samples to be used in a monte-carlo estimate of 
    possible error in the log-likelihood of the fitted model due to
    inadequacy of the Laplace approximation.
}
  \item{fixed, random, variance}{
    In the preliminary version of \code{coxme}
    the fixed and random effects were separate arguments.  These arguments
    are included for backwards compatability, but are depreciated.
    The variance argument is a depreciated alias for vfixed.
}
  \item{\dots}{
    any other arguments are passed forward to \code{coxme.control}.
}
}
\value{
An object of class \code{coxme}.
}
\references{
S Ripatti and J Palmgren, Estimation of multivariate frailty 
models using penalized partial likelihood,
Biometrics, 56:1016-1022, 2000.

T Therneau, P Grambsch and VS Pankratz,
Penalized survival models and frailty,
J Computational and Graphical Statistics, 12:156-175, 2003.
}
\author{
Terry Therneau
}
\seealso{
  \code{\link{coxmeFull}}, \code{\link{coxmeMlist}},
  \code{\link{coxme.object}}
}
\examples{
\dontrun{# Random treatment effects per institution
fit1 <- coxme(Surv(pgtime, pgstat) ~ stage + trt + (1+trt |institution),
               data=colon2)
fit2 <- coxme(Surv(pgtime, pgstat) ~ stage + trt + (trt |institution) +
              strata(institution), data=colon2)
}
# Shrinkage effects (equivalent to ridge regression)
temp <- with(lung, scale(cbind(age, wt.loss, meal.cal)))
rfit <- coxme(Surv(time, status) ~ ph.ecog + (temp | 1), data=lung)
}

\keyword{survival}
