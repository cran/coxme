\name{lmekin}
\alias{lmekin}
\title{Fit a linear mixed effects model}
\description{
  The lmekin function fits a linear mixed effects model, with random
  effects specified in the same structure as in the \code{coxme}
  function.
}
\usage{
lmekin(formula, data, weights, subset, na.action, control,
varlist, vfixed, vinit, method = c("ML", "REML"),
x = FALSE, y = FALSE, model=FALSE,
random, fixed, variance, ...)
}
\arguments{
  \item{formula}{
    a two-sided formula with the response as the left hand side of a
    \code{~} operator and the fixed and random effects on the right.
}
  \item{data}{
    an optional data frame containing the variables named in the \code{formula}.
}
  \item{subset, weights, na.action}{
    further model specifications arguments as in \code{lm}; see there
    for details.
}
  \item{control}{
   optional list of control options.  See \code{coxme.control} for details. 
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
  \item{method}{fit using either maximum likelihood or restricted
  maximum likelihood
}
  \item{x}{
    if TRUE the X matrix (fixed effects) is included in the output object
}
  \item{y}{
    if TRUE the y variable is included in the output object
}
  \item{model}{
    if TRUE the model frame is included in the output object
}
  \item{fixed, random, variance}{
    In an earlier version of \code{lmekin}
    the fixed and random effects were separate arguments.  These arguments
    are included for backwards compatability, but are depreciated.
    The variance argument is a depreciated alias for vfixed.
}
  \item{\dots}{
    any other arguments are passed forward to \code{coxme.control}.
}
}

\details{
  Let \eqn{A= \sigma^2 B}{A = sB} be the variance matrix of the random
  effects where \eqn{\sigma^2}{s} is the residual variance for the
  model.  Internally the routine solves for the parameters of
  \eqn{B}{B}, computing \eqn{A}{A} at the end.  The \code{vinit} and
  \code{vfixed} parmaters refer to \eqn{B}{B}, however.
 
  It is possible to specify certain models in \code{lmekin}
  that can not be fit with lme, in particular models with
  familial genetic effects, i.e., a \emph{kinship} matrix, and hence the
  name of the routine.  Using user-specified variance functions an even
  wider range of models is possible.
  For simple models the specification of the random effects follows the
  same form as the \code{lmer} function.  For any model which can be fit
  by both \code{lmekin} and \code{lmer}, the latter routine would
  normally be prefered due to a much wider selection of post-fit tools
  for residuals, prediction and plotting.

  Much of the underlying model code for specification and manipulation
  of the random effects is shared with the \code{coxme} routine.  In
  fact lmekin was originally written only to provide a test routine for
  those codes, and no expectation that it would find wider utility.
}
\value{
  An object of class \code{lmekin}.
}

\author{Terry Therneau}
\seealso{\code{\link{lmekin.object}}, \code{\link{coxme}}}
\examples{
  data(ergoStool, package="nlme")  # use a data set from nlme 
  fit1 <- lmekin(effort ~ Type + (1|Subject), data=ergoStool)
  \dontrun{
  # gives the same result
  require(nlme)
  fit2 <- lme(effort ~ Type, data=ergoStool, random= ~1|Subject,
            method="ML")
  } 
}
\keyword{ models }

