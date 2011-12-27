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
sparse = c(1, 0), x = FALSE, y = TRUE,
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
  maximum likelihood}
  \item{sparse}{
    rule for deciding sparsity of a random effect.  See the main
    documentation for discussion of the issue.
}
  \item{x}{
    if TRUE the X matrix (fixed effects) is included in the output object
}
  \item{y}{
    if TRUE the y variable is included in the output object
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
  This routine was originally written as a check for the \code{coxme}
  routine; it uses the same code to process input arguments and form the
  random effects, comparison of its output with \code{lme} helped
  validate those operations.  It is possible to specify some models in
  this framwork that can not be fit with lme, in particular models with
  familial genetic effects, i.e., a \emph{kinship} matrix, and hence the
  name of the routine.  Using user-specified variance functions an even
  wider range of models is possible.

  For simple models the specification of the random effects follows the
  same form as the \code{lmer} function.  For any model which can be fit
  by both \code{lmekin} and \code{lmer}, the latter routine would
  normally be prefered due to a much wider selection of post-fit tools
  for residuals, prediction, plotting, etc.  
}
\value{
  An object of class \code{lmekin}.
}

\author{Terry Therneau}
\seealso{\code{\link{lmekin.object}}, \code{\link{coxme}}}
\examples{
  fit1 <- lme(effort ~ Type, data=ergoStool, random= ~1|Subject,
            method="ML")
  fit2 <- lmekin(effort ~ Type + (1|Subject), data=ergoStool)
}
\keyword{ models }

