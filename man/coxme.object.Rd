\name{coxme.object}
\alias{coxme.object}
\title{Coxme regression output object
}
\description{
  This contains further description of the output object created by a
  \code{coxme} call.  Most components can be accessed with extractor
  functions, which is the safer route since details of the object will
  likely change over time.
}
\value{
  \item{coefficients}{the coefficients of the fixed effects.  Use the
    \code{fixef} function to extract them.}
  \item{frail}{the coefficients of the random effects.  Use the
    \code{ranef} function to extract them.  These are always stored as a
    list with one member per random effect; each parenthesised term in
    the model is a random effect.
    In a linear mixed effects model the fixed effects and the variances
    of the random effects can be obtained without explicitly computing
    the coefficients of the random effects, the latter are called BLUP
    estimates and are computed later if at all.  This is not the case for a
    Cox model, there the random effect coefficients are a required part
    of the iteration process and so are always present in the final model.}
  \item{vcoef}{the variances of the random effects.  Use the
    \code{VarCor} function to extract them.  These are always stored as a
    list with one member per random effect.}
  \item{variance}{the variance-covariance matrix of the coefficient
    vector, including both fixed and random terms.  The random effects
    are listed first.  This will often be a sparse matrix.
    The \code{vcov} function
    will extract the fixed effects portion, which is always dense.}
  \item{loglik}{the log-likelihood vector from the fit.  The first
    element is the loglik at the initial values, the second is the
    integrated partial likelihood at the solution (IPL), the third is
    the penalized partial likelihood at the solution(PPL).}
  \item{df}{degrees of freedom for the IPL and the PPL solutions.}
  \item{hmat}{sparse Cholesky factorization of the information matrix.}
  \item{iter}{outer and inner iterations performed.  For each trial
    value of the variance parameters an Cox model partial likelihood
    must be solved; the outer iterations is the reported number from the
    optim() routine which handles the variance parameters, the inner
    iterations is the cumulative number of partial likelihood
    iterations.}
  \item{control}{a copy of the \code{coxph.control} parameters used in
    the fit.}
  \item{ties}{the computational method used for ties.}
  \item{u}{the vector of first derivatives of the PPL, at the solution.}
  \item{means,scale}{means and scale for each predictor, used internally
    to scale the problem.}
  \item{linear.predictor}{the vector of linear predictors.}
  \item{n}{vector containing the number of events and the number of
    observations in the fitting data set.}
  \item{terms}{the terms object from the fixed effects of the model
    formula.  Access using the \code{terms} function.}
  \item{formulaList}{the fixed and random portions of the formula,
    separated}
  \item{na.action}{the missing value attributes of the data, if any}
  \item{x,y,model}{optional: the x matrix, response, for model frame.
    These depend on the corresponding arguments in the call.}
  \item{call}{a copy of the call to the routine}
}
\details{The structure of each element of the random effects
  coefficients (obtained with \code{ranef}) and variances
  (\code{VarCorr}) depend on the variance functions, i.e., the functions
  used in the \code{varlist} argument.
  Since users can write their own variance functions this format can
  never be completely known.
  }
\keyword{survival}
