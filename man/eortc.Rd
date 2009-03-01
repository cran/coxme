\name{eortc}
\Rdversion{1.1}
\alias{eortc}
\docType{data}
\title{Simulated data set based on an EORTC trial}
\description{
  This is a simulated surival data set for investigating random
  center effects.  To make it realistic, the number of centers and
  their sizes is based on an EORTC cancer trial. 
}
\usage{data(eortc)}
\format{
  A data frame with 2323 observations on the following 4 variables.
  \describe{
    \item{\code{y}}{survival time}
    \item{\code{uncens}}{0=alive, 1=dead}
    \item{\code{center}}{enrolling center, a number from 1 to 37}
    \item{\code{trt}}{treatment arm, 0 or 1}
  }
}
\details{ This is used in the test suite for the code. }
\source{PhD thesis work of Jose Cortinas Abrahantes}
\references{
Cortinas Abrahantes, Jose; Burzykowski, Tomasz (2002), A version of the
EM algorithm for proportional hazards models with random effects ,
Published in: Lecture Notes of the ICB Seminars. p. 15-20
}
\examples{
data(eortc)
coxme(Surv(y, uncens) ~ trt + (trt| center) + strata(center), eortc)
}
\keyword{datasets}
