\name{coxmeFull}
\Rdversion{1.1}
\alias{coxmeFull}
\title{Variance family function for coxme fits.
}
\description{
  This function sets up the default variance family information for a
  mixed effects survival model fit with \code{coxme}.
}
\usage{
coxmeFull(collapse = FALSE)
}
\arguments{
  \item{collapse}{
    Form for fitting a nested effect, either standard or collapsed.  The
    latter appears to be more numerically stable (still under research).
}
}
\details{
Coxme variance families create a list with three functions: initialize,
generate, and wrapup, that determine how the variance structure of a fit
is modeled.
}
\value{
an object of class coxvar.
}
\author{
Terry Therneau
}

\seealso{
\code{\link{coxme}}
}
\keyword{survival }

