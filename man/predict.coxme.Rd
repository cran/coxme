\name{predict.coxme}
\alias{predict.coxme}
\title{Predictions for a coxme object.}	
\description{Return predicted values from a coxme fit.}
\usage{
\method{predict}{coxme}(object, newdata, type = c("lp", "risk"), \ldots)
}
\arguments{
  \item{object}{an object of class \code{coxme}, 
     from the fit of a mixed effects survival model}
  \item{newdata}{new data set, not yet supported}
  \item{type}{type of prediction}
  \item{\ldots}{arguments for other methods}
}
\value{
  a vector of predicted values
}
\seealso{\code{\link{coxme}}}
\keyword{survival}
