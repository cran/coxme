\name{expand.nested}
\Rdversion{1.1}
\alias{expand.nested}
\title{Expand nested factors}
\description{
  Expand out the data frame for a nested factor such as (1| a/b).
  This is used by the variance function routines of \code{coxme}.
}
\usage{
expand.nested(x)
}
\arguments{
  \item{x}{A data frame containing the nesting variables}
}
\details{
  The initialize function of a \code{coxme} variance family is
  passed, as one of its arguments, a data frame \code{G} containing
  the grouping variables, each of which is a factor..
  Assume a nested factor \code{(1| a/b)} in the model formula and
  a data set whose first few rows are:
  \preformatted{
    a    b
    1    1
    1    2
    2    1
  }
  The function will replace the second column with a variable named
  \code{a/b} and values of 1/1, 1/2, 2/1, etc.
}
\value{an updated data frame}
\author{Terry Therneau}
\seealso{\code{\link{coxme}}, \code{\link{coxmeMlist}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ survival }
