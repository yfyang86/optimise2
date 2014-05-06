\docType{package}
\name{optimise2-package}
\alias{optimise2}
\alias{optimise2-package}
\title{Optimization in one dimensional problem.}
\description{
Minimizer: gsl min fminimizer quad golden This is a variant
of Brent’s algorithm which uses the safeguarded step-length
algorithm of Gill and Murray.
}
\details{
\tabular{ll}{ Package: \tab optimise2\cr Type: \tab
Package\cr Version: \tab 1.0\cr Date: \tab 2014-04-17\cr
License: \tab GPL2\cr } The usage is exactly the same as
optimise(stats) unless you want to enlarge the number of
iterations.
}
\examples{
f <- function(x) sin(x^2) + x/10
optimise(f, c(1, 4), tol = 1e-06)
optimise(f, c(1.5, 11))  # WRONG!!!
optimise2(f, c(1.5, 11), tol = 1e-06)
}
\author{
Yifan Yang

Maintainer: Yifan Yang <yifan.yang@uky.edu>
}
\references{
Richard Brent, Algorithms for minimization without
derivatives, Prentice-Hall (1973), republished by Dover in
paperback (2002), ISBN 0-486-41998-3.
}
\seealso{
optimise
}
\keyword{optimise2}
