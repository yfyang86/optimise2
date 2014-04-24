\name{optimise2}
\alias{optimise2}

\title{
OPTIMISE2
}
\description{
One-dimensional optimization problem.
}
\usage{
optimise2(f, interval, ..., lower = min(interval), upper = max(interval), maximum = FALSE, tol = .Machine$double.eps^0.25, initials = runif(1, lower, upper), trace = F, maxit = 100)
}

\arguments{
  \item{f}{
target function
}
  \item{interval}{
Interval used to identify the optimization problem, can't be inf!
}
  \item{\dots}{
Other parameters used by f.
}
  \item{lower}{
Lower bound of the interval.
}
  \item{upper}{
Upper bound of the interval.
}
  \item{maximum}{
Max number of itterations performed.
}
  \item{tol}{
|a - b| < tol + epsrel min(|a|,|b|) . epsrel=0 is my case.
}
  \item{initials}{
Testing, no use right now.
}
  \item{trace}{
Testing, no use right now.
}
  \item{maxit}{
Max number of simulations used.

}
}
\details{
Minimizer: gsl_min_fminimizer_goldensection
The golden section algorithm is the simplest method of bracketing the minimum of a function. It is the slowest algorithm provided by the library, with linear convergence.

On each iteration, the algorithm first compares the subintervals from the endpoints to the current minimum. The larger subinterval is divided in a golden section (using the famous ratio (3-\sqrt 5)/2 = 0.3189660…) and the value of the function at this new point is calculated. The new value is used with the constraint f(a') > f(x') < f(b') to a select new interval containing the minimum, by discarding the least useful point. This procedure can be continued indefinitely until the interval is sufficiently small. Choosing the golden section as the bisection ratio can be shown to provide the fastest convergence for this type of algorithm.

Minimizer: gsl_min_fminimizer_brent
The Brent minimization algorithm combines a parabolic interpolation with the golden section algorithm. This produces a fast algorithm which is still robust.

The outline of the algorithm can be summarized as follows: on each iteration Brent’s method approximates the function using an interpolating parabola through three existing points. The minimum of the parabola is taken as a guess for the minimum. If it lies within the bounds of the current interval then the interpolating point is accepted, and used to generate a smaller interval. If the interpolating point is not accepted then the algorithm falls back to an ordinary golden section step. The full details of Brent’s method include some additional checks to improve convergence.

Minimizer: gsl_min_fminimizer_quad_golden
This is a variant of Brent’s algorithm which uses the safeguarded step-length algorithm of Gill and Murray.

Notice that, in original optimise(stats), the condition "f(a') > f(x') < f(b')" will not be checked.
}
\value{
\item{maximum}{Optimization point.}
\item{minimum}{Optimization point.}
\item{objective}{Function value.}
}
\references{
Richard Brent, Algorithms for minimization without derivatives, Prentice-Hall (1973), republished by Dover in paperback (2002), ISBN 0-486-41998-3.
}
\author{
Yifan Yang
}
\note{
The function will ignore the discontinouse points.
}

\examples{
     f  <- function(x) sin(x^2)+x/10
     optimise(f, c(1, 4),tol=0.000001)  
     optimise(f, c(1.5,11 ))  #WRONG
     optimise2(f, c(1.5,11 ),tol=0.000001) 
}

\keyword{ optim2 }
