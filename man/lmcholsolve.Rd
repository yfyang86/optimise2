\name{lmcholsolve}
\alias{lmcholsolve}
\title{lmcholsolve}
\usage{
lmcholsolve(x, y)
}
\arguments{
  \item{x}{Must be a numeric matrix: integer matrix are NOT
  allowed.}

  \item{y}{Response variable}
}
\value{
\item{beta}{The estimations.} \item{Convergences}{Number of
non-zero diag-elements.}
}
\description{
This function solves \deqn{H \beta = y}, where H>=0;
}
\details{
Considering: \deqn{H \beta = y}, in each iteration newton's
scalor is \deqn{\beta_{LSE}}. The fatest algorithm is
'sweeping method'/Cholesky Decomposition. The mian reason
is that: H >= 0. If we tend not to check the condition,
then \deqn{H= LL^T}, where L is low tri-angle matrix. Hence
we derive \deqn{ [H,y]^T[H,y]=\left( \begin{array}{ll} L &
0\\ l^T & d \end{array} \right) \times \left(
\begin{array}{ll} L^T & l\\ 0 & d \end{array} \right)} and
\deqn{ L l= H^Ty , L^T \beta = l,
d^2=\|y-hat(y)\|_{\ell_2}^2} The performance cound be found
in inst/doc.
}
\note{
tol=1e-9 in Cholesky decomposition.
}
\examples{
x <- matrix(runif(100), 10)
y = rnorm(10) * 2
x <- t(x) \%*\% x
lmcholsolve(x, y)
solve(x, y)
}
\author{
Yifan Yang
}
\references{
Lange, Kenneth. Numerical analysis for statisticians. chap
7. Springer, 2010.
}
\keyword{optim2}
