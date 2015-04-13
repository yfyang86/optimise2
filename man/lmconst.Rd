\name{lmconst} 
\alias{lmconst}
\title{LSE with equal constraint}
\usage{
lmconst( y,x,x0,y0) 
}

\arguments{
\item{y}{ length n vector: response Variablei}
\item{x}{ n by p matrix: corresponding Variables}
\item{x0}{ p by l matrix: constraint beta coef}
\item{y0 }{ length l vector: value to constraints} 
}
\value{
\item{coef}{Constraint LSE/MLE(Normal assumption ). }
}
\description{
Y = X beta with x0 beta = y0 constraint. 
}
\details{
Use quadratic programming.
}
\note{
The function will only output the coefficients.
}
\examples{
# Example 1: dim =1
x=seq( 0,10,0.1) 
y = matrix( 6*x+7+runif( 101,min=-1,max=1) ,ncol=1) 
x=matrix( x,ncol=1) 
re = lmconst( y,x,2,19) 
re
re[ 1] +2*re[ 2] 

# Example 2: multi dim
# R CMD check will skipp %...
#x = matrix(runif(100),ncol=2) 
#y = x %*% c(6,2) +1 + rnorm( 50) 
#re = lmconst(y=y,x=x,x0=c(2,1),y0=15) 
#re
#re[1] +2*re[2] +re[ 3] 
}
\author{
Yifan Yang
}
\references{
Richard Brent, Algorithms for minimization without
derivatives, Prentice-Hall (1973), republished by Dover in
paperback (2002), ISBN 0-486-41998-3.
}
\keyword{lmconst} 

