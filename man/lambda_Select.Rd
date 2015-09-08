\name{lambda_Select}
\alias{lambda_Select}
\title{Select Lambda using Dr Zhou's approach}
\usage{
lambda_Select(Yvec, Xmat, lambdas = seq(from = 0.001, to = 1, length.out =
  50), nPsep = 20, perc = 0.2,iidSampletype = "rnorm", cuda = FALSE, ...)
}
\arguments{
\item{Yvec}{Y vector(n)}

\item{Xmat}{Design Matrix(n by p)}

\item{lambdas}{Proposed lambda, default value is seq(from=.001,to=1,length.out=50)}

\item{nPsep}{Number of independent simulation data.}

\item{perc}{Stopping percentage}

\item{iidSampletype}{Random number generator.
It could be any R random number generator which 1 accpet only sample size as parameter,2 should be vectorized.}

\item{cuda}{Use cudaglmnet}
}
\value{
result Result from glmnet

selection T/F of lambdas. F means
}
\description{
Select Lambda using Dr Zhou's approach
}
\examples{
set.seed(65535)
Xmat = matrix(rnorm(100*80),ncol=80)
beta0 = rnorm(80,sd=2)
beta0[sample(1:80,70)] = 0.
epsilon = rnorm(100)
Yvec = Xmat\%*\%beta0 + epsilon
lambdas = seq(from=.001,to=1,length.out=50) # a vec
lambda_Select(Yvec,Xmat)
}

