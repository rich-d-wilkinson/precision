\name{precision-package}
\alias{precision-package}
\alias{precision}
\docType{package}
\title{
Statistical tools for detecting over- or under-dispersion
}
\description{
Contains an implementation of the Meelis and James hypothesis tests for detecting under or over-dispersion in count data. Also calculates the R statistic. Implements the Bayesian model selection methods developed in Wilkinson et al 2013.
}
\details{
\tabular{ll}{
Package: \tab precision\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-10-6\cr
License: \tab GPL version 2 or newer\cr
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{

Richard Wilkinson  <r.d.wilkinson@nottingham.ac.uk>

}
\references{
Wilkinson et al. 2013.
}

\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
#demo("ModelEvidence", package="precision")
#demo("MCMC", package="precision")

library(precision)
data(GlegneriSecondary)
meelis.out <- Meelis.test(GlegneriSecondary, TwoSided=TRUE)
james.out = James.test(GlegneriSecondary, TwoSided=TRUE)


########### Run the MCMC
## prior for p is beta(a, b)
hyper$a.p <-  1 ## 20   ## Hyper parameters for p's prior distribution
hyper$b.p <-  1##100

# prior for psi - assumed to be Gaussian
hyper$mu.psi <- 0  
hyper$sd.psi <- 1

## Define a start point for the MCMC chain -
b.theta0 <-c(10, 0.1, 0.5) 
names(b.theta0) <-c("lambda", "p", "mort")

b.mcmc.out <- MCMCWithinGibbs( theta0=b.theta0,  data=GlegneriSecondary, hyper=hyper, nbatch=10^3, family="binomial", keepNM=TRUE)

# calculate the model evidence
b.log.evidence <- CalculateEvidence(mcmc.out=b.mcmc.out, data=GlegneriSecondary,  hyper=hyper, family="binomial", sd=FALSE)

}
