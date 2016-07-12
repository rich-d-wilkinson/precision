\name{rdoublebinom}
\alias{rdoublebinom}
\alias{ddoublebinom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{The Double Binomial Distribution
%%  ~~function to do ... ~~
}
\description{
Density and random generation for the double binomial distribution
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ddoublebinom(x = NULL, size, p, psi, log = FALSE)
rdoublebinom(n = 1, size, prob, psi)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
  number of observations
%%     ~~Describe \code{n} here~~
}
  \item{size}{
  number of trials
%%     ~~Describe \code{size} here~~
}

  \item{p}{
  probability parameter, same as prob
%%     ~~Describe \code{prob} here~~
}

  \item{prob}{
  probability parameter
%%     ~~Describe \code{prob} here~~
}
  \item{psi}{
  dispersion parameter. Negative values indicate over-dispersion and positive values give under-dispersed data.
%%     ~~Describe \code{psi} here~~
}
 \item{x}{
  value to evaluate density at. If null, the entire distribution is returned.
 }
 \item{log}{
 logical; if TRUE, probabilities p are given as log(p).
 }
}

\details{
Both functions first evaluate the density at all possible values of x, and then calculate the normalising constant. To generate random variables a uniform random variable is then called and the interval it falls into on the cumulative density is the value of the random variable returned.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\code{ddoublebinom} gives the density and \code{rdoublebinom} generates random deviates. 

If x is null, then the density at all values in the support is returned. Note that the indexing of this vector is out by one (i.e. the first element corresponds to x=0). 
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
See Efron 1986 for the density function.
%% ~put references to the literature/web site here ~
}
\author{Richard Wilkinson
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
rdoublebinom(n=10, size=12, p=0.1, psi=0.2)
ddoublebinom(x=1, size=12, p=0.1, psi=0.2)
ddoublebinom(size=12, p=0.1, psi=0.2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{double}
\keyword{binomial }% __ONLY ONE__ keyword per line
