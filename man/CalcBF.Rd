\name{CalcBF}
\alias{CalcBF}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calculates Bayes factors and posterior probabilities from the model evidences.
%%  ~~function to do ... ~~
}
\description{
Calculates Bayes factors and posterior probabilities from the log model evidences.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
CalcBF(log.evidence)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{log.evidence}{A vector of length three containing the log.evidences for each of the three models. Should be ordered as \code{c(binomial, multbinom, doublebinom)}
%%     ~~Describe \code{log.evidence} here~~
}
}
\details{
Given a vector of the logarithm of the model evidences, \code{CalcBF} calculates the Bayes factors for all three combinations, and the posterior probabilities assuming that each model is equally likely a priori.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{BF}{Bayes factors}
\item{probH0}{Posterior probability if the three models are equally likely a priori.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
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
log.evidence<-c(-44, -38, -41)
CalcBF(log.evidence)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Bayes factor}
\keyword{posterior probability}% __ONLY ONE__ keyword per line
