\name{Meelis.test}
\alias{Meelis.test}
\alias{James.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Meelis and James tests and R statistic
%%  ~~function to do ... ~~
}
\description{
Calculates the R statistic for count data, which is a ratio of the empirical variance to the expected binomial variance. Can take a mix of data of different sizes (i.e., different n in a Bin(n,p) distribution). Also runs the Meelis and James tests, testing for under or over dispersion. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Meelis.test(experiment, TwoSided = FALSE)
James.test(experiment, TwoSided = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{experiment}{The dataset. Should be a matrix with 2 columns. The first column contains the counts of the clutch sizes, and the second column contains the counts of the number of males in each clutch.
%%     ~~Describe \code{experiment} here~~
}
  \item{TwoSided}{
  logical; if TRUE, a two-sided test is performed. If we only want to test for under-disperion, a one-sided test could be used. 
%%     ~~Describe \code{TwoSided} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{U}{The U test statistic as described in Wilkinson et al. 2013.}
\item{p.val}{The p-value corresponding to the U statistic and whether the test is one or two-sided.}
\item{exp.table}{The experimental data organised into table format.}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Krackow et al. 2002, in 
Hardy 2002, Sex Ratios: Concepts and Research Methods, CUP.
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
data(GlegneriSecondary)
Meelis.test(GlegneriSecondary, TwoSided=TRUE)
James.test(GlegneriSecondary, TwoSided=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
