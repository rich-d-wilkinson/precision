\name{CflorusPrimary}
\alias{CflorusPrimary}
\docType{data}
\title{
Primary data on Colpoclypeus florus
}
\description{
Primary data on 55 clutches of Colpoclypeus florus, a gregarious parasitoid wasp. The data were collected by Dijkstra 1986.
}
\usage{data(CflorusPrimary)}

\format{
The components are:

\tabular{lll}{
   [,1]\tab N   \tab   counts of the number of eggs in each clutch\cr
    [,2] \tab M \tab  counts of the number of male eggs in each clutch\cr
      }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(CflorusPrimary)
Meelis.test(CflorusPrimary)
}
\keyword{datasets}
