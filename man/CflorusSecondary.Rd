\name{CflorusSecondary}
\alias{CflorusSecondary}
\docType{data}
\title{
Secondary data on Colpoclypeus florus
}
\description{
Secondary data on 53 clutches of Colpoclypeus florus, a gregarious parasitoid wasp. The data were collected by Dijkstra 1986.

}
\usage{data(CflorusSecondary)}

\format{
The components are:

\tabular{lll}{
   [,1]\tab n   \tab   counts of the number of offspring in each clutch\cr
    [,2] \tab m \tab  counts of the number of male offspring in each clutch\cr
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
data(CflorusSecondary)
Meelis.test(CflorusSecondary)
}
\keyword{datasets}
