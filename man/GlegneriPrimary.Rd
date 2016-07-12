\name{GlegneriPrimary}
\alias{GlegneriPrimary}
\docType{data}
\title{
Primary data on Goniouzus legneri
}
\description{
Primary data on 47 clutches of Goniouzus legneri, a gregarious parasitoid wasp. The data were collected by Khidr et al 2013.
}
\usage{data(GlegneriPrimary)}

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
data(GlegneriPrimary)
Meelis.test(GlegneriPrimary)
}
\keyword{datasets}
