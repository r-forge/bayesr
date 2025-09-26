\name{fatalities}
\alias{fatalities}

\title{Weekly Number of Fatalities in Austria}

\description{
  This data set includes weekly fatalities in Austria from 2000 to 2020. The data
  is taken from the Eurostat data base.
}

\usage{data("fatalities")}

\format{
The \code{fatalities} data contains the following variables:
\describe{
  \item{num:}{Integer, the number of fatalities.}
  \item{year:}{Integer, the corresponding year fatalities are recorded.}
  \item{week:}{Integer, the corresponding week fatalities are recorded.}
}
}

\details{
Previous versions of the package employed a slightly smaller version of the data
set. Since version 1.2-6 of the package the data also includes the last seven weeks
from 2020 (2020-47 to 2020-53) plus small corrections for the earlier weeks in
2020 and the addition of two previously missing weeks (2004-42 and 2014-20).
}

\references{
  Eurostat Database (2020). \emph{Population and social conditions, demography and migration,
    mortality, weekly deaths, deaths by week and NUTS 3 region, Austria}
    \url{https://ec.europa.eu/eurostat/}
}

\examples{
data("fatalities")
plot(num ~ week, data = fatalities)
}

\keyword{datasets}

