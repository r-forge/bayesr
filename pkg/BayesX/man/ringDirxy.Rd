\name{.ringDirxy}
\alias{.ringDirxy}

\title{Compute direction of points in a polygon.}

\description{
  This is (private, not exported) R code replacing the call to
  maptools:::.ringDirxy, which links to C code.
}

\usage{
.ringDirxy(xy)
}

\arguments{
  \item{xy}{a two-column matrix with x and y coordinates of at least 3 points which form a polygon}
}

\value{
  Returns the integer -1 if the direction is counter-clockwise and 1 if
  the direction is clockwise.
}

\note{
  The function uses vectorization and should not be very much slower
  than the original.
}

\author{Daniel Sabanes Bove}

\examples{\dontrun{## Create test polygon ordered clockwise.
nPoints <- 200
sequence <- seq(1, 10, length = nPoints)
firstHalf <- sequence[1:(nPoints/2)]

test <- cbind(
  "x" = c(sequence, rev(sequence)),
  "y" = c(firstHalf, rev(firstHalf), -firstHalf, -rev(firstHalf))
)

plot(test, type="n")
polygon(test)

## Compute direction with this and the original function,
BayesX:::.ringDirxy(test)
}}

\keyword{internal}

