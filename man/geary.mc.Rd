% Copyright 2001 by Roger S. Bivand
\name{geary.mc}
\alias{geary.mc}
\title{Permutation test for Geary's C statistic}
\description{
 A permutation test for Geary's C statistic calculated by using nsim random permutations of x for the given spatial weighting scheme, to establish the rank of the observed statistic in relation to the nsim simulated values.
}
\usage{
geary.mc(x, listw, nsim, zero.policy=FALSE)
}
\arguments{
  \item{x}{a numeric vector the same length as the neighbours list in listw}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{nsim}{number of permutations}
  \item{zero.policy}{if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
}

\value{
a list with
  \item{res}{nsim simulated values of statistic, final value is observed statistic}
  \item{rankres}{ranks of res}
  \item{xrank}{rank of observed statistic}
}
\references{Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 63-5.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{geary}}, \code{\link{geary.test}}}

\examples{
require(spweights)
data(oldcol)
sim1 <- geary.mc(COL.OLD$CRIME, nb2listw(COL.nb, style="W"), nsim=99)
mean(sim1$res)
var(sim1$res)
summary(sim1$res)
colold.lags <- nblag(COL.nb, 3)
sim2 <- geary.mc(COL.OLD$CRIME, nb2listw(colold.lags[[2]], style="W"),
 nsim=99)
summary(sim2$res)
sim3 <- geary.mc(COL.OLD$CRIME, nb2listw(colold.lags[[3]], style="W"),
 nsim=99)
summary(sim3$res)
}
\keyword{spatial}
