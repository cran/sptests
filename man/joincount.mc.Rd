\name{joincount.mc}
\alias{joincount.mc}

\title{Permutation test for same colour join count statistics}
\description{
A permutation test for same colour join count statistics calculated by using nsim random permutations of fx for the given spatial weighting scheme, to establish the ranks of the observed statistics (for each colour) in relation to the nsim simulated values.
}
\usage{
joincount.mc(fx, listw, nsim)
}
\arguments{
  \item{fx}{a factor of the same length as the neighbours and weights objects in listw}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{nsim}{number of permutations}
}

\value{
a list with
  \item{res}{nsim simulated values of statistic for each colour, final row are the observed statistics}
  \item{rankres}{ranks of res columns}
  \item{xrank}{ranks of observed statistics}
}
\references{Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 63-5.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{joincount}}, \code{\link{joincount.test}}}

\examples{
require(spweights)
data(oldcol)
HICRIME <- cut(COL.OLD$CRIME, breaks=c(0,35,80), labels=c("low","high"))
joincount.mc(HICRIME, nb2listw(COL.nb, style="B"), nsim=99)
joincount.test(HICRIME, nb2listw(COL.nb, style="B"))
}
\keyword{spatial}
