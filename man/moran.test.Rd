\name{moran.test}
\alias{moran.test}
\title{Moran's I test for spatial autocorrelation}
\description{
 Moran's test for spatial autocorrelation using a spatial weights matrix in weights list form. The assumptions underlying the test are sensitive to the form of the graph of neighbour relationships and other factors, and results may be checked against those of \code{moran.mc} permutations.
}
\usage{
moran.test(x, listw, randomisation=TRUE, zero.policy=FALSE)
}
\arguments{
  \item{x}{a numeric vector the same length as the neighbours list in listw}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{randomisation}{variance of I calculated under the assumption of randomisation, if FALSE normality}
  \item{zero.policy}{if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
}

\value{
A coefficient matrix with values:
  \item{I}{Moran's I statistic}
  \item{EI}{Expectation}
  \item{VI}{Variance}
  \item{ZI}{Std. deviate}
  \item{PrI}{Pr(Z)}
}
\note{Var(I) is taken from Goodchild's CATMOG 47,
see also Upton & Fingleton (1985) p. 171; it agrees with SpaceStat,
see Tutorial workbook Chapter 22; VI is as given by Cliff and Ord minus the
square of EI. Probability values one-tailed, as the hypothesis most often used is no autocorrelation against positive autocorrelation.}

\references{Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 21.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{moran}}, \code{\link{moran.mc}}}

\examples{
require(spweights)
data(oldcol)
moran.test(COL.OLD$CRIME, nb2listw(COL.nb, style="W"))
moran.test(COL.OLD$CRIME, nb2listw(COL.nb, style="B"))
moran.test(COL.OLD$CRIME, nb2listw(COL.nb, style="C"))
moran.test(COL.OLD$CRIME, nb2listw(COL.nb, style="S"))
moran.test(COL.OLD$CRIME, nb2listw(COL.nb, style="W"), randomisation=FALSE)
colold.lags <- nblag(COL.nb, 3)
moran.test(COL.OLD$CRIME, nb2listw(colold.lags[[2]], style="W"))
moran.test(COL.OLD$CRIME, nb2listw(colold.lags[[3]], style="W"))
}
\keyword{spatial}
