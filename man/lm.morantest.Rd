% Copyright 2001 by Roger S. Bivand
\name{lm.morantest}
\alias{lm.morantest}
\alias{listw2U}
\title{Moran's I test for residual spatial autocorrelation}
\description{
  Moran's I test for spatial autocorrelation in residuals from an estimated linear model (\code{lm()}). The helper function \code{listw2U()} constructs a weights list object corresponding to the sparse matrix \eqn{\frac{1}{2} ( \mathbf{W} + \mathbf{W}'}{1/2 (W + W')}
}
\usage{
lm.morantest(model, listw, zero.policy=FALSE)
listw2U(listw)
}
\arguments{
  \item{model}{an object of class \code{lm} returned by \code{lm}; it must currently include the x member from \code{lm(..., x=TRUE)}, weights and offsets should not be used}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
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
\references{Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 203, }
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{lm.LMtests}}, \code{\link{lm}}}

\examples{
data(oldcol)
oldcrime1.lm <- lm(CRIME ~ 1, data = COL.OLD, x = T, y = T)
oldcrime.lm <- lm(CRIME ~ HOVAL + INC, data = COL.OLD, x = T, y = T)
lm.morantest(oldcrime.lm, nb2listw(COL.nb, style="W"))
lm.LMtests(oldcrime.lm, nb2listw(COL.nb, style="W"))
lm.morantest(oldcrime.lm, nb2listw(COL.nb, style="S"))
lm.morantest(oldcrime1.lm, nb2listw(COL.nb, style="W"))
moran.test(COL.OLD$CRIME, nb2listw(COL.nb, style="W"), randomisation=FALSE)
}
\keyword{spatial}
