\name{geary.test}
\alias{geary.test}
\title{Geary's C test for spatial autocorrelation}
\description{
   Geary's test for spatial autocorrelation using a spatial weights matrix in weights list form. The assumptions underlying the test are sensitive to the form of the graph of neighbour relationships and other factors, and results may be checked against those of \code{geary.mc} permutations.
}
\usage{
geary.test(x, listw, randomisation=TRUE, zero.policy=FALSE)
}

\arguments{
  \item{x}{a numeric vector the same length as the neighbours list in listw}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{randomisation}{variance of I calculated under the assumption of randomisation, if FALSE normality}
  \item{zero.policy}{if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
}

\value{
A coefficient matrix with values:
  \item{C}{Geary's C statistic}
  \item{EC}{Expectation}
  \item{VC}{Variance}
  \item{ZC}{Std. deviate}
  \item{PrC}{Pr(Z)}
}
\references{Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 21.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}
\note{Probability values one-tailed, as the hypothesis most often used is no autocorrelation against positive autocorrelation.}

\seealso{\code{\link{geary}}, \code{\link{geary.mc}}}

\examples{
require(spweights)
data(oldcol)
geary.test(COL.OLD$CRIME, nb2listw(COL.nb, style="W"))
geary.test(COL.OLD$CRIME, nb2listw(COL.nb, style="W"), randomisation=FALSE)
colold.lags <- nblag(COL.nb, 3)
geary.test(COL.OLD$CRIME, nb2listw(colold.lags[[2]], style="W"))
geary.test(COL.OLD$CRIME, nb2listw(colold.lags[[3]], style="W"))

}
\keyword{spatial}
