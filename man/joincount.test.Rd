\name{joincount.test}
\alias{joincount.test}
\title{BB join count statistic for k-coloured factors}
\description{
 The BB join count test for spatial autocorrelation using a spatial weights matrix in weights list form for testing whether same-colour joins occur more frequently than would be expected if the zones were labelled in a spatially random way. The assumptions underlying the test are sensitive to the form of the graph of neighbour relationships and other factors, and results may be checked against those of \code{geary.mc} permutations.

}
\usage{
joincount.test(fx, listw)
}
\arguments{
  \item{fx}{a factor of the same length as the neighbours and weights objects in listw}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
}

}
\value{
A coefficient matrix with rows for each separate colour and columns:
  \item{1}{BB join count statistic}
  \item{2}{Expectation - nonfree sampling}
  \item{3}{Variance - nonfree sampling}
  \item{4}{Std. deviate (Z) - nonfree sampling}
  \item{5}{Pr(Z) - nonfree sampling}
}
\references{Cliff, A. D., Ord, J. K. 1981 Spatial processes, Pion, p. 20.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}
\note{Probability values one-tailed, as the hypothesis most often used is no autocorrelation against positive autocorrelation.}

\seealso{\code{\link{joincount}}, \code{\link{joincount.mc}}}

\examples{
require(spweights)
data(oldcol)
HICRIME <- cut(COL.OLD$CRIME, breaks=c(0,35,80), labels=c("low","high"))
joincount.test(HICRIME, nb2listw(COL.nb, style="B"))
joincount.test(HICRIME, nb2listw(COL.nb, style="C"))
joincount.test(HICRIME, nb2listw(COL.nb, style="S"))
joincount.test(HICRIME, nb2listw(COL.nb, style="W"))
by(card(COL.nb), HICRIME, summary)
}
\keyword{spatial}
