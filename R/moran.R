# Copyright 2001 by Roger Bivand 
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#

moran <- function(x, listw, n, S0, zero.policy=FALSE) {
	z <- scale(x, scale=F)
	zz <- sum(z^2)
	K <- (n*sum(z^4))/(zz^2)
	lz <- lag.listw(listw, z, zero.policy)
	I <- (n / S0) * ((t(z) %*% lz) / zz)
	res <- list(I=I, K=K)
	res
}

moran.test <- function(x, listw, randomisation=TRUE, zero.policy=FALSE) {
	if (class(listw) != "listw") stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if (!is.numeric(x)) stop(paste(deparse(substitute(x)),
		"is not a numeric vector"))
	if (any(is.na(x))) stop("NA in X")
	n <- length(listw$neighbours)
	if (n != length(x)) stop("objects of different length")
	wc <- spweights.constants(listw)
	S02 <- wc$S0*wc$S0
	res <- moran(x, listw, wc$n, wc$S0, zero.policy)
	I <- res$I
	K <- res$K
	EI <- (-1) / wc$n1
	if(randomisation) {
		VI <- n*(wc$S1*(wc$nn - 3*n + 3) - n*wc$S2 + 3*S02)
		tmp <- K*(wc$S1*(wc$nn - n) - 2*n*wc$S2 + 6*S02)
		VI <- (VI - tmp) / (wc$n1*wc$n2*wc$n3*S02)
		VI <- VI - EI^2
	} else {
		VI <- (wc$nn*wc$S1 - n*wc$S2 + 3*S02) / (S02*(wc$nn - 1))
		VI <- VI - EI^2
	}
	ZI <- (I - EI) / sqrt(VI)
	PrI <- 2*(1-pnorm(abs(ZI)))
	res <- matrix(0, nrow=1, ncol=5)
	rownames(res) <- deparse(substitute(x))
	colnames(res) <- c("Moran's I statistic", "Expectation", "Variance",
		"Std. deviate", "Pr(Z)")
	res[1] <- I
	res[2] <- EI
	res[3] <- VI
	res[4] <- ZI
	res[5] <- PrI
	thiscall <- match.call()
	attr(res, "call") <- thiscall
	neighbours.attrs <- names(attributes(listw$neighbours))
	attr(res, "neighbours.attrs") <- neighbours.attrs
	weights.attrs <- names(attributes(listw$weights))
	attr(res, "weights.attrs") <- weights.attrs
	cat("Moran's I test under", ifelse(randomisation, "randomisation",
		"normality"), "for", deparse(substitute(x)),
		"\nwith call:\n") 
	print(thiscall)
	cat("neighbours list attributes:", neighbours.attrs, "\n")
	cat("and weights list attributes:", weights.attrs, "\n\n")
	print.coefmat(res, signif.stars = FALSE)
	cat("\n")
	invisible(res)
}

moran.mc <- function(x, listw, nsim, zero.policy=FALSE) {
	if(class(listw) != "listw") stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if(!is.numeric(x)) stop(paste(deparse(substitute(x)),
		"is not a numeric vector"))
	if(missing(nsim)) stop("nsim must be given")
	if (any(is.na(x))) stop("NA in X")
	n <- length(listw$neighbours)
	if (n != length(x)) stop("objects of different length")
	S0 <- Szero(listw)
	res <- numeric(length=nsim+1)
	for (i in 1:nsim) res[i] <- moran(sample(x), listw, n, S0, zero.policy)$I
	res[nsim+1] <- moran(x, listw, n, S0, zero.policy)$I
	rankres <- rank(res)
	xrank <- rankres[length(res)]
	lres <- list(res=res, rankres=rankres, xrank=xrank)
	thiscall <- match.call()
	attr(lres, "call") <- thiscall
	neighbours.attrs <- names(attributes(listw$neighbours))
	attr(lres, "neighbours.attrs") <- neighbours.attrs
	weights.attrs <- names(attributes(listw$weights))
	attr(lres, "weights.attrs") <- weights.attrs
	cat("Moran's I permutation test for", deparse(substitute(x)),
		"\nObserved statistic", format(res[nsim+1]), "ranks",
		xrank, "among", length(res), "values\nwith call:\n")
	print(thiscall)
	cat("neighbours list attributes:", neighbours.attrs, "\n")
	cat("and weights list attributes:", weights.attrs, "\n\n")
	invisible(lres)
}


