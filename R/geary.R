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


geary <- function(x, listw, n, n1, S0, zero.policy=FALSE) {
	z <- scale(x, scale=F)
	zz <- sum(z^2)
	K <- (n*sum(z^4))/(zz^2)
	res <- as.numeric(rep(0, n))
	for (i in 1:n) {
		if (length(listw$neighbours[[i]]) == 0) {
			if (zero.policy) res[i] <- 0
			else res[i] <- NA
		} else {
			for (j in 1:length(listw$neighbours[[i]])) {
				jl <- listw$neighbours[[i]]
				res[i] <- res[i] + (listw$weights[[i]][j] * 
					(x[i] - x[jl[j]])^2)
			}
		}
	}
	C <- (n1 / (2*S0)) * (sum(res) / zz)
	res <- list(C=C, K=K)
	res
}

geary.test <- function(x, listw, randomisation=TRUE, zero.policy=FALSE) {
	if(class(listw) != "listw") stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if(!is.numeric(x)) stop(paste(deparse(substitute(x)),
		"is not a numeric vector"))
	if (any(is.na(x))) stop("NA in X")
	n <- length(listw$neighbours)
	if (n != length(x)) stop("objects of different length")
	wc <- spweights.constants(listw)
	S02 <- wc$S0*wc$S0
	res <- geary(x, listw, wc$n, wc$n1, wc$S0, zero.policy)
	C <- res$C
	K <- res$K
	EC <- 1
	if(randomisation) {
		VC <- (wc$n1*wc$S1*(wc$nn - 3*n + 3 - K*wc$n1))
		VC <- VC - ((1/4) * (wc$n1*wc$S2*(wc$nn + 3*n - 6 - 
			K*(wc$nn - n + 2))))
		VC <- VC + (S02*(wc$nn - 3 - K*(wc$n1^2)))
		VC <- VC / (n*wc$n2*wc$n3*S02)
	} else {
		VC <- ((2*wc$S1 + wc$S2)*wc$n1 - 4*S02) / (2*(n + 1)*S02)
	}
	ZC <- (C - EC) / sqrt(VC)
	PrC <- 2*(1-pnorm(abs(ZC)))
	res <- matrix(0, nrow=1, ncol=5)
	rownames(res) <- deparse(substitute(x))
	colnames(res) <- c("Geary C statistic", "Expectation", "Variance",
		"Std. deviate", "Pr(Z)")
	res[1] <- C
	res[2] <- EC
	res[3] <- VC
	res[4] <- ZC
	res[5] <- PrC
	thiscall <- match.call()
	attr(res, "call") <- thiscall
	neighbours.attrs <- names(attributes(listw$neighbours))
	attr(res, "neighbours.attrs") <- neighbours.attrs
	weights.attrs <- names(attributes(listw$weights))
	attr(res, "weights.attrs") <- weights.attrs
	cat("Geary's C test under", ifelse(randomisation, "randomisation",
		"normality"), "for", deparse(substitute(x)),
		"\nwith call:\n") 
	print(thiscall)
	cat("neighbours list attributes:", neighbours.attrs, "\n")
	cat("and weights list attributes:", weights.attrs, "\n\n")
	print.coefmat(res, signif.stars = FALSE)
	cat("\n")
	invisible(res)
}

geary.mc <- function(x, listw, nsim, zero.policy=FALSE) {
	if(class(listw) != "listw") stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if(!is.numeric(x)) stop(paste(deparse(substitute(x)),
		"is not a numeric vector"))
	if(missing(nsim)) stop("nsim must be given")
	if (any(is.na(x))) stop("NA in X")
	n <- length(listw$neighbours)
	if (n != length(x)) stop("objects of different length")
	wc <- spweights.constants(listw)
	res <- numeric(length=nsim+1)
	for (i in 1:nsim) res[i] <- geary(sample(x), listw, n, wc$n1, wc$S0, zero.policy)$C
	res[nsim+1] <- geary(x, listw, n, wc$n1, wc$S0, zero.policy)$C
	rankres <- rank(res)
	xrank <- rankres[length(res)]
	lres <- list(res=res, rankres=rankres, xrank=xrank)
	thiscall <- match.call()
	attr(lres, "call") <- thiscall
	neighbours.attrs <- names(attributes(listw$neighbours))
	attr(lres, "neighbours.attrs") <- neighbours.attrs
	weights.attrs <- names(attributes(listw$weights))
	attr(lres, "weights.attrs") <- weights.attrs
	cat("Geary's C simulation for", deparse(substitute(x)),
		"\nObserved statistic", format(res[nsim+1]), "ranks",
		 xrank, "among", length(res), "values\nwith call:\n")
	print(thiscall)
	cat("neighbours list attributes:", neighbours.attrs, "\n")
	cat("and weights list attributes:", weights.attrs, "\n\n")
	invisible(lres)
}

