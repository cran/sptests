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

joincount <- function(dums, listw) {
	nc <- ncol(dums)
	n <- length(listw$neighbours)
	cardnb <- card(listw$neighbours)
	res <- as.numeric(rep(0, nc))
	for (lev in 1:nc) {
		for (i in 1:n) {
			xi <- dums[i, lev]
			if (cardnb[i] > 0)
				res[lev] <- res[lev] + (dums[i, lev] *
				sum(dums[listw$neighbours[[i]], lev] *
				listw$weights[[i]]))
		}
	}
	res
}

joincount.test <- function(fx, listw) {
	if (class(listw) != "listw") stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if (!is.factor(fx)) stop(paste(deparse(substitute(x)),
		"is not a factor"))
	if (any(is.na(fx))) stop("NA in factor")
	n <- length(listw$neighbours)
	if (n != length(fx)) stop("objects of different length")
	wc <- spweights.constants(listw)
	S02 <- wc$S0*wc$S0

	dums <- lm(codes(fx) ~ fx - 1, x=TRUE)$x
	BB <- joincount(dums, listw)
	res <- matrix(0, nrow=length(BB), ncol=5)
	tab <- table(fx)
	rownames(res) <- names(tab)
	colnames(res) <- c("Same colour statistic", "Expectation", "Variance",
		"Std. deviate", "Pr(Z)")
	res[,1] <- 0.5 * BB
	ntab <- as.vector(tab)
	Ejc <- (wc$S0*(ntab*(ntab-1))) / (2*n*wc$n1)
	Vjc <- (wc$S1*(ntab*(ntab-1))) / (n*wc$n1)
	Vjc <- Vjc + (((wc$S2 - 2*wc$S1)*ntab*(ntab-1)*(ntab-2)) /
		(n*wc$n1*wc$n2))
	Vjc <- Vjc + (((S02 + wc$S1 - wc$S2)*ntab*(ntab-1)*(ntab-2)*
		(ntab-3)) / (n*wc$n1*wc$n2*wc$n3))
	Vjc <- (0.25 * Vjc) - Ejc^2
	res[,2] <- Ejc
	res[,3] <- Vjc
	res[,4] <- (res[,1] - res[,2]) / sqrt(res[,3])
	res[,5] <- 2*(1-pnorm(abs(res[,4])))
	thiscall <- match.call()
	attr(res, "call") <- thiscall
	neighbours.attrs <- names(attributes(listw$neighbours))
	attr(res, "neighbours.attrs") <- neighbours.attrs
	weights.attrs <- names(attributes(listw$weights))
	attr(res, "weights.attrs") <- weights.attrs
	cat("Join count test under nonfree sampling for",
		 deparse(substitute(fx)), "\nwith call:\n") 
	print(thiscall)
	cat("neighbours list attributes:", neighbours.attrs, "\n")
	cat("and weights list attributes:", weights.attrs, "\n\n")
	print.coefmat(res, signif.stars = FALSE)
	cat("\n")
	invisible(res)
}

joincount.mc <- function(fx, listw, nsim) {
	if(class(listw) != "listw") stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if(!is.factor(fx)) stop(paste(deparse(substitute(x)),
		"is not a factor"))
	if(missing(nsim)) stop("nsim must be given")
	if (any(is.na(fx))) stop("NA in factor")
	n <- length(listw$neighbours)
	if (n != length(fx)) stop("objects of different length")
	dums <- lm(codes(fx) ~ fx - 1, x=TRUE)$x
	nc <- ncol(dums)
	res <- matrix(0, nrow=nsim+1, ncol=nc)
	res[nsim+1,] <- 0.5 * joincount(dums, listw)
	for (i in 1:nsim) {
		fxi <- sample(fx)
		dums <- lm(codes(fxi) ~ fxi - 1, x=TRUE)$x
		res[i,] <- 0.5 * joincount(dums, listw)
	}
	rankres <- apply(res, 2, rank)
	xrank <- rankres[nrow(rankres),]
	lres <- list(res=res, rankres=rankres, xrank=xrank)
	thiscall <- match.call()
	attr(lres, "call") <- thiscall
	neighbours.attrs <- names(attributes(listw$neighbours))
	attr(lres, "neighbours.attrs") <- neighbours.attrs
	weights.attrs <- names(attributes(listw$weights))
	attr(lres, "weights.attrs") <- weights.attrs
	cat("Join count permutation test for", deparse(substitute(fx)),
		"with call:\n")
	print(thiscall)
	cat("neighbours list attributes:", neighbours.attrs, "\n")
	cat("and weights list attributes:", weights.attrs, "\n\n")
	tres <- matrix(0, nrow=nc, ncol=4)
	tab <- table(fx)
	rownames(tres) <- names(tab)
	colnames(tres) <- c("Same colour statistic", "Rank", "Expectation",
		 "Variance")
	tres[,1] <- res[nsim+1,]
	tres[,2] <- xrank
	tres[,3] <- apply(res, 2, mean)
	tres[,4] <- apply(res, 2, var)
	print.coefmat(tres, signif.stars = FALSE)
	cat("\n")
	invisible(lres)
}


