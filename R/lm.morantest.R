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

lm.morantest <- function(model, listw, zero.policy=FALSE) {
	if (class(listw) != "listw") stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if(class(model) != "lm") stop(paste(deparse(substitute(model)),
		"not an lm object"))
	if(is.null(model$x)) stop("Rerun lm with x=T")
	N <- length(listw$neighbours)
	if (N != nrow(model$x)) stop("objects of different length")

	u <- as.vector(resid(model))
	listw.U <- listw2U(listw)

	S0 <- sum(unlist(listw.U$weights))
	S1 <- 0.5 * sum((2*unlist(listw.U$weights))^2)
	lu <- lag.listw(listw.U, u, zero.policy=zero.policy)
	I <- (N/S0) * ((t(u) %*% lu) / (t(u) %*% u))
	p <- model$rank
	p1 <- 1:p
	XtXinv <- chol2inv(model$qr$qr[p1, p1, drop = FALSE])
	X <- model$x
# Cliff/Ord 1981, p. 203
	Z <- lag.listw(listw.U, X, zero.policy=zero.policy)
	C1 <- t(X) %*% Z
	trA <- -(sum(diag(XtXinv %*% C1)))
	EI <- ((N * trA) / ((N-p) * S0))
	C2 <- t(Z) %*% Z
	C3 <- XtXinv %*% C1
	trA2 <- sum(diag(C3 %*% C3))
	trB <- sum(diag(4*(XtXinv %*% C2)))
	VI <- (((N*N)/((S0*S0)*(N-p)*(N-p+2))) *
		(S1 + 2*trA2 - trB - ((2*(trA^2))/(N-p))))
	ZI <- (I - EI) / sqrt(VI)
	PrI <- 2*(1-pnorm(abs(ZI)))
	res <- matrix(0, nrow=1, ncol=5)
	rownames(res) <- "residuals"
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
	cat("Moran's I test", "\nwith call:\n") 
	print(thiscall)
	cat("neighbours list attributes:", neighbours.attrs, "\n")
	cat("and weights list attributes:", weights.attrs, "\n\n")
	print.coefmat(res, signif.stars = FALSE)
	cat("\n")
	invisible(res)
}
