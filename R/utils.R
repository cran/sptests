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

spweights.constants <- function(listw) {
	if(class(listw) != "listw") stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	n <- length(listw$neighbours)
	n1 <- n - 1
	n2 <- n - 2
	n3 <- n - 3
	nn <- n*n
	S0 <- Szero(listw)
	S1 <- 0
	S2 <- 0
	for (i in 1:n) {
		ij <- listw$neighbours[[i]]
		wij <- listw$weights[[i]]
		dm0 <- 0
		dm1 <- 0
		for (j in 1:length(ij)) {
			dij <- wij[j]
			ij.j <- ij[j]
			ij.lkup <- which(listw$neighbours[[ij.j]] == i)
			if (length(ij.lkup) == 1)
				dji <- listw$weights[[ij.j]][ij.lkup]
			else dji <- 0
			dm0 <- dm0 + dij
			dm1 <- dm1 + dji
			S1 <- S1 + (dij + dji)^2
		}
		S2 <- S2 + (dm0 + dm1)^2
	}
	S1 <- S1 * 0.5
	invisible(list(n=n, n1=n1, n2=n2, n3=n3, nn=nn, S0=S0, S1=S1, S2=S2))
}

Szero <- function(listw) {
	sum(unlist(listw$weights))
}

#lag.listw <- function(listw, x, zero.policy=FALSE) {
#	if (class(listw) != "listw") stop(paste(deparse(substitute(listw)),
#		"is not a listw object"))
#	n <- length(listw$neighbours)
#	if (length(x) != n) stop("object lengths differ")
#	if (any(is.na(x))) stop("NA in X")
#	cardnb <- card(listw$neighbours)
#	res <- as.numeric(rep(0, n))
#	for (i in 1:n) {
#		if (cardnb[i] == 0) {
#			if (zero.policy) res[i] <- 0
#			else res[i] <- NA
#		} else {
#			res[i] <- sum(x[listw$neighbours[[i]]] *
#				listw$weights[[i]])
#		}
#	}
#	res
#}

#lag.listw <- function(listw, x, zero.policy=FALSE) {
#	if (class(listw) != "listw") stop(paste(deparse(substitute(listw)),
#		"is not a listw object"))
#	n <- length(listw$neighbours)
#	if (length(x) != n) stop("object lengths differ")
#	if (any(is.na(x))) stop("NA in X")
#	cardnb <- card(listw$neighbours)
#	res <- .Call("lagw", listw$neighbours, listw$weights,
#		as.numeric(x), as.integer(cardnb), as.logical(zero.policy))
#	invisible(res)
#}

lag.listw <- function(listw, x, zero.policy=FALSE) {
	if (class(listw) != "listw") stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if (!is.numeric(x)) stop(paste(deparse(substitute(listw)),
		"not numeric"))
	if (any(is.na(x))) stop("NA in X")
	n <- length(listw$neighbours)
	cardnb <- card(listw$neighbours)
	if (is.vector(x)) {
		if (length(x) != n) stop("object lengths differ")
		res <- .Call("lagw", listw$neighbours, listw$weights,
			as.numeric(x), as.integer(cardnb),
			as.logical(zero.policy))
	} else if (is.matrix(x)) {
		if (nrow(x) != n) stop("object lengths differ")
		res <- matrix(0, nrow=nrow(x), ncol=ncol(x))
		for (i in 1:ncol(x)) {
			res[,i] <- .Call("lagw", listw$neighbours,
				listw$weights, as.numeric(x[,i]),
				as.integer(cardnb), as.logical(zero.policy))

		}
	} else {
		stop(paste(deparse(substitute(x)),
			"neither a numeric vector or matrix"))
	}
	if (any(is.na(res))) warning("NAs in lagged values")
	invisible(res)
}

listw2U <- function(listw) {
	if (class(listw) != "listw") stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	nb <- listw$neighbours
	wts <- listw$weights
	style <- paste(listw$style, "U", sep="")
	sym <- is.symmetric.nb(nb, FALSE, TRUE)
	n <- length(listw$neighbours)
	cardnb <- card(listw$neighbours)
	nlist <- vector(mode="list", length=n)
	attr(nlist, "region.id") <- attr(nb, "region.id")
	class(nlist) <- "nb"
	vlist <- vector(mode="list", length=n)
	attr(vlist, as.character(style)) <- TRUE
	if (sym) {
		nlist <- vector(mode="list", length=n)
		attr(nlist, "region.id") <- attr(nb, "region.id")
		class(nlist) <- "nb"
		for (i in 1:n) {
			inb <- nb[[i]]
			nlist[[i]] <- inb
			iwt <- wts[[i]]
			icd <- cardnb[i]
			for (j in 1:icd) {
				vlist[[i]][j] <- 0.5 *
				(iwt[j]+wts[[inb[j]]][which(nb[[inb[j]]] == i)])
			}
		}
	} else {
		nlist <- make.sym.nb(nb)
		for (i in 1:n) {
			inb <- nb[[i]]
			inl <- nlist[[i]]
			iwt <- wts[[i]]
			vlist[[i]] <- numeric(length=length(inl))
			for (j in 1:length(inl)) {
				if (inl[j] %in% inb) a <-
					iwt[which(inb == inl[j])]
				else a <- 0
				if (i %in% nb[[inl[j]]]) b <-
					wts[[inl[j]]][which(nb[[inl[j]]] == i)]
				else b <- 0
				vlist[[i]][j] <- 0.5 * (a + b)
			}
		}
	}
	res <- list(style=style, neighbours=nlist, weights=vlist)
	class(res) <- "listw"
	attr(res, "region.id") <- attr(nb, "region.id")
	attr(res, "call") <- match.call()
	attr(res, "U") <- TRUE
	invisible(res)
}


