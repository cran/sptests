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

lag.listw <- function(listw, x, zero.policy=FALSE) {
	if (class(listw) != "listw") stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	n <- length(listw$neighbours)
	if (length(x) != n) stop("object lengths differ")
	if (any(is.na(x))) stop("NA in X")
	cardnb <- card(listw$neighbours)
	res <- as.numeric(rep(0, n))
	for (i in 1:n) {
		if (cardnb[i] == 0) {
			if (zero.policy) res[i] <- 0
			else res[i] <- NA
		} else {
			res[i] <- sum(x[listw$neighbours[[i]]] *
				listw$weights[[i]])
		}
	}
	res
}


