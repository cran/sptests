lm.LMtests <- function(model, listw, zero.policy=FALSE) {
	if (class(listw) != "listw") stop(paste(deparse(substitute(listw)),
		"is not a listw object"))
	if(class(model) != "lm") stop(paste(deparse(substitute(model)),
		"not an lm object"))
	N <- length(listw$neighbours)
	if(is.null(model$y)) stop("Rerun lm with y=T, x=T")
	if (N != length(model$y)) stop("objects of different length")
	if (is.null(attr(listw$weights, "W")) || !attr(listw$weights, "W"))
		warning("Spatial weights matrix not row standardized")

	tracew <- function (listw) {
		dlmtr <- 0
		n <- length(listw$neighbours)
		for (i in 1:n) {
			dij <- listw$neighbours[[i]]
			ndij <- length(dij)
			wdij <- listw$weights[[i]]
			for (j in 1:ndij) {
				k <- dij[j]
				if (k > i) {
				    dk <- which(listw$neighbours[[k]] == i)
				    if (dk > 0 &&
					dk <= length(listw$neighbours[[k]]))
					wdk <- listw$weights[[k]][dk]
					else wdk <- 0
					dlmtr <- dlmtr + (wdk * wdk) + 2 *
					(wdij[j] * wdk) + (wdij[j] * wdij[j])
				}
			}
		}
		dlmtr
	}

	y <- model$y
	X <- model$x
	u <- as.vector(resid(model))
	yhat <- as.vector(fitted(model))
	p <- model$rank
	p1 <- 1:p
	XtXinv <- chol2inv(model$qr$qr[p1, p1, drop = FALSE])
	sigma2 <- (t(u) %*% u) / N
	T <- tracew(listw)
	Wu <- lag.listw(listw, u, zero.policy)
	Wy <- lag.listw(listw, y, zero.policy)
	Wyhat <- lag.listw(listw, yhat, zero.policy)
	XtWyhat <- t(X) %*% Wyhat
	dutWu <- (t(u) %*% Wu) / sigma2
	resa <- (dutWu ^ 2) / T
	result <- matrix(nrow=5,ncol=3)
	row.names(result) <- c("LMerr", "LMlag", "RLMerr", "RLMlag", "SARMA")
	colnames(result) <- c("statistic", "df", "Pr()")
	result[1,1] <- resa
	result[1,2] <- 1
	result[1,3] <- 1-pchisq(resa, 1)
	J <- (1/(N*sigma2)) *
		((t(Wyhat) %*% Wyhat) -
		(t(XtWyhat) %*% XtXinv %*% XtWyhat) +
		(T * sigma2))
	dutWy <- (t(u) %*% Wy) / sigma2
	res <- ((dutWu - (T*((N*J)^-1))*dutWy)^2) /
		(T * (1 - T*((N*J)^-1)))
	result[3,1] <- res
	result[3,2] <- 1
	result[3,3] <- 1-pchisq(res, 1)
	res <- (dutWy ^ 2) / (N * J)
	result[2,1] <- res
	result[2,2] <- 1
	result[2,3] <- 1-pchisq(res, 1)
	res <- ((dutWy - dutWu)^2)/ ((N*J) - T)
	result[4,1] <- res
	result[4,2] <- 1
	result[4,3] <- 1-pchisq(res, 1)
	res <- res + resa
	result[5,1] <- res
	result[5,2] <- 2
	result[5,3] <- 1-pchisq(res, 2)
	cat("\nCall:\n")
	print(model$call)
	cat("\nLagrange multiplier diagnostics for spatial dependence:\n")
	print.coefmat(result, signif.stars = FALSE)
	cat("\n")
	invisible(result)
}
