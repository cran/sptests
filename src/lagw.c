/* Copyright 2001 by Roger S. Bivand. 
*
* This program is free software; you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation; either version 2 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
**/

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Applic.h>
#define ROFFSET 1

SEXP lagw(SEXP nb, SEXP weights, SEXP x, SEXP card, SEXP zeropolicy) {
	int i, j, k, n=length(card), pc=0;
	double sum, wt;
	SEXP ans;
	PROTECT(ans = NEW_NUMERIC(n)); pc++;

	for (i=0; i < n; i++) {
	    if (INTEGER_POINTER(card)[i] == 0) {
		if (LOGICAL_POINTER(zeropolicy)[0] == TRUE)
		    NUMERIC_POINTER(ans)[i] = 0;
		else
		    NUMERIC_POINTER(ans)[i] = NA_REAL;
	    }
	    else {
		sum = 0.0;
		for (j=0; j<INTEGER_POINTER(card)[i]; j++) {
		    k = INTEGER_POINTER(VECTOR_ELT(nb, i))[j];
		    wt = NUMERIC_POINTER(VECTOR_ELT(weights, i))[j];
		    sum += NUMERIC_POINTER(x)[k-ROFFSET] * wt;
		}
		NUMERIC_POINTER(ans)[i] = sum;
	    }
        }

	UNPROTECT(pc); /* ans */
	return(ans);
}

