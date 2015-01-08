########################
## Benjamin Haibe-Kains
## All rights Reserved
## December 9, 2013
########################


## Normalize rows of a matrix to have the same quantiles, allowing for missing values.
## inpired from Gordon Smyth (normalizaeQuantiles from the limma package)
##
## input:
##  A: numeric matrix. Missing values are allowed.
##  ties: logical. If ‘TRUE’, ties in each row of ‘A’ are treated in careful way. tied values will be normalized to the mean of the corresponding pooled quantiles.
## normvector: numeric vector of values corresponding to the quantiles distribution to fit. Note that names of 'normvector' should be properly set to the colnames of 'A'
##
## output: quantil normalized matrix
`normQuant` <- 
function(A, ties=TRUE, normvector) {

  if(!missing(normvector)) {
    normvector <- sort(normvector, method="quick")
    if(length(normvector) != ncol(A)) {
      stop("length of normvector must be equal to the number of columns of A")
    }
  }
  
	n <- dim(A)
	if (is.null(n)) { return(A) }
	if (n[1] == 1) { return(A) }
	O <- S <- array( , n)
	nobs <- rep(n[2], n[1])
	i <- (0:(n[2]-1)) / (n[2]-1)
	for (j in 1:n[1]) {
		Si <- sort(A[j, ], method="quick", index.return=TRUE)
		nobsj <- length(Si$x)
		if(nobsj < n[2]) {
			nobs[j] <- nobsj
			isna <- is.na(A[j, ])
			S[j, ] <- approx((0:(nobsj-1))/(nobsj-1), Si$x, i, ties="ordered")$y
			O[j, !isna] <- ((1:n[2])[!isna])[Si$ix]
		} else {
			S[j, ] <- Si$x
			O[j, ] <- Si$ix
		}
	}
	if(!missing(normvector)) {
    m <- normvector
  } else {
    m <- colMeans(S)
  }
	for (j in 1:n[1]) {
		if (ties) { r <- rank(A[j, ]) }
		if (nobs[j] < n[2]) {
			isna <- is.na(A[j, ])
			if (ties) {
				A[j, !isna] <- approx(i, m, (r[!isna]-1)/(nobs[j]-1), ties="ordered")$y
      } else {
				A[j, O[j, !isna]] <- approx(i, m, (0:(nobs[j]-1))/(nobs[j]-1), ties="ordered")$y
      }
		} else {
			if(ties)
				A[j, ] <- approx(i, m, (r-1)/(n[2]-1), ties="ordered")$y
			else
				A[j, O[j, ]] <- m
		}
	}
	return(A)
}