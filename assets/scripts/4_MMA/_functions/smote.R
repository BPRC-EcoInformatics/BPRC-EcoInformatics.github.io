# Synthetic Minority Over-sampling TEchnique (SMOTE)
# Implement the numerical version of the algorithm

# Author  :   Takayoshi Ikeda
# Updated :   June 2009


# T: database of the minority
# N: number of synthetic instances to create
# k: knn
smote <- function(T,N,k) {

	## Initialization
	T <- t(T)
	NT <- ncol(T)
	if (NT==0) print("Check T.")
	if (NT==1) sampleset <- matrix((rep(T,2)),ncol=2)

	if (NT>1) {
		# number of nearest neighbours cannot be greater than NT-1
		if (k>(NT-1)) {
			k <- NT-1
			warning("Not so many instances in T.k is set to ", k)
		}
		# number of new examples that each example in T should generate
		NumAtt <- nrow(T)
		n <- floor(N/NT)
		remainder <- N-NT*n
		id <- sample(NT)
		No <- rep(1,NT)*n
		No[id[1:remainder]] <- No[id[1:remainder]]+1
		
		## k-NN generation
		sampleset <- NULL
		for (i in which(No!=0)) {
			d <- dist_nn(t(T[,i]),T)
			d[i] <- Inf
			if (k < log(NT)) {
				min_id <- NULL
					for (j in 1:k) {
						tmp <- min(d)
						id <- which.min(d)
						d[id] <- Inf
						min_id <- c(min_id,id)
					}
			}
			if (k > log(NT)) {
				tmp <- sort(d)
				id <- order(d)
				min_id <- id[1:k]
			}
		
			rn <- floor(runif(No[i])*k)+1
			id <- min_id[rn]
			weight <- matrix(runif(NumAtt*No[i]),NumAtt,No[i])
			D2 <- matrix(rep(T[,i],No[i]),nrow=nrow(T))
			# for numeric attributes
			D2 <- D2 + weight*(T[,id]-D2)
			sampleset <- cbind(sampleset,D2)
		}	
	}
	return(t(sampleset))
}


# 'dist' function in Matlab calculating Euclidean distance weight function.
# in the nnet toolbox.
dist_nn <- function(w,p,fp) {

	S <- nrow(w)
	R <- ncol(w)
	R2 <- nrow(p)
	Q <- ncol(p)

	if (R != R2) print("Inner matrix dimensions do not match.")
	z <- matrix(0,S,Q)
	if (Q<S) {
 		p <- t(p)
 		for (q2 in 1:Q) z[,q2] <- rowSums((w-p[q2,])^2)
  	}
	if (Q>=S) {
  		w <- t(w)
  		for (i in 1:S) z[i,] <- colSums((w[,i]-p)^2)
  	}
	return(sqrt(z))
}

