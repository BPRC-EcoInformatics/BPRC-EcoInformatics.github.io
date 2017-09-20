# Cross-Validated samples

# Author  :   Gwénaël Leday
# Updated :   May 2009



cv.samples <- function(n, k){

	k <- trunc(k)

	# Create samples
	if(k==n){
		groups <- 1:n
		leave.out <- 1
	}
	if(k<n){
		leave.out <- trunc(n/k);
		o <- sample(1:n)
		groups <- vector("list",k)
		for(j in 1:(k-1)){
			jj <- (1+(j-1)*leave.out)
			groups[[j]] <- (o[jj:(jj+leave.out-1)])
		}
		groups[[k]] <- o[(1+(k-1)*leave.out):n]
	}
	return(groups)
}