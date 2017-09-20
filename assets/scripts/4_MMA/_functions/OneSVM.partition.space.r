# Partitioning of the space of parameter for one-class SVM

# Note: it allows parallel computation

# Author  :   Gwénaël Leday
# Updated :   March 2009



OneSVM.partition.space <- function(x=4, file=T){

	gam <- seq(0.01,0.99,0.01)
	nu <- seq(0.01,0.99,0.01)
	n <- length(nu)

	## all possibilities
	cpt <- 1
	mat <- matrix(NA,n*n,2)
	for(i in 1:n){
		for(j in 1:n){
			mat[cpt,] <- c(gam[i],nu[j])
			cpt <- cpt + 1
		}
	}
	colnames(mat) <- c("gam","nu")

	## Random subspace creation
	groups <- cv.samples(n*n, x)

	## Create matrices of each subspaces
	matrices <- vector("list",x)
	for(k in 1:length(groups)){
		matrices[[k]] <- mat[groups[[k]],]
	}
	if(file){
		for(i in 1:x){
			name.file <- paste(paste("OneSVM_space_parameter_",i,sep=""),".txt",sep="")
			write.table(matrices[[i]],name.file, row.names=F)
		}
	}
	return(matrices)
}