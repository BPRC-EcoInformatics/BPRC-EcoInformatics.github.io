# Test of a model on an independant test set

# Note: This function typically read results that have
# been previously written by the function 'model.test'.
# It supposes that all those files (and only them) are
# in the same directory!

# Author  :   Gwena?E Leday
# Updated :   June 2009



read.model.test <- function(path, nrep){

	# Work done in this directory
	setwd(path)

	# List files
	names.files <- list.files()
	n.files <- length(names.files)
	split.rep <- n.files/nrep
	resampling <- unlist(strsplit(names.files[1],"_"))[2]

	# Read vectors of predictions
	cat("\n")
	cat("reading...\n\n")
	for(i in 1:split.rep){
		cat(paste("Part ",i,sep=""),"\n")
		vect <- NULL
		for(n in 1:nrep){
			name.file <- paste(paste(c("pred",resampling,n,i),sep="",collapse="_"),".txt",sep="")
			temp <- as.matrix(read.table(name.file))
			vect <- cbind(vect, temp)
		}
		if(i == 1){
			vect.mean <- as.matrix(apply(vect,1,mean))
			vect.sd <- as.matrix(apply(vect,1,sd))
		}
		if(i != 1){
			vect.mean <- rbind(vect.mean, as.matrix(apply(vect,1,mean)))
			vect.sd <- rbind(vect.sd, as.matrix(apply(vect,1,sd)))
		}
	}
	return(list("mean"=vect.mean, "sd"=vect.sd))
cat("\n\n")
cat("Done \n")
}

