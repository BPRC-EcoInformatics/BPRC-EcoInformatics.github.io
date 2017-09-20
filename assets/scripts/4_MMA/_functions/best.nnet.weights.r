# Gets best set of weights for nnet

# Author  :   Tak Ikeda
# Updated :   January 2010



best.nnet.weights <- function(sp.name,
				data.x,
				data.y,
				nnet.par,
				nreps = 1000,
				save.out)     {


	# Arguments
	#	sp.name	 = character
	#	data.x	 = x variables
	#	data.y	 = response variable
	#	nnet.par = vector of 3, parameters for nnet
	#	nreps	 = integer, how many models to test
	#	save.out = logical, T to save weights in a file


	# initialization
	ptm1 <- proc.time()
	form <- as.formula("pred ~ .")
	train <- cbind(as.data.frame(data.x), as.factor(data.y))
	colnames(train)[ncol(train)] <- "pred"
	model <- nnet(form,data=as.data.frame(train), trace=F,
			size=nnet.par[1], maxit=nnet.par[2], decay=nnet.par[3])
	nnet.len <- length(model$wts)
	nnet.Wts <- matrix(NA,nreps,nnet.len)

	# compute weights
	for (i in 1:nreps) {
		model <- nnet(form,data=as.data.frame(train), trace=F,
				size=nnet.par[1], maxit=nnet.par[2], decay=nnet.par[3])
		nnet.Wts[i,] <- model$wts 
		if ((i/(0.1*nreps) - trunc(i/(0.1*nreps))) == 0) cat("model",i,"\n") 
	}

	# output
	av.Wts <- apply(nnet.Wts,2,mean)
	n.lab <- paste(nnet.par[1],nnet.par[2],nnet.par[3],sep="_")
	if (save.out) write.table(av.Wts,paste(sp.name,"_best_weights_",n.lab,"_",bg,nb.var,"_",reduction,".txt",sep=""),row.names=F)

	# Computation time
	ptm2 <- proc.time()
	ptm <- ptm2 - ptm1

	# Time
	cat("Time:\n")
	print(ptm)
	cat("\n")


	return(av.Wts)

}

