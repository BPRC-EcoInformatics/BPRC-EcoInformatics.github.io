# Parametrization of Neural Networks 

# Author  :   GwénaE Leday
# Updated :   June 2009


param.nnet <- function(x,
                      y,
                      nnet.rep = 100,
                      s = 2:5,
                      m = c(500,750,1000,1250,1500),
                      d = c(0.01,0.001,0.0001),
                      resampling = "boot",
                      nrep = 200,
                      train.frac = .75,
                      kfold = 10,
                      file=T,
                      out.object = T){


	#-------------------- Initialization ------------------#
	## Data
	data <- cbind(x,y)
	ncol <- ncol(data)
	colnames(data)[ncol] <- "pred"

	## Number of parameters
	n.s  <- length(s)
	n.m  <- length(m)
	n.d  <- length(d)

	## Leaf
	leaf <- NULL
	if(resampling=="boot"){
		nrep2 <- nrep + 3
		for(i in 1:nnet.rep){leaf <- cbind(leaf,rbind(matrix(NA,nrep,2),matrix(NA,3,2)))}
		name.row.leaf <- c(paste("boot",1:nrep,sep=""),
                               c("Mean",".632 error",".632+ error"))
	}
	else{
		nrep2 <- nrep + 1
		for(i in 1:nnet.rep){leaf <- cbind(leaf,rbind(matrix(NA,nrep,2),matrix(NA,1,2)))}
		name.row.leaf <- c(paste("rep",1:nrep,sep=""),"Mean")
	}
	nb.row.leaf <- nrep2
      dimnames(leaf) <- list(name.row.leaf, rep(c("Train","Test"),nnet.rep))
	cpt <- 1

	# Tree of matrix with list and each node/level represents a parameter
	err <- list(rep(list(rep(list(rep(list(leaf),n.d)),n.m)),n.s))

	name.error <- ifelse(resampling=="cv","CV error", 
                           ifelse(resampling=="boot",".632+ error","holdout error"))

	# Results
	results <- matrix(NA,(n.s*n.m*n.d),4)
	rownames(results) <- paste("model ",1:(n.s*n.m*n.d),sep="")
	colnames(results) <- c("Size","Maxiter","Decay",name.error)

	# Output
	true.list <- test.id.list <- bal <- NULL
	name.root <- paste(c("param","nnet", resampling), sep="", collapse="_")

#----------------------------------- Body ----------------------------------#
	# Current CPU Time
	ptm1 <- proc.time()

	for(n in 1:nrep){
		cat(ifelse(resampling=="boot","b = ","rep = "),n,"\n")
		# Generate samples
		resample <- resample(data, method=resampling, train=train.frac, cv=kfold)
		data.train <- resample$train
		if(resampling!="cv"){
			data.test <- resample$test
			bal <- c(bal, sum(data.train[,ncol])/nrow(data.train))
			test.id.list <- c(test.id.list, list(resample$test.ind))
			true.list <- c(true.list, list(data.test[,ncol]))
		}else{
			cv.groups <- resample$groups
			test.id.list <- c(test.id.list, list(unlist(cv.groups)))
			true.list <- c(true.list, list(data.train[,ncol]))
		}
		
		for(i in 1:n.s){
			for(j in 1:n.m){
				for(k in 1:n.d){
					l <- 1
					while(l<(2*nnet.rep)){
						if(resampling=="cv") err.cv <- NULL
						if(resampling!="cv") cv.groups <- NA
						for(g in 1:length(cv.groups)){
							if(resampling!="cv"){
								train <- data.train
								test  <- data.test
							}
							else{
								train <- data.train[-cv.groups[[g]],]
								test  <- data.train[cv.groups[[g]],]
							}

							# Train and test
							res <- train.and.test.model(train=train, method="NNET",
									test=test, resp.var.name="pred",
									nnet.par=c(s[i],m[j],d[k]), prob=F,sp=FALSE)

							# Compute errors
							error.train <- res$error.train
							error.test  <- res$error.test

							if(resampling=="cv"){
								err.cv <- rbind(err.cv, c(error.train, error.test))
							}
						}
						if(resampling=="cv"){
							err[[1]][[i]][[j]][[k]][n,l] <- mean(err.cv[,1])
							err[[1]][[i]][[j]][[k]][n,(l+1)] <- mean(err.cv[,2])
						}else{
							err[[1]][[i]][[j]][[k]][n,l] <- error.train
							err[[1]][[i]][[j]][[k]][n,(l+1)] <- error.test
						}

						l <- l + 2
					}
				}
			}
		}
	}
	cat("\n")
	# Summarize results: average errors or calculate .632+ error
	for(i in 1:n.s){
		for(j in 1:n.m){
			for(k in 1:n.d){
				l <- 1
				while(l<(2*nnet.rep)){
					if(resampling=="boot"){
						err[[1]][[i]][[j]][[k]][,l:(l+1)] <- errors.632plus.model(data=data,
                                                                         method="NNET",
                                                                         y = y,
                                                                         err.mat = err[[1]][[i]][[j]][[k]][,l:(l+1)],
                                                                         B = nrep,
                                                                         resp.var.name = "pred",
                                                                         nnet.par = c(s[i],m[j],d[k]))
					}else{
						err[[1]][[i]][[j]][[k]][nrep2,l:(l+1)] <- apply(err[[1]][[i]][[j]][[k]][1:nrep,l:(l+1)],2,mean)
					}
					l <- l + 2
				}
				results[cpt,4] <- mean(err[[1]][[i]][[j]][[k]][nrep2,], na.rm=T)
				results[cpt,1:3] <- c(s[i],m[j],d[k])
				cpt <- cpt + 1
			}
		}
	}


	# Current CPU Time
	ptm2 <- proc.time()
	ptm <- ptm2 - ptm1

#---------------------------------- Output ---------------------------------#
	output <- NULL
	bool <- results[,4]==min(results[,4])
	best <- results[bool,]
	if(sum(bool)==1) best <- t(as.matrix(best))
	if(nrow(best)==1){
		rownames(best) <- "Best model"
	}else{
		rownames(best) <- paste("Best model ",1:nrow(best),sep="")
	}
	output$best.models <- best
	output$all.models  <- results
	output$rep.details <- err

	# Trace
	output$trace.rep <- test.id.list

	# List of response variable vectors of non-bootstrapped observations
	output$true.list <- true.list

	# Time
	output$time <- ptm
	cat("Time:\n\n")
	print(ptm)
	cat("\n\n")

	# Write results
	if(file){
		write.table(best,file=paste(species.name,"_",name.root,"_best_",bg,nb.var,"_",reduction,".txt",sep=""))
		write.table(results,file=paste(species.name,"_",name.root,"_all_",bg,nb.var,"_",reduction,".txt",sep=""))
	}

	# Save object
	if(out.object) save(output, file = paste(species.name,"_",name.root,"_output_object_",bg,nb.var,"_",reduction,".Rdata",sep=""))

	return(output)
}

