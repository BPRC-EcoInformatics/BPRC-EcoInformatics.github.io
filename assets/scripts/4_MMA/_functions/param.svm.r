# Parametrization of Support Vector Machines (SVM)

# Author  :   Gwenael Leday
# Updated :   June 2009


param.svm <- function(x,
                      y,
                      kernL = c("vanilladot", "rbfdot", "polydot"),
                      par.c = c(0.1,1,10,100),
                      par.gam = c(0.01,0.1,0.3,0.5,0.7,0.9),
                      par.deg = 2:4,
                      par.coef = c(0,1,4,8),
                      resampling = "boot",
                      nrep = 200,
                      train.frac = .75,
                      kfold = 10,
                      file=T,
                      out.object = T,
		      	  thresh = NA,
				  prev = NA){


	#-------------------- Initialization ------------------#
	## Data
	data <- cbind(x,y)
	ncol <- ncol(data)
	colnames(data)[ncol] <- "pred"

	## Number of parameters
	n.c    <- length(par.c)
	n.gam  <- length(par.gam)
	n.deg  <- length(par.deg)
	n.coef <- length(par.coef)

	## Leaf
	if(resampling=="boot"){
		nb.row.leaf <- nrep+3
		name.row.leaf <- c(paste("boot",1:nrep,sep=""),
                               c("Mean",".632 error",".632+ error"))
	}
	else{
		nb.row.leaf <- nrep+1
		name.row.leaf <- c(paste("rep",1:nrep,sep=""),"Mean")
	}
	leaf <- matrix(NA, nb.row.leaf, 2,
                      dimnames=list(name.row.leaf, c("Train","Test")))

	## Tree storage
	names.lin <- names.rb <- names.poly <- NULL
	nb.models <- 0
	if("vanilladot"%in%kernL){
		error.lin <- list(rep(list(leaf),n.c))
		names.lin <- rep("vanilladot",n.c)
		nb.models <- nb.models + n.c
	}
	if("rbfdot"%in%kernL){
		error.rb <- list(rep(list(rep(list(leaf),n.gam)),n.c))
		names.rb <- rep("rbfdot",(n.c*n.gam))
		nb.models <- nb.models + n.c*n.gam
	}
	if("polydot"%in%kernL){
		error.poly <- list(rep(list(rep(list(rep(list(rep(list(leaf),n.coef)),
					n.deg)),n.gam)),n.c))
		names.poly <- rep("polydot",(n.c*n.gam*n.deg*n.coef))
		nb.models <- nb.models + n.c*n.gam*n.deg*n.coef
	}
	name.error <- ifelse(resampling=="cv","CV error", 
                           ifelse(resampling=="boot",".632+ error","holdout error"))

	# Best error
	best <- matrix(NA,nb.models,5)
	rownames(best) <- c(names.lin, names.rb, names.poly)
	colnames(best) <- c("C","gamma","degree","coef0",name.error)

	# Output
	true.list <- test.id.list <- bal <- NULL
	name.root <- paste(c("param","svm", resampling), sep="", collapse="_")

	#------------------------ Body ------------------------#
	# Current CPU Time
	ptm1 <- proc.time()

	# Evaluation of each model by bootstrapping
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
		for(kern in kernL){
			switch(kern,
				"vanilladot"  = {
					len.i <- n.c
					len.j <- len.k <- len.l <- 1
					p.c    <- par.c
					p.gam  <- p.deg  <- p.coef <- NULL
				},
				"rbfdot"  = {
					len.i <- n.c
					len.j <- n.gam
					len.k <- len.l <- 1
					p.c    <- par.c
					p.gam  <- par.gam
					p.deg  <- p.coef <- NULL
				},
				"polydot"  = {
					len.i <- n.c
					len.j <- n.gam
					len.k <- n.deg
					len.l <- n.coef
					p.c    <- par.c
					p.gam  <- par.gam
					p.deg  <- par.deg
					p.coef <- par.coef
				}
			)
			for(i in 1:len.i){
				for(j in 1:len.j){
					for(k in 1:len.k){
						for(l in 1:len.l){
							if(resampling=="cv") err.cv <- NULL
							if(resampling!="cv") cv.groups <- NA
							for(m in 1:length(cv.groups)){
								if(resampling!="cv"){
									train <- data.train
									test  <- data.test
								}else{
									train <- data.train[-cv.groups[[m]],]
									test  <- data.train[cv.groups[[m]],]
								}
#cat("\n\n")
#cat("Train ",kern, ",\t c = ",p.c[i],"\t gam = ",p.gam[j],"\t deg = ", p.deg[k], "\t coef = ", p.coef[l],"\n")
								# Train and test
								res <- train.and.test.model(train=train, method="SVM",
									test=test, resp.var.name="pred",
									svm.par=c(kern, p.c[i], p.gam[j], p.deg[k], p.coef[l]), prob=F,sp=FALSE)

								# Errors
								err.train <- res$error.train
								err.test  <- res$error.test

								if(resampling=="cv"){
									err.cv <- rbind(err.cv, c(err.train, err.test))
								}
							}
							if(resampling=="cv"){
								if(kern=="vanilladot") error.lin[[1]][[i]][n,1:2] <- c(mean(err.cv[,1]),mean(err.cv[,2]))
								if(kern=="rbfdot")     error.rb[[1]][[i]][[j]][n,1:2] <- c(mean(err.cv[,1]),mean(err.cv[,2]))
								if(kern=="polydot")    error.poly[[1]][[i]][[j]][[k]][[l]][n,1:2] <- c(mean(err.cv[,1]),mean(err.cv[,2]))
							}else{
								if(kern=="vanilladot") error.lin[[1]][[i]][n,1:2] <- c(err.train,err.test)
								if(kern=="rbfdot")     error.rb[[1]][[i]][[j]][n,1:2] <- c(err.train,err.test)
								if(kern=="polydot")    error.poly[[1]][[i]][[j]][[k]][[l]][n,1:2] <- c(err.train,err.test)
							}
						}
					}
				}
			}
		}
	}

	# Determination of model that minimize error
	for(kern in kernL){
		switch(kern,
			"vanilladot"  = {
				len.i <- n.c
				len.j <- len.k <- len.l <- 1
				p.c    <- par.c
				p.gam  <- p.deg  <- p.coef <- NULL
			},
			"rbfdot"  = {
				len.i <- n.c
				len.j <- n.gam
				len.k <- len.l <- 1
				p.c    <- par.c
				p.gam  <- par.gam
				p.deg  <- p.coef <- NULL
			},
			"polydot"  = {
				len.i <- n.c
				len.j <- n.gam
				len.k <- n.deg
				len.l <- n.coef
				p.c    <- par.c
				p.gam  <- par.gam
				p.deg  <- par.deg
				p.coef <- par.coef
			}
		)
		cpt1 <- 1
		for(i in 1:len.i){
			for(j in 1:len.j){
				for(k in 1:len.k){
					for(l in 1:len.l){
						switch(kern,
							"vanilladot"  = {
								err <- error.lin[[1]][[i]]
							},
							"rbfdot"  = {
								err <- error.rb[[1]][[i]][[j]]
							},
							"polydot"  = {
								err <- error.poly[[1]][[i]][[j]][[k]][[l]]
							}
						)
						if(resampling=="boot"){
							err <- errors.632plus.model(data=data,
                                                               method="SVM",
                                                               y = y,
                                                               err.mat = err,
                                                               B = nrep,
                                                               resp.var.name = "pred",
                                                               svm.par = c(kern, p.c[i], p.gam[j], p.deg[k], p.coef[l]))
							nrep2 <- nrep+3
						}else{
							err[(nrep+1),1:2] <- apply(err[1:nrep,],2,mean)
							nrep2 <- nrep+1
						}
						switch(kern,
							"vanilladot"  = {
								error.lin[[1]][[i]] <- err
							},
							"rbfdot"  = {
								error.rb[[1]][[i]][[j]] <- err
							},
							"polydot"  = {
								error.poly[[1]][[i]][[j]][[k]][[l]] <- err
							}
						)
						best[rownames(best)==kern,][cpt1,5] <- err[nrep2,1]
						best[rownames(best)==kern,][cpt1,1:4] <- c(p.c[i], p.gam[j], p.deg[k], p.coef[l])
						cpt1 <- cpt1 + 1
					}
				}
			}
		}
	}

	# Current CPU Time
	ptm2 <- proc.time()
	ptm <- ptm2 - ptm1

	#---------------------- Output ---------------------#
	output <- NULL
	best.linear <- best.radial <- best.poly <- NULL
	names.lin <- names.rb <- names.poly <- NULL

	## Search best models for each kernel
	if("vanilladot"%in%kernL){
		xx <- best[rownames(best)=="vanilladot",1:5]
		best.linear <- xx[xx[,5]==min(xx[,5]),]
		nb.best.lin <- ifelse(is.null(nrow(best.linear)),1,nrow(best.linear))
		names.lin <- rep("vanilladot",nb.best.lin)
		output$errors.vanilladot <- error.lin
	}
	if("rbfdot"%in%kernL){
		xx <- best[rownames(best)=="rbfdot",1:5]
		best.radial <- xx[xx[,5]==min(xx[,5]),]
		nb.best.rb <- ifelse(is.null(nrow(best.radial)),1,nrow(best.radial))
		names.rb <- rep("rbfdot",nb.best.rb)
		output$errors.rbfdot <- error.rb
	}
	if("polydot"%in%kernL){
		xx <- best[rownames(best)=="polydot",1:5]
		best.poly <- xx[xx[,5]==min(xx[,5]),]
		nb.best.poly <- ifelse(is.null(nrow(best.poly)),1,nrow(best.poly))
		names.poly <- rep("polydot",nb.best.poly)
		output$errors.polydot <- error.poly
	}

	## Output object
	best.mod <- rbind(best.linear,best.radial,best.poly)
	rownames(best.mod) <- c( names.lin, names.rb, names.poly)
	colnames(best.mod) <- c("C","gamma","degree","coef0",name.error)
	output$best.models <- best.mod
	output$all.models <- best

	# Trace
	output$trace.rep <- test.id.list

	# List of response variable vectors of non-bootstrapped observations
	output$true.list <- true.list

	# Time
	output$time <- ptm
	cat("Time:\n\n")
	print(ptm)
	cat("\n\n")

	if(file){
		if (is.na(thresh) & is.na(prev)) {
			write.table(best.mod,paste(species.name,"_",name.root,"_best_models_",bg,nb.var,"_",reduction,".txt",sep=""))
			write.table(best,paste(species.name,"_",name.root,"_all_models_",bg,nb.var,"_",reduction,".txt",sep=""))
		}
		if (!is.na(thresh) | !is.na(prev)) {
			write.table(best.mod,paste(species.name,"_",name.root,"_best_models_thresh_",thresh,"_prev_",prev,"_",bg,nb.var,"_",reduction,".txt",sep=""))
			write.table(best,paste(species.name,"_",name.root,"_all_models_thresh_",thresh,"_prev_",prev,"_",bg,nb.var,"_",reduction,".txt",sep=""))
		}
	}
	if(out.object) {
		if (is.na(thresh) & is.na(prev)) save(output, file = paste(species.name,"_",name.root,"_output_object_",bg,nb.var,"_",reduction,".Rdata",sep=""))
		if (!is.na(thresh) | !is.na(prev)) save(output, file = paste(species.name,"_",name.root,"_output_object_thresh_",thresh,"_prev_",prev,"_",bg,nb.var,"_",reduction,".Rdata",sep=""))
	}
	return(output)
}
