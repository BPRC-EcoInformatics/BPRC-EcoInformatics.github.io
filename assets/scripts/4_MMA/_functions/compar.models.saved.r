# Comparison of supervised models using different performance metrics

# Author  :   Gwenael Leday
# Updated :   October 2009
# Edited: By TI, July 2012, added F-score beta option


compar.models.saved <- function(x, 
                          y,
                          nrep=200, 
                          methods=c("LDA","QDA","LOG","NB","CART",
                                    "CTREE","KNN","SVM","NNET"),
                          bagging  = NULL,
                          boosting = NULL,
                          knn.par  = NULL,
                          svm.par  = NULL, 
                          nnet.par = NULL,
                          nbtree   = 10000,
				  weights = NULL,
				  Wts  = NULL,
				  beta=1,
                          resampling="boot",
                          train.frac=0.75,
                          k=10,
                          file=FALSE, 
                          plots=FALSE,
                          plots.col=c("blue", "deepskyblue", "chocolate1",
                                      "magenta3", "green3", "chartreuse4",
                                      "darkgoldenrod1", "red", "black"),
                          plots.all=FALSE, 
                          boxplots=FALSE,
                          box.col=NULL,
                          out.object=T,sp=FALSE,coord.spind_order=NULL){

	#------------------------------ Error handling -----------------------------#
	# Check x and y
	if(is.null(dim(x))) stop("Error: dimensions of data frame 'x'")
	if(!is.null(dim(y))) stop("Error: dimensions of vector 'y'")
	if(dim(x)[1]!=length(y)) stop("Error: compatibility between 'x' and 'y'")

	# Check nrep
	if(nrep<2) stop("A minimum of two repetitions is required")

	# Check methods
	bool <- methods %in% c("LDA","QDA","LOG","NB","CART","STUMP",
                                    "CTREE","KNN","SVM","NNET","RF")

	if(sum(bool)!=length(methods)) stop("Error: do not recognize name of method")

	# Check KNN parameter
	if(is.null(knn.par) & "KNN"%in%methods){
		stop("Error: Which value for the parameter knn.par ?")}

	# Check SVM parameters
	if(is.null(svm.par) & "SVM"%in%methods) stop("Error: parameters of svm")
	if(!is.null(svm.par)){
		if(!is.vector(svm.par)) stop("Error: parameters of svm")
		if(!(svm.par[1] %in% c("vanilladot","rbfdot","polydot",
	                             "linear","radial","polynomial"))) {
			stop("Error: input kernel")
		}
		if(svm.par[1]=="linear") svm.par[1] <- "vanilladot"
		if(svm.par[1]=="radial") svm.par[1] <- "rbfdot"
		if(svm.par[1]=="polynomial") svm.par[1] <- "polydot"
		if(svm.par[1]=="vanilladot" && length(svm.par)>2) {
			stop("Error: too much parameters for the linear kernel")
		}
		if(svm.par[1]=="rbfdot" && length(svm.par)>3) {
			stop("Error: too much parameters for the radial basis kernel")
		}
		if(svm.par[1]=="polydot" && length(svm.par)>5) {
			stop("Error: too much parameters for the polynomial kernel")
		}
	}

	# Check NNET parameters
	if(is.null(nnet.par) & "NNET"%in%methods) stop("Error: parameters of nnet")

	# Check boosting
	bool <- c("KNN","SVM","RF") %in% methods
	if(!is.null(boosting) & sum(bool)!=0){
		stop("Error: boosting does not support input method(s)")
	}

	# Check nbtree, k and train.frac
	if(length(train.frac)!=1 | !is.numeric(train.frac)) stop("Error: 'train.frac' value")
	if(length(nbtree)!=1 | !is.numeric(nbtree)) stop("Error: 'nbtree' value")
	if(length(k)!=1 | !is.numeric(nbtree)) stop("Error: 'k' value")
	nbtree <- trunc(nbtree)
	k <- trunc(k)

	# Check resampling
	if(length(resampling)!=1) stop("Error: dimensions of resampling")
	if(!(resampling %in% c("holdout", "boot", "cv"))) stop("Error: 'resampling' value")

	# Check file, plots, plots.all, boxplots and out.object
	if(!is.logical(plots))     stop("Error: 'plots' has to be a logical")
	if(!is.logical(plots.all)) stop("Error: 'plots.all' has to be a logical")
	if(!is.logical(boxplots))  stop("Error: 'boxplots' has to be a logical")
	if(!is.logical(out.object))stop("Error: 'out.object' has to be a logical")

	#------------------------------ Initialization -----------------------------#
	# Number of methods
	nb.meth <- length(methods)	

	# Merge predictors and response variable
	data <- as.data.frame(cbind(x,y))
	ncol <- ncol(data)
	nrow <- nrow(data)
	colnames(data)[ncol] <- "pred"
	y <- data[,ncol]

	# Index matrix
	names.ind <- c("Accuracy","Precision","Recall","F-score","Kappa","Specificity", "TSS")
	index <- matrix(NA, nrep, length(names.ind), dimnames=list(paste("rep",1:nrep,sep=""),
                   names.ind))
	if(sp){	
	names.ind_sp <- c("sp_autocorrelation","Recall","Kappa","Specificity","TSS","AUC")
	index_sp <- matrix(NA, nrep, length(names.ind_sp), dimnames=list(paste("rep",1:nrep,sep=""),
                   names.ind_sp))}

	# Error matrix
	if(resampling=="boot"){
		nb.row.error <- nrep+3
		name.row.error <- c(paste("boot",1:nrep,sep=""),c("Mean",".632 error",".632+ error"))
	}
	else{
		nb.row.error <- nrep
		name.row.error <- paste("rep",1:nrep,sep="")
	}
	error <- matrix(NA, nb.row.error, 2,dimnames=list(name.row.error, c("Train","Test")))


	


	# Initialize probabilities, errors and indeces matrix of each method
	for(i in 1:nb.meth){
		n.prob <- paste("prob.", methods[i], sep="")
		n.err <- paste("error.", methods[i], sep="")
		n.ind <- paste("index.", methods[i], sep="")
		if(sp){ n.ind_sp <- paste("index_sp.", methods[i], sep="")}
		if(resampling=="cv") n.err.cv <- paste("error.cv.", methods[i], sep="")
		
		assign(n.prob, matrix(NA, nrow, nrep))
		
		assign(n.err, error)
		assign(n.ind, index)
		if(sp){assign(n.ind_sp, index_sp)}
		if(resampling=="cv") assign(n.err.cv, NULL)
	}
	resample1 <- resample(data, method=resampling, train=train.frac, cv=k)
	# Output
	true.list <- test.id.list <- bal <- NULL
	out <- results <- err <- uncert <- summary.index <-  pp.all <- NULL
	if(sp){summary.index_sp <- NULL}
	name.error <- ifelse(resampling=="cv","CV error", 
                           ifelse(resampling=="boot",".632+ error","holdout error"))

	# Name for files
	type.model <- ifelse(is.null(bagging), 
                           ifelse(is.null(boosting),"","boosted_"),"bagged_")
	type.resampling <- paste(resampling,"_",sep="")
	name.root <- paste(paste("compar_", type.model,      sep=""),
                         paste("models_", type.resampling, sep=""),sep="")

	#----------------------------------- Body ----------------------------------#
	# Current CPU Time
	ptm1 <- proc.time()

	## Train and test simple, bagged or boosted models
	for(n in 1:nrep){
		cat(ifelse(resampling=="boot","b = ","rep = "),n,"\n")
		# Generate samples
		resample1 <- resample(data, method=resampling, train=train.frac, cv=k)
		data.train <- resample1$train
		if(resampling!="cv"){
			data.test <- resample1$test
			bal <- c(bal, sum(data.train[,ncol])/nrow(data.train))
			test.id.list <- c(test.id.list, list(resample1$test.ind))
			true.list <- c(true.list, list(data.test[,ncol]))
		}else{
			cv.groups <- resample1$groups
			test.id.list <- c(test.id.list, list(unlist(cv.groups)))
			true.list <- c(true.list, list(data.train[,ncol]))
		}

		# Train and test each model, fill associated objects
		for(i in 1:nb.meth){
			# Names objects
			n.res  <- paste(methods[i], ".results", sep="")
			n.prob <- paste("prob."   , methods[i], sep="")
			n.err  <- paste("error."  , methods[i], sep="")
			n.ind  <- paste("index."  , methods[i], sep="")
			if(sp){n.ind_sp  <- paste("index_sp."  , methods[i], sep="")}
			if(resampling=="cv") n.err.cv <- paste("error.cv.", methods[i], sep="")
			if(resampling=="cv") n.temp   <- paste("temp.", methods[i], sep="")

			# Initialisation of inficator matrix
			#prob.mat  <- matrix(F, nrow, nrep)
			prob.mat  <- matrix(F, nrow, nrep)
			error.mat <- matrix(F, nrep+3,2)
			index.mat <- matrix(F, nrep, 5)
			if(sp){index.mat_sp <- matrix(F, nrep, 5)}
			if(resampling=="cv") assign(n.temp,NULL)

			# Set indeces
			if(resampling!="cv") prob.indeces <- which(resample1$test.ind) + nrow*(n-1)
			if(resampling=="boot") error.indeces <- c(n, n+(nrep+3))
			if(resampling!="boot") error.indeces <- c(n, n+nrep)
			index.indeces <- c(n, (1:(length(names.ind)-1)*nrep)+n)
			if(sp){index.indeces_sp <- c(n, (1:(length(names.ind_sp)-1)*nrep)+n)}

			# Set indicator matrix
			if(resampling!="cv") prob.mat  <- replace(prob.mat,  prob.indeces, T)
			if(resampling=="cv") prob.mat[,n] <- T
			error.mat <- replace(error.mat, error.indeces, T)
			index.mat <- replace(index.mat, index.indeces, T)
			if(sp){index.mat_sp <- replace(index.mat_sp, index.indeces_sp, T)}

			# Store results from training and testing
			if(resampling!="cv") cv.groups <- NA
			for(j in 1:length(cv.groups)){
				if(resampling!="cv"){
					train <- data.train
					save(train,file=paste("train_rep",n,"_",bg,nb.var,"_",reduction,".rda",sep=""))
					test  <- data.test
					save(test,file=paste("test_rep",n,"_",bg,nb.var,"_",reduction,".rda",sep=""))

				}
				else{
					train_coord <- data.train[-cv.groups[[j]],]
					save(train_coord,file=paste("train_rep",n,"_cv.group",j,"_",bg,nb.var,"_",reduction,".rda",sep=""))
					train <- data.train[-cv.groups[[j]],-(1:2)]

					test_coord  <- data.train[cv.groups[[j]],]
					save(test_coord,file=paste("test_rep",n,"_cv.group",j,"_",bg,nb.var,"_",reduction,".rda",sep=""))
					test <- test_coord [,-(1:2)]
				}
				if(is.null(bagging) & is.null(boosting)){

					#print('assign n.res')
					assign(n.res, train.and.test.model(
							     train=train,
							     method=methods[i],
							     test=test,
							     resp.var.name="pred",
							     knn.par=knn.par, svm.par=svm.par, 
                                               nnet.par=nnet.par, nbtree=nbtree,
							     Wts=Wts,
								 beta=beta,test_coord=test_coord,sp,coord.spind_order=coord.spind_order[cv.groups[[j]],]))
								 #print('get n.res')
								 result<-get(n.res)
								 save(result,file=paste("res_compar_"  , methods[i],"_rep",n,"_cv.group",j,"_",bg,nb.var,"_",reduction,".rda",sep=""))
				}
				if(!is.null(bagging)){
					assign(n.res, train.and.test.bagged.model(
							     train=train, 
							     methods=methods[i],
							     test=test, 
							     B=bagging, resp.var.name="pred",
							     knn.par=knn.par, svm.par=svm.par, 
                                               nnet.par=nnet.par, nbtree=nbtree,
							     Wts=Wts,
								 beta=beta))
				}
				if(!is.null(boosting)){
					assign(n.res, train.and.test.boosted.model(
							     train=train, 
							     methods=methods[i],
						  	     test=test, 
							     M=boosting, resp.var.name="pred",
							     knn.par=knn.par, svm.par=svm.par, 
                                               nnet.par=nnet.par, nbtree=nbtree,
							     Wts=Wts,
								 beta=beta))
				}
				if(resampling=="cv"){
					#print('get n.err')
					assign(n.err.cv, rbind(get(n.err.cv), 
                                                     c(get(n.res)$error.train,
                                                       get(n.res)$error.test)))
					#print('get n.temp')
					assign(n.temp, rbind(get(n.temp),  cbind(get(n.res)$prob.test, cv.groups[[j]])))
						
			}
			}

			# Update prob matrix of the associated method
			if(resampling!="cv") assign(n.prob,  replace(get(n.prob), 
                                              prob.mat, get(n.res)$prob.test))
			if(resampling=="cv"){
				proba.res.temp <- get(n.temp)[order(get(n.temp)[,2]),1]
				assign(n.prob,  replace(get(n.prob), prob.mat, proba.res.temp))
			}

			# Train & test error as cv average if cv
			if(resampling!="cv") trainNtest <- c(get(n.res)$error.train,get(n.res)$error.test)
			if(resampling=="cv") trainNtest <- apply(get(n.err.cv),2,mean)

			# Update error matrix of the associated method
			assign(n.err, replace(get(n.err), error.indeces, trainNtest))

			# Update index matrix of the associated method
			assign(n.ind, replace(get(n.ind), index.indeces,get(n.res)$indexes.test))
			if(sp){assign(n.ind_sp, replace(get(n.ind_sp), index.indeces_sp,get(n.res)$indexes.test.spcorrected))}
			
		}
	#print('end method loop')
	}
	#print('end rep loop')

cat("\n\n")
	## Arrange results of each method
	#print('Arrange results of each method')
	for(i in 1:nb.meth){
		# Create name of object for access
		n.err  <- paste("error."  , methods[i], sep="")
		n.ind  <- paste("index."  , methods[i], sep="")
		if(sp){n.ind_sp <- paste("index_sp."  , methods[i], sep="")}
		n.prob <- paste("prob."   , methods[i], sep="")

		# Calculate the .632+ estimate of error for each model
		if(resampling=="boot"){
			cat("Calculation .632+ error \n")
			assign(n.err, errors.632plus.model(data=data, method=methods[i], 
                                    y=y, err.mat=get(n.err), B=nrep, resp.var.name="pred",
                                    knn.par=knn.par, svm.par=svm.par, nnet.par=nnet.par, 
                                    nbtree=nbtree, bagging=bagging, boosting=boosting))
		}
		else{
			assign(n.err, rbind(get(n.err), "Mean"=apply(get(n.err),2,mean)))
		}

		# Change structure of index for output
		assign(n.ind, summary.indexes(get(n.ind), sp=FALSE,  boxplots=T))
		if(sp){assign(n.ind_sp, summary.indexes(get(n.ind_sp), sp=TRUE,  boxplots=T))}

		# Change structure of prob for output
		assign(n.prob, summary.post.prob(get(n.prob),  nrep, y))

		# Objects for output
		if(sp){results <- c(results, list(c(get(n.prob), "errors"=list(get(n.err)), get(n.ind), get(n.ind_sp))))}
		else{results <- c(results, list(c(get(n.prob), "errors"=list(get(n.err)), get(n.ind))))}
		temp    <- ifelse(resampling!="cv",get(n.err)[(nrep+3)],get(n.err)[nrep+1,2])
		err     <- rbind(err, temp)
		uncert  <- rbind(uncert, results[[i]]$uncertainty)
		summary.index <- rbind(summary.index, get(n.ind)$summary[1,])
		if(sp){summary.index_sp <- rbind(summary.index_sp, get(n.ind_sp)$summary[1,])}

		if(plots) pp.all <- c(pp.all,list(cbind(get(n.prob)$post.prob[,1], y)))

		# Write results for predictions obtained by averaging posterior probabilities
		if(file){
			name.file <- paste(paste(paste(species.name,"_",name.root,"prob_",sep=""),
                                                 methods[i],sep=""),"_",bg,nb.var,"_",reduction,".txt",sep="")
			write.table(get(n.prob)$post.prob, file=name.file,  row.names=F)
		}
	}
	#print('End Arrange results of each method')
	# Add names to output objects
	names(results) <- rownames(err) <- rownames(uncert) <- rownames(summary.index)  <- methods
	if(sp){rownames(summary.index_sp) <- methods}
	if(plots) names(pp.all) <- methods
	colnames(err) <- name.error
	colnames(uncert) <- "Uncertainty"

	# Current CPU Time
	ptm2 <- proc.time()
	ptm <- ptm2 - ptm1

cat("done\n\n")
	#---------------------------------- Output ---------------------------------#
	## All results for each classifier, i.e:
	# probabilities and sd, confusion matrix, indexes, errors
	out$results   <- results

	## General results
	# .632+ error
	if(resampling=="boot")    out$err.632plus <- err
	if(resampling=="cv")      out$err.cv      <- err
	if(resampling=="holdout") out$err.holdout <- err

	# Uncertainty
	out$uncertainty <- uncert

	# Table of ranks and score for decision
	ranks <- cbind(err, rank(err), uncert, rank(uncert), (rank(err)+rank(uncert)))
	colnames(ranks) <- c(name.error, "Rank1", "Uncertainty", "Rank2", "Score")
	out$ranking.table <- ranks

	# Averaged indexes only
	out$indexes <- summary.index
	if(sp){out$indexes_sp <- summary.index_sp}

	# Trace
	out$trace.rep <- test.id.list

	# List of response variable vectors of non-bootstrapped observations
	out$true.list <- true.list

	# Balanceness of training set
	if(resampling!="cv"){
		bal <- as.matrix(bal)
		dimnames(bal) <- list(paste("nrep",1:nrep,sep=""), c("% of 1"))
		out$balanceness <- bal
	}

	# Time
	out$time <- ptm
	cat("Time:\n")
	print(ptm)
#------------------------------ Error handling -----------------------------#
	# Check x and y
	if(is.null(dim(x))) stop("Error: dimensions of data frame 'x'")
	if(!is.null(dim(y))) stop("Error: dimensions of vector 'y'")
	if(dim(x)[1]!=length(y)) stop("Error: compatibility between 'x' and 'y'")

	# Check nrep
	if(nrep<2) stop("A minimum of two repetitions is required")

	# Check methods
	bool <- methods %in% c("LDA","QDA","LOG","NB","CART","STUMP",
                                    "CTREE","KNN","SVM","NNET","RF")

	if(sum(bool)!=length(methods)) stop("Error: do not recognize name of method")

	# Check KNN parameter
	if(is.null(knn.par) & "KNN"%in%methods){
		stop("Error: Which value for the parameter knn.par ?")}

	# Check SVM parameters
	if(is.null(svm.par) & "SVM"%in%methods) stop("Error: parameters of svm")
	if(!is.null(svm.par)){
		if(!is.vector(svm.par)) stop("Error: parameters of svm")
		if(!(svm.par[1] %in% c("vanilladot","rbfdot","polydot",
	                             "linear","radial","polynomial"))) {
			stop("Error: input kernel")
		}
		if(svm.par[1]=="linear") svm.par[1] <- "vanilladot"
		if(svm.par[1]=="radial") svm.par[1] <- "rbfdot"
		if(svm.par[1]=="polynomial") svm.par[1] <- "polydot"
		if(svm.par[1]=="vanilladot" && length(svm.par)>2) {
			stop("Error: too much parameters for the linear kernel")
		}
		if(svm.par[1]=="rbfdot" && length(svm.par)>3) {
			stop("Error: too much parameters for the radial basis kernel")
		}
		if(svm.par[1]=="polydot" && length(svm.par)>5) {
			stop("Error: too much parameters for the polynomial kernel")
		}
	}

	# Check NNET parameters
	if(is.null(nnet.par) & "NNET"%in%methods) stop("Error: parameters of nnet")

	# Check boosting
	bool <- c("KNN","SVM","RF") %in% methods
	if(!is.null(boosting) & sum(bool)!=0){
		stop("Error: boosting does not support input method(s)")
	}

	# Check nbtree, k and train.frac
	if(length(train.frac)!=1 | !is.numeric(train.frac)) stop("Error: 'train.frac' value")
	if(length(nbtree)!=1 | !is.numeric(nbtree)) stop("Error: 'nbtree' value")
	if(length(k)!=1 | !is.numeric(nbtree)) stop("Error: 'k' value")
	nbtree <- trunc(nbtree)
	k <- trunc(k)

	# Check resampling
	if(length(resampling)!=1) stop("Error: dimensions of resampling")
	if(!(resampling %in% c("holdout", "boot", "cv"))) stop("Error: 'resampling' value")

	# Check file, plots, plots.all, boxplots and out.object
	if(!is.logical(plots))     stop("Error: 'plots' has to be a logical")
	if(!is.logical(plots.all)) stop("Error: 'plots.all' has to be a logical")
	if(!is.logical(boxplots))  stop("Error: 'boxplots' has to be a logical")
	if(!is.logical(out.object))stop("Error: 'out.object' has to be a logical")

	#------------------------------ Initialization -----------------------------#
	# Number of methods
	nb.meth <- length(methods)	

	# Merge predictors and response variable
	data <- as.data.frame(cbind(x,y))
	ncol <- ncol(data)
	nrow <- nrow(data)
	colnames(data)[ncol] <- "pred"
	y <- data[,ncol]

	# Index matrix
	names.ind <- c("Accuracy","Precision","Recall","F-score","Kappa","Specificity", "TSS")
	index <- matrix(NA, nrep, length(names.ind), dimnames=list(paste("rep",1:nrep,sep=""),
                   names.ind))
	if(sp){names.ind_sp <- c("sp_autocorrelation","Recall","Kappa","Specificity","TSS","AUC")
	index_sp <- matrix(NA, nrep, length(names.ind_sp), dimnames=list(paste("rep",1:nrep,sep=""),
                   names.ind_sp))}

	# Error matrix
	if(resampling=="boot"){
		nb.row.error <- nrep+3
		name.row.error <- c(paste("boot",1:nrep,sep=""),c("Mean",".632 error",".632+ error"))
	}
	else{
		nb.row.error <- nrep
		name.row.error <- paste("rep",1:nrep,sep="")
	}
	error <- matrix(NA, nb.row.error, 2,dimnames=list(name.row.error, c("Train","Test")))

	# Initialize probabilities, errors and indeces matrix of each method
	for(i in 1:nb.meth){
		n.prob <- paste("prob.", methods[i], sep="")
		n.err <- paste("error.", methods[i], sep="")
		n.ind <- paste("index.", methods[i], sep="")
		if(sp){n.ind_sp <- paste("index_sp.", methods[i], sep="")}
		if(resampling=="cv") n.err.cv <- paste("error.cv.", methods[i], sep="")
		
		assign(n.prob, matrix(NA, nrow, nrep))
		
		assign(n.err, error)
		assign(n.ind, index)
		if(sp){assign(n.ind_sp, index_sp)}
		if(resampling=="cv") assign(n.err.cv, NULL)
	}
	resample1 <- resample(data, method=resampling, train=train.frac, cv=k)
	# Output
	true.list <- test.id.list <- bal <- NULL
	out <- results <- err <- uncert <- summary.index <- summary.index_sp <- pp.all <- NULL
	name.error <- ifelse(resampling=="cv","CV error", 
                           ifelse(resampling=="boot",".632+ error","holdout error"))

	# Name for files
	type.model <- ifelse(is.null(bagging), 
                           ifelse(is.null(boosting),"","boosted_"),"bagged_")
	type.resampling <- paste(resampling,"_",sep="")
	name.root <- paste(paste("compar_", type.model,      sep=""),
                         paste("models_", type.resampling, sep=""),sep="")

	#----------------------------------- Body ----------------------------------#
	# Current CPU Time
	ptm1 <- proc.time()

	## Train and test simple, bagged or boosted models
	for(n in 1:nrep){
		cat(ifelse(resampling=="boot","b = ","rep = "),n,"\n")
		# Generate samples
		resample1 <- resample(data, method=resampling, train=train.frac, cv=k)
		data.train <- resample1$train
		if(resampling!="cv"){
			data.test <- resample1$test
			bal <- c(bal, sum(data.train[,ncol])/nrow(data.train))
			test.id.list <- c(test.id.list, list(resample1$test.ind))
			true.list <- c(true.list, list(data.test[,ncol]))
		}else{
			cv.groups <- resample1$groups
			test.id.list <- c(test.id.list, list(unlist(cv.groups)))
			true.list <- c(true.list, list(data.train[,ncol]))
		}

		# Train and test each model, fill associated objects
		for(i in 1:nb.meth){
			# Names objects
			n.res  <- paste(methods[i], ".results", sep="")
			n.prob <- paste("prob."   , methods[i], sep="")
			n.err  <- paste("error."  , methods[i], sep="")
			n.ind  <- paste("index."  , methods[i], sep="")
			if(sp){n.ind_sp  <- paste("index_sp."  , methods[i], sep="")}
			if(resampling=="cv") n.err.cv <- paste("error.cv.", methods[i], sep="")
			if(resampling=="cv") n.temp   <- paste("temp.", methods[i], sep="")

			# Initialisation of inficator matrix
			#prob.mat  <- matrix(F, nrow, nrep)
			prob.mat  <- matrix(F, nrow, nrep)
			error.mat <- matrix(F, nrep+3,2)
			index.mat <- matrix(F, nrep, 5)
			if(sp){index.mat_sp <- matrix(F, nrep, 5)}
			if(resampling=="cv") assign(n.temp,NULL)

			# Set indeces
			if(resampling!="cv") prob.indeces <- which(resample1$test.ind) + nrow*(n-1)
			if(resampling=="boot") error.indeces <- c(n, n+(nrep+3))
			if(resampling!="boot") error.indeces <- c(n, n+nrep)
			index.indeces <- c(n, (1:(length(names.ind)-1)*nrep)+n)
			if(sp){index.indeces_sp <- c(n, (1:(length(names.ind_sp)-1)*nrep)+n)}

			# Set indicator matrix
			if(resampling!="cv") prob.mat  <- replace(prob.mat,  prob.indeces, T)
			if(resampling=="cv") prob.mat[,n] <- T
			error.mat <- replace(error.mat, error.indeces, T)
			index.mat <- replace(index.mat, index.indeces, T)
			if(sp){index.mat_sp <- replace(index.mat_sp, index.indeces_sp, T)}

			# Store results from training and testing
			if(resampling!="cv") cv.groups <- NA
			for(j in 1:length(cv.groups)){
				if(resampling!="cv"){
					train <- data.train
					save(train,file=paste("train_rep",n,"_",bg,nb.var,"_",reduction,".rda",sep=""))
					test  <- data.test
					save(test,file=paste("test_rep",n,"_",bg,nb.var,"_",reduction,".rda",sep=""))

				}
				else{
					train_coord <- data.train[-cv.groups[[j]],]
					save(train_coord,file=paste("train_rep",n,"_cv.group",j,"_",bg,nb.var,"_",reduction,".rda",sep=""))
					train <- data.train[-cv.groups[[j]],-(1:2)]

					test_coord  <- data.train[cv.groups[[j]],]
					save(test_coord,file=paste("test_rep",n,"_cv.group",j,"_",bg,nb.var,"_",reduction,".rda",sep=""))
					test <- data.train[cv.groups[[j]],-(1:2)]

				}
				if(is.null(bagging) & is.null(boosting)){
					assign(n.res, train.and.test.model(
							     train=train,
							     method=methods[i],
							     test=test,
							     resp.var.name="pred",
							     knn.par=knn.par, svm.par=svm.par, 
                                               nnet.par=nnet.par, nbtree=nbtree,
							     Wts=Wts,
								 beta=beta,test_coord=test_coord,sp,coord.spind_order=coord.spind_order))
								 result<-get(n.res)
								 save(result,file=paste("res_compar_"  , methods[i],"_rep",n,"_cv.group",j,"_",bg,nb.var,"_",reduction,".rda",sep=""))
				}
				if(!is.null(bagging)){
					assign(n.res, train.and.test.bagged.model(
							     train=train, 
							     methods=methods[i],
							     test=test, 
							     B=bagging, resp.var.name="pred",
							     knn.par=knn.par, svm.par=svm.par, 
                                               nnet.par=nnet.par, nbtree=nbtree,
							     Wts=Wts,
								 beta=beta))
				}
				if(!is.null(boosting)){
					assign(n.res, train.and.test.boosted.model(
							     train=train, 
							     methods=methods[i],
						  	     test=test, 
							     M=boosting, resp.var.name="pred",
							     knn.par=knn.par, svm.par=svm.par, 
                                               nnet.par=nnet.par, nbtree=nbtree,
							     Wts=Wts,
								 beta=beta))
				}
				if(resampling=="cv"){
					assign(n.err.cv, rbind(get(n.err.cv), 
                                                     c(get(n.res)$error.train,
                                                       get(n.res)$error.test)))
					assign(n.temp, rbind(get(n.temp),  cbind(get(n.res)$prob.test, cv.groups[[j]])))
						
			}
			}

			# Update prob matrix of the associated method
			if(resampling!="cv") assign(n.prob,  replace(get(n.prob), 
                                              prob.mat, get(n.res)$prob.test))
			if(resampling=="cv"){
				proba.res.temp <- get(n.temp)[order(get(n.temp)[,2]),1]
				assign(n.prob,  replace(get(n.prob), prob.mat, proba.res.temp))
			}

			# Train & test error as cv average if cv
			if(resampling!="cv") trainNtest <- c(get(n.res)$error.train,get(n.res)$error.test)
			if(resampling=="cv") trainNtest <- apply(get(n.err.cv),2,mean)

			# Update error matrix of the associated method
			assign(n.err, replace(get(n.err), error.indeces, trainNtest))

			# Update index matrix of the associated method
			assign(n.ind, replace(get(n.ind), index.indeces,get(n.res)$indexes.test))
			if(sp){assign(n.ind_sp, replace(get(n.ind_sp), index.indeces_sp,get(n.res)$indexes.test.spcorrected))}
		}
	}

cat("\n\n")
	## Arrange results of each method
	for(i in 1:nb.meth){
		# Create name of object for access
		n.err  <- paste("error."  , methods[i], sep="")
		n.ind  <- paste("index."  , methods[i], sep="")
		if(sp){n.ind_sp <- paste("index_sp."  , methods[i], sep="")}
		n.prob <- paste("prob."   , methods[i], sep="")

		# Calculate the .632+ estimate of error for each model
		if(resampling=="boot"){
			cat("Calculation .632+ error \n")
			assign(n.err, errors.632plus.model(data=data, method=methods[i], 
                                    y=y, err.mat=get(n.err), B=nrep, resp.var.name="pred",
                                    knn.par=knn.par, svm.par=svm.par, nnet.par=nnet.par, 
                                    nbtree=nbtree, bagging=bagging, boosting=boosting))
		}
		else{
			assign(n.err, rbind(get(n.err), "Mean"=apply(get(n.err),2,mean)))
		}

		# Change structure of index for output
		assign(n.ind, summary.indexes(get(n.ind), sp=FALSE,  boxplots=T))
		if(sp){assign(n.ind_sp, summary.indexes(get(n.ind_sp), sp=TRUE,  boxplots=T))}

		# Change structure of prob for output
		assign(n.prob, summary.post.prob(get(n.prob),  nrep, y))

		# Objects for output
		if(sp){results <- c(results, list(c(get(n.prob), "errors"=list(get(n.err)), get(n.ind), get(n.ind_sp))))}
		else{results <- c(results, list(c(get(n.prob), "errors"=list(get(n.err)), get(n.ind))))} 
		temp    <- ifelse(resampling!="cv",get(n.err)[(nrep+3)],get(n.err)[nrep+1,2])
		err     <- rbind(err, temp)
		uncert  <- rbind(uncert, results[[i]]$uncertainty)
		summary.index <- rbind(summary.index, get(n.ind)$summary[1,])
		if(sp){summary.index_sp <- rbind(summary.index_sp, get(n.ind_sp)$summary[1,])}

		if(plots) pp.all <- c(pp.all,list(cbind(get(n.prob)$post.prob[,1], y)))

		# Write results for predictions obtained by averaging posterior probabilities
		if(file){
			name.file <- paste(paste(paste(species.name,"_",name.root,"prob_",sep=""),
                                                 methods[i],sep=""),"_",bg,nb.var,"_",reduction,".txt",sep="")
			write.table(get(n.prob)$post.prob, file=name.file,  row.names=F)
		}
	}
	# Add names to output objects
	names(results) <- rownames(err) <- rownames(uncert) <- rownames(summary.index) <- methods
	if(sp){rownames(summary.index_sp) <- methods}
	if(plots) names(pp.all) <- methods
	colnames(err) <- name.error
	colnames(uncert) <- "Uncertainty"

	# Current CPU Time
	ptm2 <- proc.time()
	ptm <- ptm2 - ptm1

cat("done\n\n")
	#---------------------------------- Output ---------------------------------#
	## All results for each classifier, i.e:
	# probabilities and sd, confusion matrix, indexes, errors
	out$results   <- results

	## General results
	# .632+ error
	if(resampling=="boot")    out$err.632plus <- err
	if(resampling=="cv")      out$err.cv      <- err
	if(resampling=="holdout") out$err.holdout <- err

	# Uncertainty
	out$uncertainty <- uncert

	# Table of ranks and score for decision
	ranks <- cbind(err, rank(err), uncert, rank(uncert), (rank(err)+rank(uncert)))
	colnames(ranks) <- c(name.error, "Rank1", "Uncertainty", "Rank2", "Score")
	out$ranking.table <- ranks

	# Averaged indexes only
	out$indexes <- summary.index
	if(sp){out$indexes_sp <- summary.index_sp}

	# Trace
	out$trace.rep <- test.id.list

	# List of response variable vectors of non-bootstrapped observations
	out$true.list <- true.list

	## If results must be save in a file
	if(file){
		# General results
		write.table(out$err,          file=paste(species.name,"_",name.root,"general_error_",bg,nb.var,"_",reduction,".txt",sep=""))
		write.table(out$uncertainty,  file=paste(species.name,"_",name.root,"general_uncertainty_",bg,nb.var,"_",reduction,".txt",sep=""))
		write.table(out$indexes,      file=paste(species.name,"_",name.root,"general_indexes_",bg,nb.var,"_",reduction,".txt",sep=""))
		if(sp){write.table(out$indexes_sp,      file=paste(species.name,"_",name.root,"general_indexes_sp_",bg,nb.var,"_",reduction,".txt",sep=""))}
		write.table(out$ranking.table,file=paste(species.name,"_",name.root,"general_ranks_table_",bg,nb.var,"_",reduction,".txt",sep=""))
	}

	## Plots
	if(plots.all){
		# ROC plots from all bootstrapped
		for(i in 1:nb.meth){
			n.prob <- paste("prob.",methods[i],sep="")
			plot.roc.all(prob=get(n.prob)$prob.list, true=true.list, 
                              name.file=paste(name.root, methods[i], sep=""))
		}
	}
	if(plots){
		# ROC plots from average posterior probabilities
		perfs <- plot.perf(pp = pp.all,
                               names = methods,
                               colors = plots.col,
                               models.plots = 1:nb.meth,
                               file = name.root)
		out$aucROC <- perfs$aucROC
		out$aucLift <- perfs$aucLift
	}
	if(boxplots){
		X11(width=15,height=10)
		boxplot(summary.index)
		savePlot(paste(name.root,"boxplot_indexes_",bg,nb.var,"_",reduction,".png",sep=""),type="png")
		# Boxplots of indexes
		#index <- NULL
		#for(i in 1:nb.meth){
		#	n.ind <- paste("index.",methods[i],sep="")
		#	assign("index", c(index, list(get(n.ind))))
		#}
		#boxplots.indexes(index = index,
          #                   names = methods,
           #                  colors= box.col,
            #                 file=name.root)

		if(sp){
		X11(width=15,height=10)
		boxplot(summary.index_sp)
		savePlot(paste(name.root,"boxplot_indexes_sp_",bg,nb.var,"_",reduction,".png",sep=""),type="png")}
		# Boxplots of indexes
		#index_sp <- NULL
		#for(i in 1:nb.meth){
		#	n.ind_sp <- paste("index_sp.",methods[i],sep="")
		#	assign("index_sp", c(index_sp, list(get(n.ind_sp))))
		#}
		#boxplots.indexes(index = index_sp,
           #                  names = methods,
            #                 colors= box.col,
			#		    sp=TRUE,
	}

	if(out.object) save(out, file = paste(species.name,"_",name.root,"output_object_",bg,nb.var,"_",reduction,".Rdata",sep=""))
	return(out)
}

