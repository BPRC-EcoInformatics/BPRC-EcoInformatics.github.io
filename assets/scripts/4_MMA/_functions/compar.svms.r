# Comparison of supervised svms using different performance metrics

# Author  :   Tak Ikeda
# Updated :   January 2010



compar.svms <- function(x,y,
                          nrep=200, 
                          methods="svm",
			  svm.par.mat,
                          bagging  = NULL,
                          boosting = NULL,
                          knn.par  = NULL,
                          svm.par  = NULL, 
                          nnet.par = NULL,
                          nbtree   = 10000,
			  weights = NULL,
			  Wts  = NULL,
                          resampling="cv",
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
                          out.object=T){

	#------------------------------ Error handling -----------------------------#

	# Check number of SVMs to compare is greater than 1
	if(nrow(svm.par.mat)<=1) stop("Error: Must compare 2 or more SVMs!")
	

	# Check x and y
	if(is.null(dim(x))) stop("Error: dimensions of data frame 'x'")
	if(!is.null(dim(y))) stop("Error: dimensions of vector 'y'")
	if(dim(x)[1]!=length(y)) stop("Error: compatibility between 'x' and 'y'")

	# Check nrep
	if(nrep<2) stop("A minimum of two repetitions is required")



	#------------------------------ Initialization -----------------------------#
	# Number of methods
	nb.meth <- 1	
	
	# Number of SVMs
	nb.svm <- nrow(svm.par.mat)


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

	# Error matrix
	nb.row.error <- nrep
	name.row.error <- paste("rep",1:nrep,sep="")
	error <- matrix(NA, nb.row.error, 2,
                      dimnames=list(name.row.error, c("Train","Test")))

	# Initialize probabilities, errors and indeces matrix of each method
	svm.names <- paste(svm.par.mat[,1],svm.par.mat[,2],sep="_")
	for(i in 1:nb.svm){
		n.prob <- paste("prob.", svm.names[i], sep="")
		n.err <- paste("error.", svm.names[i], sep="")
		n.ind <- paste("index.", svm.names[i], sep="")
		n.err.cv <- paste("error.cv.", svm.names[i], sep="")
		assign(n.prob, matrix(NA, nrow, nrep))
		assign(n.err, error)
		assign(n.ind, index)
		assign(n.err.cv, NULL)
	}

	# Output
	true.list <- test.id.list <- bal <- NULL
	out <- results <- err <- uncert <- summary.index <- pp.all <- NULL
	name.error <- "CV error"
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
#nreps <- 1
	for(n in 1:nrep){
		cat(ifelse(resampling=="boot","b = ","rep = "),n,"\n")
		# Generate samples
		resampl <- resample(data, method=resampling, train=train.frac, cv=k)
		data.train <- resampl$train
		cv.groups <- resampl$groups
		test.id.list <- c(test.id.list, list(unlist(cv.groups)))
		true.list <- c(true.list, list(data.train[,ncol]))
		

		# Train and test each model, fill associated objects
#i<-1
		for(i in 1:nb.svm){
			# Names objects
			n.res  <- paste(svm.names[i], ".results", sep="")
			n.prob <- paste("prob."   , svm.names[i], sep="")
			n.err  <- paste("error."  , svm.names[i], sep="")
			n.ind  <- paste("index."  , svm.names[i], sep="")
			n.err.cv <- paste("error.cv.", svm.names[i], sep="")
			n.temp   <- paste("temp.", svm.names[i], sep="")

			# Initialisation of inficator matrix
			prob.mat  <- matrix(F, nrow, nrep)
			error.mat <- matrix(F, nrep+3,2)
			index.mat <- matrix(F, nrep, 5)
			assign(n.temp,NULL)

			# Set indeces
			error.indeces <- c(n, n+nrep)
			index.indeces <- c(n, (1:(length(names.ind)-1)*nrep)+n)

			# Set indicator matrix
			prob.mat[,n] <- T
			error.mat <- replace(error.mat, error.indeces, T)
			index.mat <- replace(index.mat, index.indeces, T)

			# Store results from training and testing
#j<-1
			for(j in 1:length(cv.groups)){
				train <- data.train[-cv.groups[[j]],]
				test  <- data.train[cv.groups[[j]],]
				if(is.null(bagging) & is.null(boosting)){
					assign(n.res, train.and.test.model(
							     train=train,
							     method="svm",
							     test=test,
							     resp.var.name="pred",
							     knn.par=knn.par, 
						svm.par=c("rbfdot",svm.par.mat$c1[j],svm.par.mat$c2[j]), 
                                               nnet.par=nnet.par, nbtree=nbtree,
							     Wts=Wts))
				}
				if(resampling=="cv"){
					assign(n.err.cv, rbind(get(n.err.cv), 
                                                     c(get(n.res)$error.train,
                                                       get(n.res)$error.test)))
					assign(n.temp, rbind(get(n.temp), 
                                                   cbind(get(n.res)$prob.test, cv.groups[[j]])))
				}
			}

			# Update prob matrix of the associated method
			assign(n.prob,  replace(get(n.prob), prob.mat, 
                                               get(n.temp)[order(get(n.temp)[,2]),1]))

			# Train & test error as cv average if cv
			trainNtest <- apply(get(n.err.cv),2,mean)

			# Update error matrix of the associated method
			assign(n.err, replace(get(n.err), error.indeces, trainNtest))

			# Update index matrix of the associated method
			assign(n.ind, replace(get(n.ind), index.indeces,get(n.res)$indexes))
		}
	}	

cat("\n\n")
	## Arrange results of each method
	for(i in 1:nb.svm){
		# Create name of object for access
		n.err  <- paste("error."  , svm.names[i], sep="")
		n.ind  <- paste("index."  , svm.names[i], sep="")
		n.prob <- paste("prob."   , svm.names[i], sep="")

		# Calculate the .632+ estimate of error for each model
		assign(n.err, rbind(get(n.err), "Mean"=apply(get(n.err),2,mean)))

		# Change structure of index for output
		assign(n.ind, summary.indexes(get(n.ind),  boxplots=T))

		# Change structure of prob for output
		assign(n.prob, summary.post.prob(get(n.prob),  nrep, y))

		# Objects for output
		results <- c(results, list(c(get(n.prob), "errors"=list(get(n.err)), get(n.ind))))
		temp    <- ifelse(resampling!="cv",get(n.err)[(nrep+3)],get(n.err)[nrep+1,2])
		err     <- rbind(err, temp)
		uncert  <- rbind(uncert, results[[i]]$uncertainty)
		summary.index <- rbind(summary.index, get(n.ind)$summary[1,])

		if(plots) pp.all <- c(pp.all,list(cbind(get(n.prob)$post.prob[,1], y)))

		# Write results for predictions obtained by averaging posterior probabilities
		if(file){
			name.file <- paste(paste(paste(name.root,"prob_",sep=""),
                                                 svm.names[i],sep=""),".txt",sep="")
			write.table(get(n.prob)$post.prob, file=name.file,  row.names=F)
		}
	}
	# Add names to output objects
	names(results) <- rownames(err) <- rownames(uncert) <- rownames(summary.index) <- svm.names
	if(plots) names(pp.all) <- svm.names
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
 	out$err.cv      <- err


	# Uncertainty
	out$uncertainty <- uncert

	# Table of ranks and score for decision
	ranks <- cbind(err, rank(err), uncert, rank(uncert), (rank(err)+rank(uncert)))
	colnames(ranks) <- c(name.error, "Rank1", "Uncertainty", "Rank2", "Score")
	out$ranking.table <- ranks

	# Averaged indexes only
	out$indexes <- summary.index

	# Trace
	out$trace.rep <- test.id.list

	# List of response variable vectors of non-bootstrapped observations
	out$true.list <- true.list


	# Time
	out$time <- ptm
	cat("Time:\n")
	print(ptm)

	## If results must be save in a file
	if(file){
		# General results
		write.table(out$err,          file=paste(name.root,"general_error.txt",sep=""))
		write.table(out$uncertainty,  file=paste(name.root,"general_uncertainty.txt",sep=""))
		write.table(out$indexes,      file=paste(name.root,"general_indexes.txt",sep=""))
		write.table(out$ranking.table,file=paste(name.root,"general_ranks_table.txt",sep=""))
	}

	## Plots
	if(plots.all){
		# ROC plots from all bootstrapped
		for(i in 1:nb.svm){
			n.prob <- paste("prob.",svm.names[i],sep="")
			plot.roc.all(prob=get(n.prob)$prob.list, true=true.list, 
                              name.file=paste(name.root, svm.names[i], sep=""))
		}
	}
	if(plots){
		# ROC plots from average posterior probabilities
		perfs <- plot.perf(pp = pp.all,
                               names = svm.names,
                               colors = plots.col,
                               models.plots = 1:nb.svm,
                               file = name.root)
		out$aucROC <- perfs$aucROC
		out$aucLift <- perfs$aucLift
	}
	if(boxplots){
		windows()
		# Boxplots of indexes
		index <- NULL
		for(i in 1:nb.svm){
			n.ind <- paste("index.",svm.names[i],sep="")
			assign("index", c(index, list(get(n.ind))))
		}
		boxplots.indexes(index = index,
                             names = svm.names,
                             colors= box.col,
                             file=name.root)
	}

	if(out.object) save(out, file = paste(name.root,"output_object.Rdata",sep=""))
	return(out)
}

