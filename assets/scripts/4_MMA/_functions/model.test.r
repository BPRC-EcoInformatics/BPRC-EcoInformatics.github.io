# Test of a supervised model on an independent test set

# Note: this function is made for handling high dimensional 
# test set such as the worldclim database. It splits the vector
# of posterior probabilities and write it in small files. It
# allows the function 'read.model.test' to read them quite fastly.


# Author  :   Gwénaël Leday
# Updated :   June 2009



model.test <- function(x,
                       y,
                       test, 
                       method,
                       nrep,
                       resampling="boot",
                       train.frac=0.75,
                       k=10,
                       bagging  = NULL,
                       boosting = NULL,
                       knn.par  = NULL,
                       svm.par  = NULL,
                       nnet.par = NULL,
                       nbtree   = 10000,
                       n.split = 20000,
                       file=F,
                       out.object=F){

	## Error handling
	# Check x and y
	if(is.null(dim(x))) stop("Error: dimensions of data frame 'x'")
	if(!is.null(dim(y))) stop("Error: dimensions of vector 'y'")
	if(dim(x)[1]!=length(y)) stop("Error: compatibility between 'x' and 'y'")

	# Check method
	bool <- methods %in% c("LDA","QDA","LOG","NB","CART","STUMP",
                                    "CTREE","KNN","SVM","NNET","RF")

	if(sum(bool)!=length(method)) stop("Error: do not recognize name of method")

	# Check KNN parameter
	if(is.null(knn.par) & "KNN"%in%method){
		stop("Error: Which value for the parameter knn.par ?")
	}

	# Check SVM parameters
	if(is.null(svm.par) & "SVM"%in%method) stop("Error: parameters of svm")
	if(!is.null(svm.par)){
		if(!is.vector(svm.par)) stop("Error: parameters of svm")
		if(!(svm.par[1] %in% c("vanilladot","rbfdot","polydot",
                             "linear","radial","polynomial"))){
			stop("Error: input kernel")
		}
		if(svm.par[1]=="linear") svm.par[1] <- "vanilladot"
		if(svm.par[1]=="radial") svm.par[1] <- "rbfdot"
		if(svm.par[1]=="polynomial") svm.par[1] <- "polydot"
		if(svm.par[1]=="vanilladot" && length(svm.par)>2){
			stop("Error: too much parameters for the linear kernel")
		}
		if(svm.par[1]=="rbfdot" && length(svm.par)>3){
			stop("Error: too much parameters for the radial basis kernel")
		}
		if(svm.par[1]=="polydot" && length(svm.par)>5){
			stop("Error: too much parameters for the polynomial kernel")
		}
	}

	# Check NNET parameters
	if(is.null(nnet.par) & "NNET"%in%method) stop("Error: parameters of nnet")
	if(!is.null(nnet.par) & !is.vector(nnet.par)) stop("Error: three values are required for nnet.par")

	# Check boosting
	bool <- c("KNN","SVM","RF") %in% method
	if(!is.null(boosting) & sum(bool)!=0){
		stop("Error: boosting does not support input method(s)")
	}

	# Check n.split
	if(is.null(n.split)) stop("Error: no value for 'n.split'")
	n.split <- trunc(n.split)

	# Check nbtree, k and train.frac
	if(length(nbtree)!=1 | !is.numeric(nbtree)) stop("Error: 'nbtree' value")
	nbtree <- trunc(nbtree)

	## Initialization
	# Trainning set
	data <- as.data.frame(cbind(x,y))
	ncol <- ncol(data)
	colnames(data)[ncol] <- "pred"

	## Splits
	split <- c(1,1+n.split*1:floor(nrow(test)/n.split), nrow(test)+1)
	split.rep <- length(split)-1

	## Body
	for(n in 1:nrep){
		cat(ifelse(resampling=="boot","b = ","rep = "),n,"\n")
		# Generate samples
		resample <- resample(data, method=resampling, train=train.frac, cv=k)
		data.train <- resample$train
		if(resampling!="cv"){
			cv.groups <- NA
		}else{
			cv.groups <- resample$groups
		}

		for(j in 1:length(cv.groups)){
			if(resampling!="cv"){
				train <- data.train
			}
			else{
				train <- data.train[-cv.groups[[j]],]
			}

			if(is.null(bagging) & is.null(boosting)){
				res <- train.and.test.model(
						train=train,
						method=method,
						test=test,
						resp.var.name="pred",
						knn.par=knn.par, svm.par=svm.par, 
     		                   nnet.par=nnet.par, nbtree=nbtree)
			}
			if(!is.null(bagging)){
				res <- train.and.test.bagged.model(
						train=train, 
						method=method,
						test=test, 
						B=bagging, resp.var.name="pred",
						knn.par=knn.par, svm.par=svm.par, 
     		                   nnet.par=nnet.par, nbtree=nbtree)
			}
			if(!is.null(boosting)){
				res <- train.and.test.boosted.model(
						train=train, 
						methods=method,
					  	test=test, 
						M=boosting, resp.var.name="pred",
						knn.par=knn.par, svm.par=svm.par, 
						nnet.par=nnet.par, nbtree=nbtree)
			}
			if(resampling=="cv"){
				vect <- NULL
				vect <- cbind(vect, res$prob.test)
			}
		}

		# Write vector of predictions
		if(resampling=="cv") vect <- apply(vect, 1, mean)
		if(resampling!="cv") vect <- res$pred.test
		for(i in 1:split.rep){
			name.file <- paste(paste(c("pred",resampling,n,i),sep="",collapse="_"),".txt",sep="")
	   		write.table(vect[split[i]:(split[i+1]-1)],
                              name.file,append=F,row.names=F, col.names=F)
		}
	}
cat("\n\n")
cat("Done \n")
}

