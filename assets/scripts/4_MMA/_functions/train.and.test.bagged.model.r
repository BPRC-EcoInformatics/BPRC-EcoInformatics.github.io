# Bagging on models

# Author  :   Gwénaël Leday
# Updated :   May 2009



train.and.test.bagged.model <- function(train, 
                                        methods=NULL,
                                        test=NULL, 
                                        B=200, 
                                        resp.var.name, 
                                        knn.par=NULL, 
                                        svm.par=NULL, 
                                        nnet.par=NULL,
                                        nbtree=NULL,
                                        file=F,
                                        file.name=NULL){

	## Error handling
	# Check B
	if(B<2) stop("Error: A minimum of two booststrap samples is required")

	# Check train and resp.var.name
	if(is.null(train)) stop("Error: no input data")
	if(is.null(resp.var.name)) stop("Error: no input name of response variable")
	if(!is.character(resp.var.name)) stop("Error: variable name is not a character!")

	# Check methods
	bool <- methods %in% c("LDA","QDA","LOG","NB","CART", "STUMP", 
                        "CTREE","KNN","SVM","NNET","RF")
	if(sum(bool)!=length(methods)) stop("Error: do not recognize name of method")

	# Check KNN parameter
	if(is.null(knn.par) & "KNN"%in%methods){
		stop("Error: Which value for the parameter knn.par ?")}

	# Check SVM parameters
	if(is.null(svm.par) & "SVM"%in%methods) stop("Error: parameters of svm")
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

	# Check NNET parameters
	if(is.null(nnet.par) & "NNET"%in%methods) stop("Error: parameters of nnet")
	if(!is.vector(nnet.par)) stop("Error: three values are required for nnet.par")

	# Check nbtree
	if(length(nbtree)!=1 | !is.numeric(nbtree)) stop("Error: 'nbtree' value")

	## Initialization
	# Number of methods
	nb.meth <- length(methods)

	# Merge predictors and response variable for training set
	data.train <- as.data.frame(train)
	nrow1 <- nrow(data.train)
	ncol1 <- ncol(data.train)
	colnames(data.train)[ncol1] <- "pred"

	# Merge predictors and response variable for test set
	if(!is.null(test)){
		data.test <- as.data.frame(test)
		nrow2 <- nrow(data.train)
		ncol2 <- ncol(data.test)
		colnames(data.test)[ncol2] <- "pred"
	}
	else{data.test <- NULL}

	# Predictions of each model
	for(i in 1:nb.meth){
		n.pred1 <- paste(methods[i],".pred.train", sep="")
		n.pred2 <- paste(methods[i],".pred.test", sep="")
		assign(n.pred1, NULL)
		assign(n.pred2, NULL)
		assign(methods[i], NULL)
	}

	## Body
	cat("Bagging ",paste(methods, sep="", collapse=", "),"\t b =  ")
	# Learn models on each bootstrap sample and test it on
	# the input training and test sets
	for(b in 1:B){
	   if(b!=1){if(b>9){if(b>99){if(b>999){cat("\b\b\b\b")}else{
         cat("\b\b\b")}}else{cat("\b\b")}}else{cat("\b")}}
	   cat(b) 
		boot <- sample(nrow1, nrow1, replace=TRUE)
		data.train.boot <- data.train[boot,]

		# Train and test each model, fill associated objects
		for(i in 1:nb.meth){

			# Names objects
			n.res <- paste(methods[i], ".results", sep="")
			n.pred1 <- paste(methods[i],".pred.train", sep="")
			n.pred2 <- paste(methods[i],".pred.test", sep="")

			# Store results from training and testing
			assign(n.res, train.and.test.model(train=data.train.boot,
                                                     method=methods[i],
                                                     test=data.test,
                                                     test2=data.train,
                                                     resp.var.name="pred",
                                                     knn.par=knn.par,
                                                     svm.par=svm.par,
                                                     nnet.par=nnet.par,
                                                     nbtree=nbtree))
			if(!is.null(test)) assign(n.pred2, cbind(get(n.pred2), get(n.res)$pred.test))
			assign(n.pred1, cbind(get(n.pred1), get(n.res)$pred.test2))
		}
	}
	cat("\n")

	## Vote of models
	# Decision function
	f <- function(x){ifelse(sum(x)>=ceiling(B/2),1,0)}

	for(i in 1:nb.meth){

		# Names objects
		n.pred1 <- paste(methods[i],".pred.train", sep="")
		n.pred2 <- paste(methods[i],".pred.test", sep="")
		n.mod.pred1 <- paste(paste("models.",methods[i],sep=""),".pred.train", sep="")
		n.mod.pred2 <- paste(paste("models.",methods[i],sep=""),".pred.test", sep="")
		n.mod.prob2 <- paste(paste("models.",methods[i],sep=""),".prob.test", sep="")
		n.err1 <- paste("error.train.", methods[i], sep="")
		n.err2 <- paste("error.test.", methods[i], sep="")
		n.conf <- paste(paste("conf.",methods[i],sep=""),".test", sep="")
		n.index <- paste(paste("index.",methods[i],sep=""),".test", sep="")

		# Store information on training
		assign(n.mod.pred1, apply(get(n.pred1),1,f))
		assign(n.err1, 1- compute.indexes(confusion.matrix(get(n.mod.pred1),
                                                    as.integer(data.train$pred)))[1,1])

		# Prepare output object
		assign(methods[i], c(get(methods[i]),list("pred.train"=get(n.mod.pred1), 
                                                      "error.train"=get(n.err1))))

		# Write in a file if required
		if(file) write.table(get(n.err1), paste(file.name,"train_err.txt",sep=""), 
                                              row.names=F, col.names=F)

		# Store information on testing (if testing)
		if(!is.null(test)){ 
			# Store
			assign(n.mod.prob2, apply(get(n.pred2),1,sum)/B)
			assign(n.mod.pred2, apply(get(n.pred2),1,f))
			assign(n.conf, confusion.matrix(get(n.mod.pred2),as.integer(data.test$pred)))
			assign(n.index, compute.indexes(get(n.conf)))
			assign(n.err2, 1 - get(n.index)[1,1])

			# Write in a file if required
			if(file){
				write.table(get(n.mod.prob2), paste(file.name,"test_prob.txt",sep=""), 
                                                            row.names=F, col.names=F)
				write.table(get(n.err2), paste(file.name,"test_err.txt",sep=""), 
                                                            row.names=F, col.names=F)
				write.table(get(n.index), paste(file.name,"test_indexes.txt",sep=""), 
                                                            row.names=F, col.names=F)
			}

			# Prepare output object
			assign(methods[i], c(get(methods[i]),list("pred.test"=get(n.mod.pred2),
                                                            "prob.test"=get(n.mod.prob2),
                                                            "conf.test"=get(n.conf),
                                                            "indexes.test"=get(n.index),
                                                            "error.test"=get(n.err2))))
		}
	}

	## Output
	out <- NULL
	for(i in 1:nb.meth){
		assign("out", c(out, get(methods[i])) )
	}
	if(nb.meth>1) names(out) <- methods
	return(out)
}

