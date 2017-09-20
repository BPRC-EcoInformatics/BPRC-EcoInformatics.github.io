# Train and test a specified model

# Author  :   Gwena?E Leday
# Updated :   May 2009
# Edited: July 2012, by TI, to handle different betas for f-score, default is 1


train.and.test.model_test <- function(train,
                                 method=NULL,
                                 test=NULL,
                                 test2=NULL,
                                 test3=NULL,
                                 resp.var.name, 
                                 knn.par=NULL, 
                                 svm.par=NULL, 
                                 nnet.par=NULL,
                                 nbtree=NULL,
                                 Wts=NULL,
                                 priors=NULL,
                                 prob=T,
								 beta=1){

	## Error handling
	# Check train and resp.var.name
	if(is.null(train)) stop("Error: no input data")
	if(is.null(resp.var.name)) stop("Error: no input name of response variable")
	if(!is.character(resp.var.name)) stop("Error: variable name is not a character!")

	# Check method
	bool <- method %in% c("LDA","QDA","LOG","NB","CART", "STUMP", 
                        "CTREE","KNN","SVM","NNET","RF")
	if(!bool) stop("Error: wrong method")

	# Check KNN parameter
	if(is.null(knn.par) & "KNN"%in%method){
		stop("Error: Which value for the parameter knn.par ?")}

	# Check SVM parameters
	if(is.null(svm.par) & "SVM"%in%method) stop("Error: parameters of svm")
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
	if(is.null(nnet.par) & "NNET"%in%method) stop("Error: parameters of nnet")
	if(!is.null(nnet.par)){
		if(!is.vector(nnet.par)) stop("Error: three values are required for nnet.par")
	}

	# Check nbtree
	if(!is.null(nbtree)&(length(nbtree)!=1 | !is.numeric(nbtree))) stop("Error: 'nbtree' value")
	## Learn chosen model
	model <- learn.model_test(data=train,
                           resp.var.name = resp.var.name,
                           method = method,
                           svm.par = svm.par,
                           nnet.par = nnet.par,
                           nbtree = nbtree,
                           Wts = Wts,
                           priors = priors)
	model.pred <- predict.model_test(model = model,
                                  method = method,
                                  train = train,
                                  test = train, 
                                  resp.var.name = resp.var.name,
                                  knn.par = knn.par,
                                  svm.par = svm.par,
                                  nnet.par = nnet.par,
	                            Wts = Wts,
                                  prob = prob)$pred
	conf <- confusion.matrix(model.pred,train[,resp.var.name])
	error.train <- (conf[1,2] + conf[2,1])/sum(conf)

	## Test on test set: posterior probabilities, 
	## predictions, confusion matrix, error and indexes
	if(!is.null(test)){
		predictions <- predict.model_test(model = model,
                                         method = method,
                                         train = train,
                                         test = test,
                                         resp.var.name = resp.var.name,
                                         knn.par = knn.par,
                                         svm.par = svm.par,
                                         nnet.par = nnet.par,
                                         prob = prob)
		if(prob) prob.test <- predictions$prob
		pred.test <- predictions$pred
		if(ncol(train)==ncol(test)){
			conf.test <- confusion.matrix(pred.test,test[,resp.var.name])
			indexes.test <- compute.indexes(conf.test,beta=beta)
			error.test <- 1 - indexes.test[1,1]
		}
	}

	# Test on test set 2 (Optional)
#	if(!is.null(test2)){
#		predictions <- predict.model_test(model = model,
#                                         method = method,
#                                         train = train,
#                                         test = test2,
#                                         resp.var.name = resp.var.name,
#                                         knn.par = knn.par,
#                                         svm.par = svm.par,
#                                         nnet.par = nnet.par,
#                                         prob = prob)
#		if(prob) prob.test2 <- predictions$prob
##		if(ncol(train)==ncol(test2)){
#			conf.test2 <- confusion.matrix(pred.test2,test2[,resp.var.name])
#			indexes.test2 <- compute.indexes(conf.test2,beta=beta)
#			error.test2 <- 1 - indexes.test2[1,1]
#		}
#	}
	# Test on test set 3 (Optional)
#	if(!is.null(test3)){
#		predictions <- predict.model_test(model = model,
#                                         method = method,
#                                         train = train,
#                                         test = test3,
#                                         resp.var.name = resp.var.name,
#                                         knn.par = knn.par,
#                                         svm.par = svm.par,
#                                         nnet.par = nnet.par,
#                                         prob = prob)
#		if(prob) prob.test3 <- predictions$prob
#		pred.test3 <- predictions$pred
#		if(ncol(train)==ncol(test3)){
#			conf.test3 <- confusion.matrix(pred.test3,test3[,resp.var.name])
#			indexes.test3 <- compute.indexes(conf.test3,beta=beta)
#			error.test3 <- 1 - indexes.test3[1,1]
#		}
#	}

	## Output
	output <- NULL
	output$pred.train <- model.pred
	output$error.train <- error.train
	if(!is.null(test)){
		output$pred.test <- pred.test
		if(prob) output$prob.test <- prob.test
		if(ncol(train)==ncol(test)){
			output$conf.test <- conf.test
			output$indexes.test <- indexes.test
			output$error.test <- error.test
		}
	}
	if(!is.null(test2)){
		output$pred.test2 <- pred.test2
		if(prob) output$prob.test2 <- prob.test2
		if(ncol(train)==ncol(test2)){
			output$conf.test2 <- conf.test2
			output$indexes.test2 <- indexes.test2
			output$error.test2 <- error.test2
		}
	}
	if(!is.null(test3)){
		output$pred.test3 <- pred.test3
		if(prob) output$prob.test3 <- prob.test3
		if(ncol(train)==ncol(test3)){
			output$conf.test3 <- conf.test3
			output$indexes.test3 <- indexes.test3
			output$error.test3 <- error.test3
		}
	}

	return(output)
}