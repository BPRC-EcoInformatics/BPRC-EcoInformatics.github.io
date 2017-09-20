# Train and test boosted model(s)

# Note: This Adaboost.M1 algorithm

# Author  :   Gwénaël Leday
# Updated :   May 2009



train.and.test.boosted.model <- function(train, 
                                         methods=NULL,
                                         test=NULL, 
                                         M=200, 
                                         resp.var.name, 
                                         knn.par=NULL, 
                                         svm.par=NULL, 
                                         nnet.par=NULL,
                                         nbtree=NULL,
                                         file=F,
                                         file.name=NULL){

	## Error handling
	# Check M
	if(M<2) stop("Error: M<2")

	# Check train and resp.var.name
	if(is.null(train)) stop("Error: no input data")
	if(is.null(resp.var.name)) stop("Error: no input name of response variable")
	if(!is.character(resp.var.name)) stop("Error: variable name is not a character!")

	# Check methods
	bool <- methods %in% c("LDA","QDA","LOG","NB","CART", "STUMP", 
                        "CTREE","NNET")
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



	## Initialization
	# Number of methods
	nb.meth <- length(methods)

	# Merge predictors and response variable for training set
	data.train <- as.data.frame(train)
	nrow1 <- nrow(data.train)
	ncol1 <- ncol(data.train)
	colnames(data.train)[ncol1] <- "pred"
	y <- data.train[,ncol1]

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
		# Names of objects to create
		n.res   <- paste("res."       , methods[i], sep="")
		n.pred1 <- paste("pred.train.", methods[i], sep="")
		n.pred2 <- paste("pred.test." , methods[i], sep="")
		n.alpha <- paste("alpha."     , methods[i], sep="")
		n.w     <- paste("w."         , methods[i], sep="")
		n.err   <- paste("error."     , methods[i], sep="")
		n.stop  <- paste("stop."      , methods[i], sep="")
		n.m     <- paste("m."         , methods[i], sep="")
		n.prior <- paste("prior."     , methods[i], sep="")

		# Initialization
		assign(n.res, NULL)
		assign(n.pred1, NULL)
		assign(n.pred2, NULL)
		assign(n.alpha, NULL)
		assign(n.prior, NULL)
		assign(n.w, matrix(rep(1/nrow(data.train),nrow(data.train))))
		assign(n.err, NULL)
		assign(methods[i], NULL)
		assign(n.stop, FALSE)
		assign(n.m, 0)
	}

	# Counter
	cpt <- 0

	#----------------------------------- Body ----------------------------------#
	cat("Boosting ",paste(methods, sep="", collapse=", "),"\t m =\t ")
	# Adaboost algorithm that minimizes exponential loss function
	for(m in 1:M){
		# Print progress
		if(m!=1){if(m>9){if(m>99){if(m>999){cat("\b\b\b\b")}else{
		cat("\b\b\b")}}else{cat("\b\b")}}else{cat("\b")}}
		cat(m)

		# Train and test each model, fill associated objects
		for(i in 1:nb.meth){
			# Names of objects
			n.res   <- paste("res."       , methods[i], sep="")
			n.pred1 <- paste("pred.train.", methods[i], sep="")
			n.pred2 <- paste("pred.test." , methods[i], sep="")
			n.alpha <- paste("alpha."     , methods[i], sep="")
			n.prior <- paste("prior."     , methods[i], sep="")
			n.err   <- paste("error."     , methods[i], sep="")
			n.stop  <- paste("stop."      , methods[i], sep="")
			n.m     <- paste("m."         , methods[i], sep="")
			n.w     <- paste("w."         , methods[i], sep="")

			# Calculate priors
			assign(n.prior, rbind(get(n.prior),c(sum(get(n.w)[y==0,m]), 
                                            sum(get(n.w)[y==1,m]))))

			# Priors and Weights stopping criteria
			if(methods[i]%in%c("LDA","QDA","NB","SVM")){
				if(m>1 & (sum(get(n.prior)[m,]==get(n.prior)[m-1,])==2)){
					assign(n.stop, TRUE)
					cpt <- cpt + 1 
				}
			}else{
				if(m>1 & (sum(get(n.w)[m,]==get(n.w)[m-1,])==nrow1) ){
					assign(n.stop, TRUE)
					cpt <- cpt + 1 
				}
			}

			if(!get(n.stop)){
				
				# Results from training and testing
				assign(n.res, train.and.test.model(train = data.train, 
                                                           method = methods[i],
							                 test = data.test,
                                                           resp.var.name = "pred", 
							                 weights = get(n.w)[,m],
                                                           priors = get(n.prior)[m,],
							                 knn.par = knn.par,
                                                           svm.par = svm.par, 
							                 nnet.par = nnet.par,
                                                           nbtree = nbtree,
                                                           prob=F))

				# Store predictions on training and test set
				assign(n.pred1, cbind(get(n.pred1), matrix(get(n.res)$pred.train)))
				if(!is.null(test)) assign(n.pred2, cbind(get(n.pred2), 
							                       matrix(get(n.res)$pred.test)))

				# Calculation of training (weighted) error
				disagree <- get(n.pred1)[,m]!=y
				assign(n.err, rbind(get(n.err), sum(get(n.w)[disagree,m])))

				# Update alpha & weights if 0 < error < 0.5
				if((get(n.err)[m]<0.5) & (get(n.err)[m]!=0)){
					# Calculate updated weight for next iteration
					assign(n.alpha, rbind(get(n.alpha), log((1-get(n.err)[m])/get(n.err)[m])))
					updatedW <- get(n.w)[,m]
					updatedW[disagree] <- updatedW[disagree] * exp(get(n.alpha)[m])
					assign(n.w, cbind(get(n.w), updatedW/sum(updatedW)))
					assign(n.m, get(n.m)+1)
				}
				else{
					assign(n.stop, TRUE)
					cpt <- cpt + 1 
				}
			}
		}
		if(cpt==length(methods)) break
	}
	cat("\n")

	## Vote of models
	# Decision function
	f1 <- function(x){ifelse(x[1]>=x[2],1,0)}
	f2 <- function(x){x[1]/x[2]}

	for(i in 1:nb.meth){
		# Names objects
		n.pred1 <- paste("pred.train."   , methods[i], sep="")
		n.pred2 <- paste("pred.test."    , methods[i], sep="")
		n.w     <- paste("w.", methods[i], sep="")
		n.vote1 <- paste("vote.train."   , methods[i], sep="")
		n.vote2 <- paste("vote.test."    , methods[i], sep="")
		n.prob2 <- paste("prob.test."    , methods[i], sep="")
		n.err1  <- paste("error.train."  , methods[i], sep="")
		n.err2  <- paste("error.test."   , methods[i], sep="")
		n.conf  <- paste(paste("conf."   , methods[i],sep=""),".test", sep="")
		n.index <- paste(paste("index."  , methods[i],sep=""),".test", sep="")
		n.m     <- paste("m."            , methods[i], sep="")

		# Weighted vote on training (if not perfect model at iteration 1)
		if(sum(get(n.pred1)[,1]!=y)!=0){
			alpha1 <- apply(get(n.pred1)[,1:get(n.m)] %*% get(n.alpha),1,sum)
			alpha0 <- apply((1-get(n.pred1))[,1:get(n.m)] %*% get(n.alpha),1,sum)
			assign(n.vote1, apply(cbind(alpha1, alpha0), 1, f1))
		}
		else{
			assign(n.vote1, get(n.pred1)[,1])
		}
		assign(n.err1, 1- compute.indexes(confusion.matrix(get(n.vote1),
                                              as.integer(y)))[1,1])

		# Prepare output object
		assign(methods[i], c(get(methods[i]),list("pred.train"=get(n.vote1), 
                                                      "error.train"=get(n.err1))))

		# Write in a file if required
		if(file) write.table(get(n.err1), paste(file.name,"train_err.txt",sep=""), 
                                              row.names=F, col.names=F)

		# Store information on testing (if testing)
		if(!is.null(test)){
			# Sum of alpha of each class, weighted vote and posterior probability as
			# the ratio: sum(weights of chosen class)/sum(all alpha)
			if(sum(get(n.pred1)[,1]!=y)!=0){
				alpha1 <- apply(get(n.pred2)[,1:get(n.m)] %*% get(n.alpha),1,sum)
				alpha0 <- apply((1-get(n.pred2))[,1:get(n.m)] %*% get(n.alpha),1,sum)
				assign(n.vote2, apply(cbind(alpha1, alpha0), 1, f1))
				assign(n.prob2, apply(cbind(alpha1, alpha1+alpha0), 1, f2))
			}
			else{
				assign(n.vote2, get(n.pred2)[,1])
				assign(n.prob2, get(n.vote2))
			}

			# Calculate vote, confusion matrix, indeces and error
			assign(n.conf, confusion.matrix(get(n.vote2),as.integer(data.test$pred)))
			assign(n.index, compute.indexes(get(n.conf)))
			assign(n.err2, 1 - get(n.index)[1,1])

			# Write in a file if required
			if(file){
				write.table(get(n.prob2), paste(file.name,"test_prob.txt",sep=""), 
                                                            row.names=F, col.names=F)
				write.table(get(n.err2), paste(file.name,"test_err.txt",sep=""), 
                                                            row.names=F, col.names=F)
				write.table(get(n.index), paste(file.name,"test_indexes.txt",sep=""), 
                                                            row.names=F, col.names=F)
			}

			# Prepare output object
			assign(methods[i], c(get(methods[i]),list("pred.test"=get(n.vote2),
                                                            "prob.test"=get(n.prob2),
                                                            "conf.test"=get(n.conf),
                                                            "indexes.test"=get(n.index),
                                                            "error.test"=get(n.err2))))
		}
	}

	#---------------------------------- Output ---------------------------------#
	## Creation of output object
	out <- NULL
	for(i in 1:nb.meth){
		assign("out", c(out, get(methods[i])) )
	}
	if(nb.meth>1) names(out) <- methods
	return(out)
}

