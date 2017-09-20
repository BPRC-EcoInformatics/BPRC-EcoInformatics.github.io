# Learning a specified model
# This function harmonizes output and avoid redundancy of code

# Author  :   Gwenael Leday
# Updated :   May 2009


learn.model <- function(data, 
                        resp.var.name, 
                        method=NULL, 
                        svm.par=NULL, 
                        nnet.par=NULL,
                        nbtree=NULL,
				weights=NULL,
                        Wts=NULL, 
                        priors=NULL){

	## Error handling
	if(is.null(data)) stop("Error: no input data")
	if(is.null(resp.var.name)) stop("Error: no input name of response variable")

	## Data
	train <- as.data.frame(data)

	## Formula
	form <- as.formula(paste(resp.var.name," ~ .",sep=""))

	## Predictors
	predictors <- !names(train) %in% resp.var.name

	## Train the chosen model
	switch(method,
		"LDA"  = {
			if(!is.null(priors)) model <- lda(form, data=train, prior=priors)
			if(is.null(priors))  model <- lda(form, data=train)
		},
		"QDA"  = {
			if(!is.null(priors)) model <- qda(form, data=train, prior=priors)
			if(is.null(priors))  model <- qda(form, data=train)
		},
		"LOG"  = {
			if(!is.null(weights)) model <- glm(form, data=train, family="quasibinomial", weights=weights)
			if(is.null(weights))  model <- glm(form, data=train, family="binomial")
		},
		"NB"   = {
			train <- cbind(train[,predictors], as.factor(train[,resp.var.name]))
			names(train)[ncol(train)] <- resp.var.name
			if(!is.null(priors)) model <- NaiveBayes(form, data=train, prior=priors)
			if(is.null(priors))  model <- NaiveBayes(form, data=train)
		},
		"CART" = {
			if(!is.null(weights)) cart <- rpart(form, data=train, method="class",cp=0.00001, weights=weights)
			if(is.null(weights))  cart <- rpart(form, data=train, method="class",cp=0.00001)
			if(dim(cart$cptable)[2]==5){
				cp.opt <- cart$cptable[which.min(cart$cptable[,4]),1]
				model <- prune(cart,cp=cp.opt)
			}else{
				model <- cart
			}
		},
		"STUMP" = {
			model <- rpart(form, data=train, weights=weights, method="class", cp=0.00001, 
                                 control=rpart.control(maxdepth=1))
		},
		"CTREE"= {
			if(!is.null(weights)) weights <- as.integer(weights * 10^6)
			if(!is.null(weights)) model <- ctree(form,data=as.data.frame(train), weights=weights)
			if(is.null(weights))  model <- ctree(form,data=as.data.frame(train))
		},
		"KNN"  = {model <- NULL},
		"RF"  = {model <- randomForest(train[,predictors], as.factor(train[,resp.var.name]), 
                               importance=F, ntree=nbtree)
		},
		"SVM"  = {
			switch(svm.par[1],
				"vanilladot" ={
					if(!is.null(priors)){
						model <- ksvmMinstep(form,train, type="C-svc", kernel="vanilladot", 
                                             prob.model=F, C=as.double(svm.par[2]), class.weights=priors)
					}
					if(is.null(priors)){
						model <- ksvmMinstep(form,train, type="C-svc", kernel="vanilladot", 
                                             prob.model=T,C=as.double(svm.par[2]))
					}
				},
				"rbfdot" ={
					if(!is.null(priors)){
						model <- ksvmMinstep(form, train, type="C-svc", kernel="rbfdot", 
                                             prob.model=F, C=as.double(svm.par[2]), 
                                             kpar=list(sigma=as.numeric(svm.par[3])), class.weights=priors, scaled=F)
					}
					if(is.null(priors)){
						model <- ksvmMinstep(form, train, type="C-svc", kernel="rbfdot", 
                                             prob.model=T, C=as.double(svm.par[2]), kpar=list(sigma=as.numeric(svm.par[3])))
					}
				},
				"polydot" ={
					if(!is.null(priors)){
						model <- ksvmMinstep(form, train, type="C-svc", kernel="polydot", 
                                             prob.model=F,
						         cost=as.integer(svm.par[2]), kpar=list(scale=as.numeric(svm.par[3]), 
						         degree=as.integer(svm.par[4]), offset=as.integer(svm.par[5])), class.weights=priors)
					}
					if(is.null(priors)){
						model <- ksvmMinstep(form, train, type="C-svc", kernel="polydot", 
                                             prob.model=T, cost=as.integer(svm.par[2]), 
                                             kpar=list(scale=as.numeric(svm.par[3]), 
						         degree=as.integer(svm.par[4]), offset=as.integer(svm.par[5])))
					}
				}
			)
		},
		"NNET" = {
			# Response variable hsa to be a factor
			train <- cbind(train[,predictors], as.factor(train[,resp.var.name]))
			names(train)[ncol(train)] <- resp.var.name
			if(!is.null(Wts)) model <- nnet(form,data=as.data.frame(train), trace=F, Wts=Wts,
							size=nnet.par[1], 
							maxit=nnet.par[2], 
							decay=nnet.par[3])
			if(is.null(Wts)) model <- nnet(form,data=as.data.frame(train), trace=F,
							size=nnet.par[1], 
							maxit=nnet.par[2], 
							decay=nnet.par[3])
		}
	)

	# Output
	return(model)
}
