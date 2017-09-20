# Evaluation of a one-class SVM model by cross-validation (CV)
# Note: CV is derived from package bootstrap

# Author  :   Gwena?E Leday
# Updated :   March 2009



OneSVM.CrossVal <- function(data1, data2, NbGroup = 10, OneSVM.par){

	# Error Handling
	x <- as.matrix(data1)
	x2 <- as.matrix(data2)
	n <- nrow(data1)
	ngroup <- trunc(NbGroup)
	if( ngroup < 2){
		stop ("The number of groups should be greater than or equal to 2")
	}
	if(ngroup > n){
		stop ("The number of groups should be less than or equal to the number of observations")
	}
	if(!is.vector(OneSVM.par)) stop("Parameters of One-Class SVM should be a vector")
	if(length(OneSVM.par)!=2) stop("Exactly two parameters for One-Class SVM")

  
	# Samples
	if(ngroup==n) {groups <- 1:n; leave.out <- 1}
	if(ngroup<n){
		leave.out <- trunc(n/ngroup);
		o <- sample(1:n)
		groups <- vector("list",ngroup)
		for(j in 1:(ngroup-1)){
			jj <- (1+(j-1)*leave.out)
			groups[[j]] <- (o[jj:(jj+leave.out-1)])
		}
		groups[[ngroup]] <- o[(1+(ngroup-1)*leave.out):n]
	}

	# One-Class SVM training and test
	err.test  <- rep(NA,ngroup)
	err.train <- rep(NA,ngroup)
	for(j in 1:ngroup){
		model.svm 	<- ksvm(x[-groups[[j]],], type="one-svc", kernel="rbfdot",
					kpar=list(sigma=OneSVM.par[1]), nu=OneSVM.par[2])
		pred 		<-  predict(model.svm,x[groups[[j]],])
		err.train[j] <- model.svm@error
		err.test[j]  <- (length(pred) - sum(pred)) / length(pred)	
	}

	model.svm 	<- ksvm(x, type="one-svc", kernel="rbfdot",
				kpar=list(sigma=OneSVM.par[1]), nu=OneSVM.par[2])

# divide into 6 prediction vectors to avoid vector allocation limit 
# if data2 has 600000 rows.
# x2 is a part of wordclim, so dim(x2) is always between 500000 and 600000.
	if (nrow(data2)>500001) {
		pred1 <- predict(model.svm,x2[     1:100000,])
		pred2 <- predict(model.svm,x2[100001:200000,])
		pred3 <- predict(model.svm,x2[200001:300000,])
		pred4 <- predict(model.svm,x2[300001:400000,])
		pred5 <- predict(model.svm,x2[400001:500000,])
		pred6 <- predict(model.svm,x2[500001:nrow(x2),])
		pred  <- c(pred1,pred2,pred3,pred4,pred5,pred6)
		nb.pred	<-  sum(pred)
	}
	if (nrow(data2)<500001) {
		pred <- predict(model.svm,x2)
		nb.pred	<-  sum(pred)
	}
	
	# Output
	if(leave.out==1) groups <- NULL
	out <- NULL
	out$error.test 	<- mean(err.test)
	out$error.train 	<- mean(err.train)
	out$nb.pred		<- nb.pred
	out$pred		<- pred
	errors <- matrix(NA, ngroup, 2)
	colnames(errors) 	<- c("Training","Test")
	rownames(errors) 	<- paste("Cross",1:ngroup, sep="")
	out$errors 		<- errors
	
	return(out)
}