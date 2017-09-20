# Compute bootstraps, .632 and .632+ errors a specified model

# Author  :   Gwénaël Leday
# Updated :   May 2009



errors.632plus.model <- function(data,
                                 method=NULL,
                                 y,
                                 err.mat,
                                 B,
                                 resp.var.name, 
                                 knn.par=NULL, 
                                 svm.par=NULL, 
                                 nnet.par=NULL,
                                 nbtree=NULL,
                                 bagging=NULL,
                                 boosting=NULL){

	# Formula
	form <- as.formula(paste(resp.var.name," ~ .",sep=""))

	# Function Mean removing NA
	f <- function(x)mean(x, na.rm=T)

	# Underestimated resubstitution error and
	# Overestimated boot error
	err.mat[(B+1),] <- apply(err.mat[1:B,],2,f)

	# .632 error
	err.mat[(B+2),1] <- sum(err.mat[(B+1),]*c(0.368,0.632))

	# Gamma
	if(is.null(bagging) & is.null(boosting)){
		model.pred <- train.and.test.model(train=data,
                                               method=method,
                                               resp.var.name="pred",
                                               knn.par=knn.par,
                                               svm.par=svm.par,
                                               nbtree=nbtree,
                                               nnet.par=nnet.par)$pred.train
	}
	if(!is.null(bagging)){
		model.pred <- train.and.test.bagged.model(train=data,
                                                      methods=method,
                                                      B=bagging,
                                                      resp.var.name="pred",
                                                      knn.par=knn.par,
                                                      svm.par=svm.par,
                                                      nnet.par=nnet.par)$pred.train
	}
	if(!is.null(boosting)){
		model.pred <- train.and.test.boosted.model(train=data,
                                                       methods=method,
                                                       M=boosting,
                                                       resp.var.name="pred",
                                                       knn.par=knn.par,
                                                       svm.par=svm.par,
                                                       nnet.par=nnet.par)$pred.train
	}
	gamma <- (sum(y)*sum(model.pred))/length(y)+
               (sum(1-y)*sum(1-model.pred))/length(y)

	# No-Information Rate
	NIR <- (err.mat[(B+1),2] - err.mat[(B+1),1])/(gamma - err.mat[(B+1),1])

	# Ajusted weight
	w <- 0.632/(1-0.368*NIR)

	# .632+ error
	err.mat[(B+3),1] <- w*err.mat[(B+1),2] + (1 - w)*err.mat[(B+1),1]

	return(err.mat)
}