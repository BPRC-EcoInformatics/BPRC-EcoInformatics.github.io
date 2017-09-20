# Prediction with a specified model
# This function harmonizes output and avoid redundancy of code

# Author  :   Gwenael Leday
# Updated :   May 2009

# Edited  :   Tak Ikeda
# Updated :   Jan 2010
# Reason  :   To handle large test sets in predictions. 
# 		'6000' comes from the max size of a species dataset.


predict.model <- function(model,
                          method, 
                          train=NULL, 
                          test, 
                          resp.var.name, 
                          knn.par=NULL, 
                          svm.par=NULL, 
                          nnet.par=NULL,
				  Wts=NULL, 
                          prob=T){

	## Error handling
	if(is.null(model) & method!="KNN") stop("Error: no input model")
	if(is.null(resp.var.name)) stop("Error: no input name of response variable")

	## Data
	test <- as.data.frame(test)
	if(!is.null(train)) train <- as.data.frame(train)

	## Formula
	form <- as.formula(paste(resp.var.name," ~ .",sep=""))

	## Predictors
	predictors <- !names(train) %in% resp.var.name

	## Predict test set with the chosen model
	switch(method,
		"LDA"  = {
			if (nrow(test)<6000) {
				pred <- predict(model,newdata=test)
				model.pred  <- as.numeric(pred$class)-1
				if(prob) model.prob  <- pred$posterior[,2]
			}
			if (nrow(test)>6000) { # testing on worldclim or nz
				pred1 <- predict(model,newdata=test[     1:300000,])
				pred2 <- predict(model,newdata=test[300001:nrow(test),])
				pred <- c(pred1,pred2)
				model.pred <- c(as.numeric(pred1$class)-1,as.numeric(pred2$class)-1)
				if(prob) model.prob  <- c(pred1$posterior[,2],pred2$posterior[,2])
			}
		},
		"QDA"  = {
			if (nrow(test)<6000) {
				pred <- predict(model,newdata=test)
				model.pred  <- as.numeric(pred$class)-1
				if(prob) model.prob  <- pred$posterior[,2]
			}
			if (nrow(test)>6000) { # testing on worldclim or nz
				pred1 <- predict(model,newdata=test[     1:300000,])
				pred2 <- predict(model,newdata=test[300001:nrow(test),])
				pred <- c(pred1,pred2)
				model.pred <- c(as.numeric(pred1$class)-1,as.numeric(pred2$class)-1)
				if(prob) model.prob  <- c(pred1$posterior[,2],pred2$posterior[,2])
			}
		},
		"LOG"  = {
			if (nrow(test)<6000) model.prob  <- predict(model,newdata=test, type="response")
			if (nrow(test)>6000) { # testing on worldclim or nz
				pred1 <- predict(model,newdata=test[     1:300000,], type="response")
				pred2 <- predict(model,newdata=test[300001:nrow(test),], type="response")
				model.prob <- c(pred1,pred2)
			}
			model.pred  <- as.numeric(model.prob>0.5)
		},
		"NB"   = {
			if (nrow(test)<6000) model.prob  <- predict(model,newdata=test)$posterior[,2]
			if (nrow(test)>6000) { # testing on worldclim or nz
				pred1 <- predict(model,newdata=test[     1:300000,], type="response")
				pred2 <- predict(model,newdata=test[300001:nrow(test),], type="response")
				#model.prob <- c(pred1,pred2) #error
        model.prob  <- c(pred1$posterior[,2],pred2$posterior[,2])
        
        
			}
			#model.pred  <- as.numeric(model.prob>0.5)  #error see MMA2step
      	    model.pred<-as.numeric(model.prob>0.5)
		},
		"CART" = {
			if (nrow(test)<6000) {
				model.pred  <- as.numeric(predict(model, newdata=test, type="class"))-1
				if(prob) model.prob  <- predict(model , newdata=test, type="prob")[,2]
			}
			if (nrow(test)>6000) { # testing on worldclim or nz
				pred1 <- as.numeric(predict(model, newdata=test[     1:300000,], type="class"))-1
				pred2 <- as.numeric(predict(model, newdata=test[300001:nrow(test),], type="class"))-1
				model.pred <- c(pred1,pred2)
				if(prob) {
					model.prob1  <- predict(model , newdata=test[     1:300000,], type="prob")[,2]
					model.prob2  <- predict(model , newdata=test[300001:nrow(test),], type="prob")[,2]
					model.prob   <- c(model.prob1,model.prob2)	
				}
			}
		},
		"STUMP" = {
			if (nrow(test)<6000) {
				model.pred  <- as.numeric(predict(model, newdata=test, type="class"))-1
				if(prob) model.prob  <- predict(model , newdata=test, type="prob")[,2]
			}
			if (nrow(test)>6000) { # testing on worldclim or nz
				pred1 <- as.numeric(predict(model, newdata=test[     1:300000,], type="class"))-1
				pred2 <- as.numeric(predict(model, newdata=test[300001:nrow(test),], type="class"))-1
				model.pred <- c(pred1,pred2)
				if(prob) {
					model.prob1  <- predict(model , newdata=test[     1:300000,], type="prob")[,2]
					model.prob2  <- predict(model , newdata=test[300001:nrow(test),], type="prob")[,2]
					model.prob   <- c(model.prob1,model.prob2)	
				}
			}
		},
		"CTREE"= {
			if (nrow(test)<6000) model.prob  <- predict(model ,newdata=test)
			if (nrow(test)>6000) { # testing on worldclim or nz
				pred1 <- predict(model,newdata=test[     1:300000,])
				pred2 <- predict(model,newdata=test[300001:nrow(test),])
				model.prob <- c(pred1,pred2)
			}
			model.pred  <- as.numeric(model.prob>0.5)
		},
		"RF"= {
			if (nrow(test)<6000) model.prob  <- predict(model, newdata=test, type="prob")[,2]
			if (nrow(test)>6000) { # testing on worldclim or nz
				pred1 <- predict(model,newdata=test[     1:300000,], type="prob")[,2]
				pred2 <- predict(model,newdata=test[300001:nrow(test),], type="prob")[,2]
				model.prob <- c(pred1,pred2)
			}
			model.pred  <- as.numeric(model.prob>0.5)
		},
		"KNN" = {
			if (nrow(test)<6000) {
				model.pred  <- as.integer(knn(train[,predictors], test[,predictors],
							as.factor(train[,resp.var.name]), 
							k=knn.par, prob=F))-1
				if(prob){
					model.prob <- attr(knn(train[,predictors], test[,predictors],
							as.factor(train[,resp.var.name]), 
							k=knn.par, prob=T),"prob")
					model.prob[model.pred==0] <- 1 - model.prob[model.pred==0]
				}
			}
			if (nrow(test)>6000) { # testing on worldclim or nz
				model.pred1  <- as.integer(knn(train[,predictors], test[1:300000,predictors],
							as.factor(train[,resp.var.name]), 
							k=knn.par, prob=F))-1
				model.pred2  <- as.integer(knn(train[,predictors], test[300001:nrow(test),predictors],
							as.factor(train[,resp.var.name]), 
							k=knn.par, prob=F))-1
				model.pred <- c(model.pred1,model.pred2)
			
			
				if(prob){
	
					model.prob1 <- attr(knn(train[,predictors], test[1:300000,predictors],
							as.factor(train[,resp.var.name]), 
							k=knn.par, prob=T),"prob")
					model.prob2 <- attr(knn(train[,predictors], test[300001:nrow(test),predictors],
							as.factor(train[,resp.var.name]), 
							k=knn.par, prob=T),"prob")
					model.prob <- c(model.prob1,model.prob2)
					model.prob[model.pred==0] <- 1 - model.prob[model.pred==0]
				}
			}
		},
		"SVM"  = {
			if (nrow(test)<6000) { # testing on training set
				print('i am ok in SVM')
				#print(model)
				#print(test)
				model.pred  <- predict(model, test, type="response")
				print('I model pred ok')
				if(prob) {#print(predict(model, test, type="probabilities"))
				model.prob  <- predict(model, test, type="probabilities")[,2]}
				print("DONE SVM")
			}
			if (nrow(test)>6000) { # testing on worldclim or nz
				pred1 <- predict(model,test[     1:300000,], type="response")
				pred2 <- predict(model,test[300001:nrow(test),], type="response")
				model.pred  <- c(pred1,pred2)
				if(prob) {
					pred1 <- predict(model,test[     1:300000,], type="probabilities")[,2]
					pred2 <- predict(model,test[300001:nrow(test),], type="probabilities")[,2]
					model.prob <- c(pred1,pred2)
				}
			}
		},
		"NNET" = {
			if (nrow(test)<6000) { # testing on training set
				model.pred  <- as.numeric(predict(model, newdata=test, type="class"))
				if(prob) model.prob  <- predict(model, newdata=test, type="raw",Wts=Wts)
			}
			if (nrow(test)>6000) { # testing on worldclim or nz
				pred1 <- as.numeric(predict(model, newdata=test[     1:300000,], type="class"))
				pred2 <- as.numeric(predict(model, newdata=test[300001:nrow(test),], type="class"))
				model.pred  <- c(pred1,pred2)
				if(prob) {
					pred1 <- predict(model,newdata=test[     1:300000,], type="raw",Wts=Wts)
					pred2 <- predict(model,newdata=test[300001:nrow(test),], type="raw",Wts=Wts)
					model.prob <- c(pred1,pred2)
				}
			}


		}

	)

	print("output")
	# Output
	ifelse(prob,out <- list("prob"=model.prob, "pred"=model.pred),out <- list("pred"=as.integer(model.pred>0.5)))
	return(out)
}
