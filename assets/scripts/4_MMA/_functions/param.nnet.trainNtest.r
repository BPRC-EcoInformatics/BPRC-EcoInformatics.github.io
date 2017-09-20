# param.nnet.trainNtest trains and tests neural networks for a given bootstrap

# dat = = database for the species in data.frame format (if variable selection is required, it has to be done before)
# bootVect = boolean vector of the bootstrap
# n = number of times a model has to be repeated
# name = name of the output text files

param.nnet.trainNtest <- function(dat,bootVect,n=100, name=NULL){

#------------------------------ Initialization -----------------------------#
	# Determine boolean vector for bootstrapped and non-bootstrapped data
	boot <- bootVect[,1]
	notboot <- !bootVect

	# Redo database so that response variable is factor and
	# splt databse into bootstrapped and non-bootstrapped observations
	data <- cbind(dat[,1:(ncol(dat)-1)],as.factor(dat[,ncol(dat)]))
	ncol <- ncol(data)
	colnames(data)[ncol] <- "pred"
	data.boot <- data[boot,]
	data.notboot <- data[notboot,]
	# Keep response variable not in a factor format for both databases
	y.boot <- dat[boot,ncol(dat)]
	y.not.boot <- dat[notboot,ncol(dat)]

	# Storage matrix of errors (resubstitution and test)
	error <- matrix(NA,length(s)*length(m)*length(d),3,
			dimnames=list(NULL,c("Size","Maxiter","Decay")))
	for(i in 1:n){
		error <- cbind(error,matrix(NA,length(s)*length(m)*length(d),2,
			dimnames=list(NULL,c(paste("Resub",i,sep=""),paste("Boot",i,sep="")))))
	}

	# Parameters of NNET
	s <- 2:5
	m <- c(500,750,1000,1250,1500)
	d <- c(0.01,0.001,0.0001)


#----------------------------------- Body ----------------------------------#
	# Initialize counter
	cpt <- 1

	# Loops, one for each parameter 
	for(i in 1:length(s)){
		for(j in 1:length(m)){
			for(k in 1:length(d)){
				cat("i=",i," ,j=",j," ,k=",k,"\n")

				# Storage combination of parameters
				error[cpt,1] <- i
				error[cpt,2] <- j
				error[cpt,3] <- k

				# Repeat training and test of a model
				l <- 1
				while(l<(2*n)){
					# Training
					model <- nnet(pred~.,data=as.data.frame(data.boot), trace=F,size=s[i], maxit=m[j], decay=d[k])
					# Test model on training data
					model.pred <- as.numeric(predict(model,newdata=as.data.frame(data.boot), type="class"))
					# Test model on test data
					model.pred.test <- as.numeric(predict(model,newdata=as.data.frame(data.notboot), type="class"))

					# If predictions and true vectors are comparable then compute both
					# resubstitution errors and bootstrap error
					if( (sum(model.pred)!=0) && (max(model.pred)==max(y.boot)) &&
					  	(sum(model.pred.test)!=0) && (max(model.pred.test)==max(y.not.boot))){
						conf <- confusion.matrix(model.pred,y.boot)
						error[cpt,3 + l] <- (conf[1,2] + conf[2,1])/sum(conf)
						conf <- confusion.matrix(model.pred.test,y.not.boot)
						error[cpt,3 + l + 1] <- (conf[1,2] + conf[2,1])/sum(conf)
					}
					l <- l + 2
				}
				cpt <- cpt + 1
			}
		}
	}


#---------------------------------- Output ---------------------------------#
	write.table(error,paste(name,"_results.txt",sep=""),row.names=F)
	return(error)
}