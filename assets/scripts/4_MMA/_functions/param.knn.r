# Parametrization of the k-nearest neighbours

# Author  :   GwénaE Leday
# Updated :   May 2009


param.knn <- function(x,
                      y,
                      knn.rep = 100,
                      resampling = "boot",
                      nrep = 200,
                      train.frac = .75,
                      kfold = 10,
					  out.object = T){


	#-------------------- Initialization ------------------#
	## Data
	data <- cbind(x,y)
	ncol <- ncol(data)
	colnames(data)[ncol] <- "pred"

	## Error matrices
	err.train <- err.test <- rep(list(matrix(NA,17,knn.rep)),nrep)
	if(resampling=="boot"){
		error.632 <- noInfRate <- error.632plus <- matrix(NA,17,knn.rep)
		gamma <- rep(NA,17)
	}else{
		error <- matrix(NA,17,knn.rep)
	}
	temp1 <- temp2 <- matrix(0,17,knn.rep)

	# Output
	true.list <- test.id.list <- bal <- NULL
	name.root <- paste(c("param","knn", resampling), sep="", collapse="_")

	#------------------------ Body ------------------------#
	# Current CPU Time
	ptm1 <- proc.time()

	for(n in 1:nrep){
		cat(ifelse(resampling=="boot","b = ","rep = "),n,"\n")
		# Generate samples
		resample <- resample(data, method=resampling, train=train.frac, cv=kfold)
		data.train <- resample$train
		if(resampling!="cv"){
			data.test <- resample$test
			bal <- c(bal, sum(data.train[,ncol])/nrow(data.train))
			test.id.list <- c(test.id.list, list(resample$test.ind))
			true.list <- c(true.list, list(data.test[,ncol]))
		}else{
			cv.groups <- resample$groups
			test.id.list <- c(test.id.list, list(unlist(cv.groups)))
			true.list <- c(true.list, list(data.train[,ncol]))
		}
		
		# Compute error for different k
		for(j in 1:knn.rep){
			for(i in 4:20){
				if(resampling=="cv") err.cv <- NULL
				if(resampling!="cv") cv.groups <- NA
				for(l in 1:length(cv.groups)){
					if(resampling!="cv"){
						train <- data.train
						test  <- data.test
					}
					else{
						train <- data.train[-cv.groups[[l]],]
						test  <- data.train[cv.groups[[l]],]
					}
					# Train and test
					res <- train.and.test.model(train=train, method="KNN",
							test=test, test2=train, resp.var.name="pred",
							knn.par=i, prob=F,sp=FALSE)

					# Compute errors
					error.train <- res$error.train
					error.test <- res$error.test

					if(resampling=="cv"){
						err.cv <- rbind(err.cv, c(error.train, error.test))
					}
				}

				if(resampling=="cv"){
					err.train[[n]][i-3,j] <- mean(err.cv[,1])
					err.test[[n]][i-3,j]  <- mean(err.cv[,2])
				}else{
					err.train[[n]][i-3,j] <- error.train
					err.test[[n]][i-3,j]  <- error.test
				}
			}
		}
		temp1 <- temp1 + err.train[[n]]
		temp2 <- temp2 + err.test[[n]]
	}

	# Averaged errors
	temp1 <- temp1/nrep
	error <- temp2 <- temp2/nrep

	# BOOT
	if(resampling=="boot"){
		error.632 <- 0.632*temp2 + 0.368*temp1
		for(i in 4:20){
			all.pred <- as.integer(knn(as.data.frame(data[,-ncol]),as.data.frame(data[,-ncol]),
					data[,ncol],k=i,prob=F))-1
			gamma[i-3] <- (sum(y)*sum(all.pred))/length(y)+
				   (sum(1-y)*sum(1-all.pred))/length(y)
		}
		noInfRate <- (temp2 - temp1)/(matrix(rep(gamma,knn.rep),17,knn.rep) - temp1)
		w <- 0.632/(1-0.368*noInfRate)
		err632plus <- w*temp2 + (1 - w)*temp1
		error <- err632plus
	}

	# Current CPU Time
	ptm2 <- proc.time()
	ptm <- ptm2 - ptm1

	#---------------------- Output ---------------------#
	out <- NULL
	# Calculate k that minimize error (.632+ or cv.err or test error)
	print(error)
	f1 <- function(x){x==min(x)}
	f2 <- function(x){which(x==TRUE)}
	k.min.err <- unlist(apply(apply(error,2,f1),2,f2)+3)

	# Plot for decision
	k.table <- matrix(0,1,17)
	colnames(k.table) <- 4:20
	for(i in 1:17) k.table[i] <- sum(k.min.err==(i+3))

	## Plot
	X11(width = 9, height = 9)
	par(mar = c(5, 4, 2, 2) + 0.1)
	barplot(k.table/sum(k.table),col=grey(0.8),xlab="k",ylab="Frequency", ylim=c(0,1))
	#savePlot(paste(species.name,"_",name.root,"_plot.jpg",sep=""),type="jpg")
	savePlot(paste(species.name,"_",name.root,"_plot_",bg,nb.var,"_",reduction,".png",sep=""),type="png")
	#savePlot(paste(species.name,"_",name.root,"_plot.pdf",sep=""),type="pdf")

	## Output object
	out$k.min.err <- k.min.err
	out$k.distribution <- k.table
	rownames(error) <- as.character(4:20)
	colnames(error) <- paste("rep",1:knn.rep,sep="")
	out$errors <- error

	# Trace
	out$trace.rep <- test.id.list

	# List of response variable vectors of non-bootstrapped observations
	out$true.list <- true.list

	# Time
	out$time <- ptm
	cat("Time:\n\n")
	print(ptm)
	cat("\n\n")

	if(out.object) save(out, file = paste(species.name,"_",name.root,"_output_object_",bg,nb.var,"_",reduction,".Rdata",sep=""))
	return(out)
}
