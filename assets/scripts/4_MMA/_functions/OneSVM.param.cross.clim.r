# Parametrization of one-class SVM
# i.e. explore space of parameters

# Note: This function is made for high dimensional climate data
# instead of storing predictions it splits the vector and write it.
# An associated function read files in order to combine results

# Author  :   Gwenael Leday
# Updated :   March 2009



OneSVM.param.cross.clim <- function(data1.x,
                                    data0.x,
                                    k=5,
                                    nb.models=100,
                                    prior=NULL,
                                    file=T,
                                    save.out=F,bg=bg){

	## Error handling
	if(is.character(nb.models) & nb.models!="all"){
		stop("Error: wrong character for 'nb.models'!")
	}
	if(!is.null(prior)){
		prior <- as.matrix(prior)
		if(ncol(prior)!=2) stop("Error: wrong input prior!")
	}

	## Current CPU Time
	ptm1 <- proc.time()

	## Data
	data.x     <- as.matrix(data1.x)
	dataTest.x <- as.matrix(data0.x)
	DimTest    <- dim(dataTest.x)

	## Space of parameters for one-class SVM
	gam <- seq(0.01,0.99,0.01)
	nu  <- seq(0.01,0.99,0.01)
#	gam <- seq(0.1,0.9,0.1)
#	nu  <- seq(0.1,0.9,0.1)
	indeces <- matrix(1:(length(gam)*length(nu)),length(gam),length(nu))
	indx <- NULL
	if(!is.null(prior)){
		for(i in 1:nrow(prior)) indx <- c(indx,indeces[prior[i,1]*100,prior[i,2]*100])
		indeces[!(indeces %in% indx)] <- NA
	}
	list.mod <- indeces[!is.na(indeces)]
	if(nb.models=="all") nb.models <- length(list.mod)

	## Error matrices
	error <- matrix(NA, nb.models, 5)
	colnames(error) <- c("Sigma", "nu", "Err.Train", "Err.Test", "Nb_Pred_1")
	error2 <- NULL

	## Predictions
	pred  <- NULL

	## Splits
	split <- c(1,1+20000*1:floor(DimTest[1]/20000), DimTest[1]+1)
	nsplit <- length(split)
	

	## Body
	nb <- nb2 <- cpt <- 0
	while((nb<nb.models) & (cpt!=9801)){
		cat("\n")
		cat(cpt, "\t")
		# Random choice of parameters
		continue <- T
		while(continue){
			# Random selection of a 'gam' and 'nu' value
			x <- which(indeces==sample(list.mod,1), arr.ind=T)
			g <- gam[x[1]]
			n <- nu[x[2]]

			# Check if the combination of parameters has already been done
			if(!is.na(indeces[x[1],x[2]])){
				indeces[x[1],x[2]] <- NA
				list.mod <- indeces[!is.na(indeces)]
				continue <- F
			}
		}

		# Evaluation of one-class SVM by k-fold cross-validation
		
		# PErhAPS ERROR HERE #
		model <- OneSVM.CrossVal(data1=data.x, data2=dataTest.x, NbGroup = k, OneSVM.par = c(g, n))
		if(model$error.test < .5){
			nb <- nb + 1
			error[nb,1] <- g
			error[nb,2] <- n
			error[nb,3] <- model$error.train
			error[nb,4] <- model$error.test
			error[nb,5] <- model$nb.pred
# Added nsplit for different row numbers of background data
# Condition related to nrow(data0.x) not needed
#			if (nrow(data0.x)>500001) {
				for(i in 1:(nsplit-1)){
					name.file <- paste(paste(c("pred","model",bg,g*100,n*100,i),sep="",collapse="_"),".txt",sep="")
	   				write.table(model$pred[split[i]:(split[i+1]-1)],
                  	                  name.file,append=F,row.names=F, col.names=F)
				}
#			}
# Removed 
#			if (nrow(data0.x)<500001) {
#				for(i in 1:4){
#					name.file <- paste(paste(c("pred","model",g*100,n*100,i),sep="",collapse="_"),".txt",sep="")
#	   				write.table(model$pred[split[i]:(split[i+1]-1)],
#                 	 			name.file,append=F,row.names=F, col.names=F)
#				}
#			}
			write.table(error, paste("model_in_",bg,nb.var,".txt",sep=""))
			cat("ok \t nb.models = ",nb,"\n")
		}else{
			nb2 <- nb2 + 1
			error2 <- rbind(error2, c(g, n, model$error.train, model$error.test, model$nb.pred))
			write.table(error2, paste("model_out_",bg,nb.var,".txt",sep=""))
		}
		cpt <- cpt + 1
	}

	## Output
	out <- NULL
	rownames(error) <- paste("model_", 1:nb.models, sep="")
	error.ord <- error[order(error[,4]),]
#	if (is.null(dim(error2))) error2 <- t(as.matrix(error2))
	colnames(error2) <- colnames(error)
	rownames(error2) <- paste("model_", 1:nrow(error2), sep="")
	out$models.in  <- error.ord
	out$models.out <- error2
	if(file){
		write.table(error.ord, paste("model_in_",bg,nb.var,".txt",sep=""))
		write.table(error2, paste("model_out_",bg,nb.var,".txt",sep=""))
	}

	## Current CPU Time
	ptm2 <- proc.time()
	ptm <- ptm2 - ptm1
	out$time <- ptm
	cat("\nTime:\n")
	print(ptm)

	if(save.out) save(out, file=paste("one_class_svm_out_",bg,nb.var,".RData",sep=''))
	return(out)
}
