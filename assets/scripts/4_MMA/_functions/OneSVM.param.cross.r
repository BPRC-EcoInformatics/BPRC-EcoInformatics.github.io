# Parametrization of one-class SVM
# i.e. explore space of parameters

# Author  :   Gwénaël Leday
# Updated :   March 2009



OneSVM.param.cross <- function(data1.x, data0.x, k=5, nb.models=100, file=T, save.out=F){

	# Data
	data.x     <- as.matrix(data1.x)
	dataTest.x <- as.matrix(data0.x)
	DimTest    <- dim(dataTest.x)

	## Space of parameters for one-class SVM
	gam <- seq(0.01,0.99,0.01)
	nu  <- seq(0.01,0.99,0.01)
	indeces <- matrix(1:(length(gam)*length(nu)),length(gam),length(nu))
	if(is.character(nb.models) & nb.models!="all") stop("Error: wrong character for 'nb.models'")
	if(nb.models=="all") nb.models <- length(gam)*length(nu)

	# Radial basis kernel
	error <- matrix(NA, nb.models, 5)
	colnames(error) <- c("Sigma", "nu", "Err.Train", "Err.Test", "Nb_Pred_1")
	error2 <- NULL

	# Predictions
	pred  <- NULL

	# Current CPU Time
	ptm1 <- proc.time()

	# Body
	nb <- nb2 <- cpt <- 0
	while((nb<nb.models) & (cpt!=9801)){
		cat("\n")
		cat(cpt, "\t")
		# Random choice of parameters
		continue <- T
		while(continue){
			# Random selection of a 'gam' value
			g   <- sample(gam, 1)
			ind1 <- which(gam==g)

			# Random selection of a 'nu' value
			n <- sample(nu, 1)
			ind2 <- which(nu==n)

			# Check if the combination of parameters has already been done
			if(!is.na(indeces[ind1,ind2])){
				indeces[ind1,ind2] <- NA
				continue <- F
			}
		}

		# Evaluation of one-class SVM by k-fold cross-validation
		model <- OneSVM.CrossVal(data.x, dataTest.x, NbGroup = k, OneSVM.par = c(g, n))
		if(model$error.test < .5){
			nb <- nb + 1
			error[nb,1] <- g
			error[nb,2] <- n
			error[nb,3] <- model$error.train
			error[nb,4] <- model$error.test
			error[nb,5] <- model$nb.pred
			pred <- cbind(pred, model$pred)
			cat("ok \t nb.models = ",nb,"\n")
		}else{
			nb2 <- nb2 + 1
			error2 <- rbind(error2, c(g, n, model$error.train, model$error.test, model$nb.pred))
		}
		cpt <- cpt + 1
	}

	## Calculation of proportions
	prop <- apply(pred,1,sum) / nrow(pred)
	prop <- as.matrix(prop)
	rownames(prop) <- rownames(dataTest.x)

	# Current CPU Time
	ptm2 <- proc.time()
	ptm <- ptm2 - ptm1

	## Output
	out <- NULL
	rownames(error) <- paste("model_", 1:nb.models, sep="")
	rownames(pred) <- rownames(dataTest.x)
	colnames(pred) <- paste("model_", 1:nb.models, sep="")
	colnames(error2) <- colnames(error)
	rownames(error2) <- paste("model_", 1:nrow(error2), sep="")
	out$proportions <- as.matrix(prop[order(prop,decreasing=T),])
	error.ord <- error[order(error[,4]),]
	out$models.in  <- error.ord
	out$models.out <- error2
	out$pred <- pred
	out$time <- ptm
	if(file){
		write.table(error.ord, "model_in.txt")
		write.table(prop, "proportions.txt", col.names=F)
	}
	cat("\nTime:\n")
	print(ptm)

	if(save.out) save(out, "one_class_svm_out.Rdata")
	return(out)
}