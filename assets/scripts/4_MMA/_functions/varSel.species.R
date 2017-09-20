varSel.species <- function(species,thresh,preva,bg,nb.var) {

# Variable selection with RF and stepwise

	data <- read.table(paste(species,"_data_pres_abs_thresh_",thresh,"_prev_",preva,"_",bg,nb.var,".txt",sep=""), header=T)
	data.x <- data[,1:(ncol(data)-1)]
	data.y <- data[,ncol(data)]	

	# STEP 1: Using random forest
	#----------------------------#
	selRF <- varSelRF(data.x, as.factor(data.y), mtryFactor=1,
	ntree=10000, ntreeIterat=5000, vars.drop.frac = NULL,vars.drop.num = 1)

	# Variables selection
	varSel <- selRF$selected.vars

	# Variables selection
	selhis <- selRF$selec.history

	# Database with selected normalized variables of random forest
	dataSel <- as.data.frame(cbind(data.x[,varSel],data.y))
	colnames(dataSel)[length(varSel)+1] <- species

	# Only predictors
	#dataSel.x <- dataSel[,colnames(dataSel)!=species.code]
	# Only response variable
	#dataSel.y <- dataSel[,species.code]

	# STEP 2: Stepwise selection
	#---------------------------#
	model <- step(glm(as.formula(paste(species,"~.",sep="")),
	data=dataSel, family="binomial"),direction="both", k=2)

	# Variable selection 2
	varSel2 <- colnames(model$model)[2:length(model$model)]

	# Database with final selected variables
	#dataSel2 <- as.data.frame(cbind(dataSel[,varSel2],data.y))

	# Only predictors
	#dataSel2.x <- dataSel[,varSel2]
	# Only response variable
	#dataSel2.y <- dataSel[,species.code]

	vs <- NULL
#	vs$model <- model
	vs$selhis <- selhis
	vs$rf_initial <- varSel
	vs$step_initial <- varSel2


	return(vs)
}

