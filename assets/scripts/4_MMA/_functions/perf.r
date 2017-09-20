# Calculate ROC and Lift Chart Area Under Curve (AUC)
# and associated objects for plotting

# Author  :   Gwénaël Leday
# Updated :   May 2009

perf <- function(pp){

	# Basic computation with the confusion matrix
	conf <- confusion.matrix(as.numeric(pp[,1]>0.5), pp[,2])
	TP <- conf[1,1]
	FN <- conf[1,2]
	FP <- conf[2,1]
	TN <- conf[2,2]
	total <- sum(conf)
	pred    <- prediction(pp[,1], pp[,2])

	# Calculation and output
	out <- NULL
	out$perfROC   <- performance(pred, measure="tpr", x.measure="fpr")
	out$perfLift   <- performance(pred, measure="tpr", x.measure="rpp")
	out$aucROC  <- unlist(attributes(performance(pred, measure="auc"))$y.values)
	out$aucLift <- ((1/total)*(((TP+FN)^2/2)+(TP+FN)*(TN+FP)*out$aucROC))/(TP+FN)

	return(out)
}