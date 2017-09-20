# ROC plots for resampling methods

# Author  :   GwénaE Leday
# Updated :   May 2009



plot.roc.all <- function(prob, true, name.file=NULL){

	## Error handling
	if(!is.list(prob) | !is.list(true)) stop("Error: list format required")
	if(!is.character(name.file)) stop("Error: 'name.file' has to be a string!")
	name <- paste(name.file, "_", sep="")
	
	# Number of repetitions/bootstraps
	n <- length(prob)

	# AUC
	auc <- matrix(NA,n+2,1)
	rownames(auc) <- c(paste("iter",1:n,sep=""),"Mean","Sd")
	colnames(auc) <- "AUC"

	# ROC AUC
	pred <- prediction(prob,true)
	perf <- performance(pred, measure="tpr", x.measure="fpr")
	auc[1:n]  <- unlist(attributes(performance(pred, measure="auc"))$y.values)#erase yvalues
	auc[n+1] <- mean(auc[1:n])
	auc[n+2] <- sd(auc[1:n])
	aucROC <- round(auc[n+1],3)

	# Plot1: all ROC curves in grey plus the averaged one in black
	par(mar=c(5, 5, 2, 2) + 0.1)
	plot(perf,col=grey(0.8))
	plot(perf, avg="vertical", spread.estimate="none", lwd=2, col="black",
		xlab="False Positive Rate", ylab="True Positive Rate",add=T)
	abline(h=seq(0,1,by=0.1), v=seq(0,1,by=0.1), col = grey(0.8), lty=3)
	lines(seq(0,1,length=50),seq(0,1,length=50),  lwd=2, col="black")
	text(0.7,0.35, paste("AUC = ", aucROC), font=2, col="black")
	savePlot(paste(name,"plot_roc_all1_",bg,nb.var,"_",reduction,".png",sep=""),type="png")
	#savePlot(paste(name,"plot_roc_all1.pdf",sep=""),type="pdf")

	# Plot2: averaged ROC curve in black plus boxplots at regular space
	par(mar=c(5, 5, 2, 2) + 0.1)
	plot(perf, avg="vertical", spread.estimate="boxplot", lwd=2, col="black",
		xlab="False Positive Rate", ylab="True Positive Rate")
	abline(h=seq(0,1,by=0.1), v=seq(0,1,by=0.1), col = grey(0.8), lty=3)
	lines(seq(0,1,length=50),seq(0,1,length=50),  lwd=2, col="black")
	text(0.7,0.35, paste("AUC = ",aucROC), font=2, col="black")
	savePlot(paste(name,"plot_roc_all2_",bg,nb.var,"_",reduction,".png",sep=""),type="png")
	#savePlot(paste(name,"plot_roc_all2.pdf",sep=""),type="pdf")

	return(auc)
}
