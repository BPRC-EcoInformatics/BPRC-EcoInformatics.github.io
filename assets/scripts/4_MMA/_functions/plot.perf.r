# Plot the performance of models : ROC plot and Lift chart

# Author  :   GwÈnaÅE Leday
# Updated :   May 2009


plot.perf <- function(pp = NULL,
                      names = NULL,
                      colors = NULL,
                      models.plots = 1,
                      file=NULL
                     ){

	## Error handling
	# on 'pp'
	if(is.null(pp))  stop("Error: no input")
	if(!is.list(pp)) stop("Error: 'pp' has to be a list")
	# Number of input models
	nb.mod <- length(pp)
	# on 'names'
	if(length(names)!=nb.mod) stop("Error: Length of 'names' has to 
                                           be the number of input models")
	# on 'models.plots'
	if(!is.vector(models.plots)) stop("Error: 'models.plots' has to be a vector")
	if(!is.numeric(models.plots)) stop("Error: 'models.plots' has to be numeric")
	if(max(models.plots)>nb.mod) stop("Error: max of 'models.plots' is greater than number of models")


	## Default names and colors
	if(is.null(names)) names <- paste("mod",1:nb.mod,sep="")
	if(is.null(colors)) colors <- rainbow(nb.mod)

	## Performances
	perfs <- NULL
	aucROC <- aucLift <- matrix(NA,nb.mod,1)
	rownames(aucROC) <- rownames(aucLift) <- names
	for(i in 1:nb.mod){
		# Calculate performance values for both kind of plots
		perfs <- c(perfs, list(perf(pp[[i]])))

		# Merge ROC/Lift AUC of models
		aucROC[i] <- perfs[[i]]$aucROC
		aucLift[i] <- perfs[[i]]$aucLift
	}

	## Rank models compared to their area under curves
	aucROC.ord <- as.matrix(aucROC[order(aucROC[,1],decreasing=T),])
      aucLift.ord <- as.matrix(aucLift[order(aucLift[,1],decreasing=T),])
	colnames(aucROC) <- colnames(aucROC.ord) <- "auc ROC"
	colnames(aucLift) <- colnames(aucLift.ord) <- "auc Lift"

  #------------------------ Plots -----------------------#
  if(nb.mod<=10){
	# Coordinates for display AUC
	ycoor <- seq(0.5,0.05,-0.05)

	## Lift charts
	for(j in sort(unique(models.plots))){
		# File's name that has to be saved
		if(j==nb.mod){
			name <- paste(file,"final_plot_lift_all",sep="")
		}
		else{
			name <- paste(paste(file,"final_plot_lift_best",sep=""),j,sep="")
		}

		# Plot j best models
		par(mar = c(5, 4, 2, 2) + 0.1)
		plot(seq(0,1,length=50),seq(0,1,length=50),type="l", lwd=1, col="white",
		xlab="Positive Rate", ylab="True Positive Rate")
		abline(h=seq(0,1,by=0.1), v=seq(0,1,by=0.1), col = grey(0.8), lty=3)
		for(i in 1:nb.mod){
			if(i==1) lines(seq(0,1,length=50),seq(0,1,length=50),type="l",lwd=1.9,lty=1)
			if(names[i] %in% rownames(aucLift.ord)[1:j]){
				lines(unlist(perfs[[i]]$perfLift@x.values),
     	                        unlist(perfs[[i]]$perfLift@y.values),
     	                        lwd=2, col=colors[i])
				text(0.8, ycoor[i], paste(paste(paste("AUC", names[i]),"= "),
     	                        round(aucLift[i],3)), font=2, col=colors[i])
			}
		}

		# Save plot in different format
		savePlot(paste(name,"_",bg,nb.var,"_",reduction,".png",sep=""),type="png")
		#savePlot(paste(name,".pdf",sep=""),type="pdf")
	}

	# New window
	X11()

	## ROC plots
	for(j in sort(unique(models.plots))){
		# File's name that has to be saved
		if(j==nb.mod){
			name <- paste(file,"final_plot_roc_all",sep="")
		}
		else{
			name <- paste(paste(file,"final_plot_roc_best",sep=""),j,sep="")
		}

		# Plot j best models
		par(mar = c(5, 4, 2, 2) + 0.1)
		plot(seq(0,1,length=50),seq(0,1,length=50),type="l", lwd=1, col="white",
		xlab="False Positive Rate", ylab="True Positive Rate")
		abline(h=seq(0,1,by=0.1), v=seq(0,1,by=0.1), col = grey(0.8), lty=3)
		for(i in 1:nb.mod){
			if(i==1) lines(seq(0,1,length=50),seq(0,1,length=50),type="l",lwd=1.9,lty=1)
			if(names[i] %in% rownames(aucROC.ord)[1:j]){
				lines(unlist(perfs[[i]]$perfROC@x.values),
     	                        unlist(perfs[[i]]$perfROC@y.values),
     	                        lwd=2, col=colors[i])
				text(0.8, ycoor[i], paste(paste(paste("AUC", names[i]),"= "),
     	                        round(aucROC[i],3)), font=2, col=colors[i])
			}
		}

		# Save plot in different format
		savePlot(paste(name,"_",bg,nb.var,"_",reduction,".png",sep=""),type="png")
		#savePlot(paste(name,".pdf",sep=""),type="pdf")
	}
  }
  else{
	warning("Cannot plot more than 10 models")
  }

	#------------------------ Output -----------------------#
	out <- NULL
	out$aucROC <- aucROC.ord
	out$aucLift <- aucLift.ord
	return(out)
}
