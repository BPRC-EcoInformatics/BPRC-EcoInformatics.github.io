# Variable selection based on Wilcoxon Cross Validation method

varSelWilcox <- function(sp.data,B) {

	# sp.data = presence/absence data for species
	# B = number of repetitions

	# Error handling
	if(is.null(sp.data)) stop("No input data!")
	if(is.null(B) | B<2) stop("B must be greater than 1!")

	# Set process time
	ptm <- proc.time()[3]

	# WilcoxCV
	data.x <- sp.data[,1:(ncol(sp.data)-1)]
	data.y <- sp.data[,ncol(sp.data)]
	my.n <- length(data.y)
	my.ntest <- round(my.n/5)
	my.split <- generate.split(niter=B, n=my.n, ntest=my.ntest)
	my.wilcox <- wilcox.selection.split(x=data.x,y=data.y,split=my.split,algo="new",pvalue=TRUE)
	my.wilcox.os <- my.wilcox$ordering.split
	my.wilcox.p <- my.wilcox$pvalue.split
	my.wilcox.pmax <- round(apply(my.wilcox.p,2,max),3) # is max pvalue less than 0.05?
	my.wilcox.pmean <- round(apply(my.wilcox.p,2,mean),3) # is mean pvalue less than 0.05?
	my.wilcox.psd <- round(apply(my.wilcox.p,2,sd),3) # is mean pvalue less than 0.05?
	m <- nrow(my.wilcox.os)
	n <- ncol(my.wilcox.os)
	#cat("Wilcox CV done \n")

	# Initialization
	prop.v <- rank.v <- output <- NULL
	rank.v <- matrix(NA,n,2)
	
	# Get rankings
	for (i in 1:n) {
		count.v <- NULL
		uniq.v <- unique(my.wilcox.os[,i])
		for (j in 1:length(uniq.v)) {
			count.v[j] <- sum(my.wilcox.os[,i]==uniq.v[j])
		}
		rank.v[i,] <- c(colnames(data)[uniq.v[which.max(count.v)]],max(count.v)/m)
	}
	
	rank.v.m <- cbind(as.double(rank.v[,2]),
			as.double(my.wilcox.pmax[rank.v[,1]]),
			as.double(my.wilcox.pmean[rank.v[,1]]),
			as.double(my.wilcox.psd[rank.v[,1]]))
	rownames(rank.v.m) <- rank.v[,1]
	colnames(rank.v.m) <- c("Proportion","Max_P","Mean_P","SD_p")
	colnames(my.wilcox.os) <- colnames(data.x)

	# Output, basic info
	output$ordering.split <- my.wilcox.os
	output$pvalue.split <- my.wilcox.p
	output$n.rep <- B
	output$time <- proc.time()[3]-ptm

	# Output, ranked variables
	output$ranked.variables <- rank.v.m
	
	# Output of variable sets:
	# with last variable with proportion greater than 0.75
		output$set.prop <- sort(rownames(rank.v.m[1:(min(which((rank.v.m[,1]>0.75)==F))-1),])) 
	# variables with max p less than 0.05
		output$set.pmax <- sort(rownames(rank.v.m[rank.v.m[,2]<0.05,])) 
	# variables with mean p less than 0.05
		output$set.pmean <- sort(rownames(rank.v.m[rank.v.m[,3]<0.05,])) 

	return(output)
}	
		