# Summary of posterior probabilities

summary.post.prob <- function(prob, B, y){

	# Functions to apply
	f  <- function(x) abs(x-0.5)
	f1 <- function(x) mean(x, na.rm=T)
	f2 <- function(x) sd(x, na.rm=T)

	# Mean and Standard deviation of posterior probabilities for each observation
	post.prob <- cbind(apply(prob,  1,f1), apply(prob,  1,f2))

	# Add column of predictions obtained by averaging posterior probabilities 
	post.prob <- cbind(post.prob,  as.integer(post.prob[,1] > 0.5))

	# Add column for uncertainty of predictions
	post.prob <- cbind(post.prob,   as.numeric((post.prob[,1] + post.prob[,2])>0.5  
                                               & (post.prob[,1] - post.prob[,2])<0.5))

	# Add column for uncertainty's weight of predictions
	# (weight = |mean-0.5|/sd)
	weight <- rep(1,length(y))
	bool.uncert <- post.prob[,4] == 1
#	if(sum(as.integer(bool.uncert),na.rm=T)!=0){ # Edited 26/8/10 by Tak
	if(sum(as.integer(bool.uncert))!=0){ 
		weight[bool.uncert] <- abs(post.prob[bool.uncert,1] - 0.5)
     		weight[bool.uncert] <- weight[bool.uncert] / post.prob[bool.uncert,2]
		post.prob <- cbind(post.prob,   weight)
	}
	else{
		post.prob <- cbind(post.prob,   weight)
	}

	# Same number of digits
	post.prob <- round(post.prob,6)

	# Set names of columns
	colnames(post.prob)   <- c("Mean","Sd","Pred","Uncert","UncertW")

	# Transform prob matrices as list so that it can directly be used with ROCR package
	prob.list <- NULL
	for(i in 1:B){
		prob.list <- c(prob.list, list(prob[!is.na(prob[,i]),i]))
	}

	# Confusion matrix with averaged posterior probabilities
	conf   <- confusion.matrix(post.prob[,3] ,y)

	# Associated measures of performance
	index   <- compute.indexes(conf)

	# Summary of uncertainty for output
	uncert <- sum(post.prob[,4]) / length(post.prob[,4])

	# Output object
	output <- NULL
	output$post.prob   <- post.prob
	output$prob.list   <- prob.list
	output$conf        <- conf
	output$indexes     <- index
	output$uncertainty <- uncert

	return(output)
}