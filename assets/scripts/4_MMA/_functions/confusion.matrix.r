# Compute the confusion matrix
# (in a binary classification problem only)

# Author : Gwénaël Leday
# Date   : March 2008


confusion.matrix <- function(class.predict,class){

	## Error handling
	if(length(class.predict)!=length(class)) stop("Vectors must have the same length")
	
	## Initialization
	confusion <- matrix(0,ncol=2,,nrow=2)

	## Fill matrix
	confusion[1,1] <- sum(class.predict[class==1] ==  1)
	confusion[1,2] <- sum(class.predict[class==1] ==  0)
	confusion[2,1] <- sum(class.predict[class==0] ==  1)
	confusion[2,2] <- sum(class.predict[class==0] ==  0)

	## Output
	rownames(confusion) <- c("TRUE 1","TRUE 0")
	colnames(confusion) <- c("PRED 1","PRED 0")
	return(round(confusion,3))
}
