# Resample data: bootstrap, hold-out

# Author  :   Gwénaël Leday
# Updated :   May 2009


resample <- function(data, method="holdout", train=NULL, cv=NULL){

	## Error handling
	#if(method=="holdout" & is.null(train)) stop("Error: train fraction?")
	#if(!is.null(train) & method=="holdout") warning("train.fraction set at 0.75!!")
	#if(train>1 | train<0) stop("Error: probability required!")
	#if(!is.null(train) & (train==0 | train==1)) stop("Error: need train/test set!")

	## Initialization
	nrow <- nrow(data)

	## Boostrap
	if(method=="boot"){
		id  <- sample(nrow, nrow, replace=T)
	}

	## Hold-out
	if(method=="holdout"){
		ind <- sample(c(T,F), nrow, replace=T, prob=c(train, 1-train))
		id  <- which(ind)
	}

	## Cross-Validation
	if(method=="cv"){
		groups  <- cv.samples(n=nrow, k=cv)
	}

	## Train/test sets creation
	if(method=="cv"){
		train <- data
	}
	else{
		train <- data[id,]
		test.ind <- !(1:nrow %in% id)
		test <- data[test.ind,]
	}

	## Output
	out <- list("train"=train)
	if(method!="cv") out <- c(out, list("test"=test))
	if(method=="cv") out <- c(out, list("groups"=groups))
	if(method!="cv") out <- c(out, list("test.ind"=test.ind))
	return(out)
}
