# Combine results from function OneSVM.param.cross.clim
# i.e make vote of ensemble of one-class SVM

# Author  :   Gwénaël Leday
# Updated :   March 2009



OneSVM.param.cross.clim.read2 <- function(err, path=NULL, file=T){

	d <- dim(err)
	bool <- err[,4]<.5
	cpt <- 0

	## Calculation of vote and proportion
	for(j in which(bool)){
cat(j,"\n")
		vect <- NULL
		nn <- (length(dir())-3)/nrow(err)
    nn2<-1:nn
		nn2<-nn2[-4]
		for(i in nn2){
			name.file <- paste(paste(c("pred","model",err[j,1]*100,err[j,2]*100,i),
                               sep="",collapse="_"),".txt",sep="")
			if(!is.null(path)) name.file <- paste(c(path,name.file),sep="",collapse="/")
			vect <- rbind(vect, read.table(name.file, header=F))
		}

		vect <- as.numeric(as.matrix(vect))
		if(j==1){
			vote <- vect
		}
		else{
			vote <- vote + vect
		}
		cpt <- cpt + 1
	}
	out <- cbind(vote, vote / cpt)
	colnames(out) <- c("Vote", "%")
	if(file) write.table(out,"votes.txt", row.names=F)
	return(out)
}


