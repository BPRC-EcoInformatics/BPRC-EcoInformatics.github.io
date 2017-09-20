# Make bootstrap files for Condor directory "boot"

# Load functions
bootstrap.condor <- function(species.name,nreps) {
	data <- read.table(paste(species.name,"_data_pres_abs_best.txt",sep=""), header=TRUE)
	data <- data[data$pres==1,]
	setwd(paste("C:\\MMA\\",species.name,"\\_condor\\boot",sep=""))
	boot.mat <- matrix(NA,dim(data)[1],nreps)
	colnames(boot.mat) <- c(paste(rep("Boot",nreps),1:nreps,sep=""))

	for (i in 1:nreps) {
		boot.mat[,i] <- !resample(data,method="boot")$test.ind
		write.table(as.data.frame(boot.mat[,i]),paste("Boot",i,".txt",sep=""),
			row.names=F,col.names=paste("Boot",i,sep=""))
	}
}
