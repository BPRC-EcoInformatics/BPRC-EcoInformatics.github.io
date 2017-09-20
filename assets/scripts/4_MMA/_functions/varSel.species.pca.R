varSel.species.pca <- function(species.name,thresh,preva,bg,nb.var) {

# Variable selection with RF and stepwise

	data <- read.table(paste(species.name,"_data_pres_abs_thresh_",thresh,"_prev_",preva,"_",bg,nb.var,".txt",sep=""), header=T)
	data.x <- data[,1:(ncol(data)-1)]
	data.y <- data[,ncol(data)]	


	# Run PCA
	selPCA <- dudi.pca(data.x, center = TRUE, scale = TRUE, scannf = FALSE, nf = 3)

	# Prints of the statistics of inertia in a one-table analysis 
	# Contributions are printed in 1/10000 and the sign is the sign of the coordinate
	inertia<-inertia.dudi(selPCA,col.inertia=TRUE)
	var.contrib <- inertia$col.abs # absolute contributions of the decomposition of inertia for the column
		
	tot_contrib<-var.contrib[,1]*selPCA$eig[1]+var.contrib[,2]*selPCA$eig[2]+var.contrib[,3]*selPCA$eig[3]

	names_var<-rownames(var.contrib)
	res_all<-data.frame(names_var,tot_contrib)
	res_all[with(res_all, order(tot_contrib)), ]

	threshold_3comp <- 80
	res<-res_all[which(res_all[,2]>threshold_3comp),]
	colnames(res)<-c("Variable","Contribution")
	dim(res)

	if(dim(res)[1]==0) {
		threshold_3comp <- median(tot_contrib)
		res<-res_all[which(res_all[,2]>threshold_3comp),]
	}
	
	return(res)
}

