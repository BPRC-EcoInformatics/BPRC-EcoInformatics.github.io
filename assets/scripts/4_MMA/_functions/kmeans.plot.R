# Plotting relevant absence clusters, closest representative absence locations from kmeans clustering

# Author  :   Tak Ikeda
# Updated :   January 2010


kmeans.plot <- function(sp.name, bg, thresh=0, save.fig, file.out, prev=50, nb.var) {

	# Start to measure process time
	ptm1 <- proc.time()

	# Load data
	setwd(paste(main_wd,species.name,"_dist_data_OCSVM",sep="/"))   #####UT modified sep/
	load(paste(sp.name,"_kmeans_undersampling_object_thresh_",thresh,"_prev_",prev,"_",bg,nb.var,".RData",sep=""))
	species_abs <- read.table(paste(sp.name,"_data_absences_thresh_",thresh,"_",bg,nb.var,".txt",sep=""),header=T)
	prevalence <- length(data.abs.kmeans$size)
	n_points_cluster <- data.abs.kmeans$size 
	which_cluster <- data.abs.kmeans$cluster
	n_absences_prob0 <- length(data.abs.kmeans$cluster)
	centroid_mat <- data.abs.kmeans$centers 

cat("Data loaded \n")

	# Cluster colors plotted, k=prevalence
	X11(width=15,height=10)
	map('world', fill = TRUE, col = "white")
	#plot(border.no.validation.country, cex=0.2, xlab="Longitude", ylab="Latitude")
	points(species_abs$xcoord,species_abs$ycoord,col=colors()[1+which_cluster],pch=16,cex=0.2)
	savePlot(paste(sp.name,"_cluster_colours_for_relevant_absences_thresh_",thresh,"_prev_",prev,"_",bg,nb.var,".png",sep=""), type="png")

cat("Clusters plotted \n")

	# Where are the centroids closest to 
	# calculate eucl_dist based on the variables

	min_eucl_dist <- eucl_dist <- NULL; temp <- Inf; which.m <- 1
	out <- matrix(NA,prevalence,2)
	colnames(out) <- list("min_dist","index") # index is the point of species_abs that is closest

cat("Calculating minimum distances... \n")

	for (k in 1:prevalence) {
		species_abs_k <- species_abs[which_cluster==k,]
		for (m in 1:n_points_cluster[k]) {
			eucl_dist <- sqrt(sum((centroid_mat[k,]-species_abs_k[m,-(1:4)])^2))
			min_eucl_dist <- min(eucl_dist,temp)
			if (eucl_dist < temp) which.m <- as.integer(rownames(species_abs_k[m,]))
			temp <- min_eucl_dist
		}
	
		out[k,] <- c(min_eucl_dist,which.m )
cat("k =",k,"\n")
		min_eucl_dist <- eucl_dist <- NULL; temp <- Inf; which.m <- 1
	}

	if (file.out) write.table(out,file=paste(sp.name,"_centroid_closest_data_thresh_",thresh,"_prev_",prev,"_",bg,nb.var,".txt",sep=""),row.names=F)

cat("Minimum distance done \n")


	# Plot closest data points to centroids on map
	#plot(border.no.validation.country, cex=0.2, xlab="Longitude", ylab="Latitude")
	map('world', fill = TRUE, col = "white")	
	points(species_abs$xcoord[out[,2]],species_abs$ycoord[out[,2]],cex=0.5,pch=16,col=3)
	if (save.fig) savePlot(paste(sp.name,"_closest_points_to_centroids_thresh_",thresh,"_prev_",prev,"_",bg,nb.var,".png",sep=""), type="png")
	write.table(species_abs[out[,2],], paste(sp.name,"_kmeans_abs_thresh_",thresh,"_prev_",prev,"_",bg,nb.var,".txt",sep=""), row.names=F)
		
	# Plot presences and points closest to representative absences
	#plot(border.no.validation.country, cex=0.2, xlab="Longitude", ylab="Latitude")
	map('world', fill = TRUE, col = "white")		
	points(species$xcoord,species$ycoord,cex=0.5,pch=16,col=2)
	points(species_abs$xcoord[out[,2]],species_abs$ycoord[out[,2]],cex=0.5,pch=16,col=3)
	legend(-150,-40,col=c(2,3),legend=c("Presence","Absence"),pch=16)
	if (save.fig) savePlot(paste(sp.name,"_closest_points_to_centroids_with_pres_thresh_",thresh,"_prev_",prev,"_",bg,nb.var,".png",sep=""), type="png")

	# calculate total process time
	ptm2 <- proc.time()
	ptm <- (ptm2 - ptm1)[3]

	# Print process time
	cat("\nTime:\n")
		print(ptm)

}
