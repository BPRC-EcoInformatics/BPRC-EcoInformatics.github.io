# to create overlay of predictions for several species 
# with votes for predictions greater than a threshold of choice. 


overlay <- function(species, methods, p.lab, thresh, col.set, save.fig=T, save.out=T, save.file=T) {

	# --------------------- # Info # --------------------- #

	# Arguments:
	# species	= a vector of names, ex: c("ac","atm")
	# methods	= corresponding vector of methods, ex: c("svm","knn")
	# p.lab	= corresponding vector of parameter values, ex: c("10_0.1", "4")
	# thresh	= between (0,1)
	# col.set	= "heat" for yellow to red, "bluered" for blue to red
	# save.fig	= logical, T for saving figure 
	# save.out	= logical, T for saving Rdata
	# save.file	= logical, T to output ArcGIS file of votes

	# Value:
	# returns overlaid data of list of species

	# --------------------- # Initialization # --------------------- #

	# Start to measure process time
	ptm1 <- proc.time()

	nb.species <- length(species)
	f <- function(x) sum(x>=55537)==0 # distinguish missing values
#	sp.lab <- NULL
#	for (i in 1:nb.species) sp.lab <- paste(sp.lab,species[i],sep="_")
	sp.lab <- paste(nb.species,"species",sep="")
	

	# Set colors for plots
	if (col.set=="heat") 	col.vec <- rev(heat.colors(nb.species))
	if (col.set=="bluered")	col.vec <- c("#27408B",rev(rainbow(nb.species+nb.species-1)))[-(2:(nb.species+1))]


	# Default vector of length bioclim_nz, L=417606
	votes <- rep(0,nrow(bioclim_nz)) # empty vector of votes, 0 to start, length of L

	# --------------------- # Main Part # --------------------- #

	# Take votes of species predictions above threshold
	for (i in 1:nb.species) {
#for (i in 2:2) { # i<-2
		# Load predictions vector, length<L, and make T/F vector >= thresh
		setwd(paste("C:\\MMA\\",species.name,"\\_prediction_maps\\_predict_NZ",sep=""))
		cat("Counting votes for ",species[i],", ",methods[i],", ",p.lab[i]," ..... ", sep="")
		load(paste("NZ_prediction_",species[i],"_",methods[i],"_",p.lab[i],".Rdata",sep=""))
		res <- train.and.test.model.predict$prob.test

		# Read in data of species presence/absence and create boolean vector
		setwd(paste("C:\\MMA\\",species.name,"\\_dist_data_OCSVM",sep=""))
		data <- read.table(paste(species[i],"_data_pres_abs_best.txt",sep=""), header=TRUE)
		bool <- apply(bioclim_nz[,colnames(data)[1:(ncol(data)-1)]],1,f) 


		# Update votes vector, in correct rows, by adding new votes 
cat(i, length(bool), length(res), sum(bool))
		votes[bool] <- votes[bool] + as.integer(res>=thresh)
cat(" Done \n")
		# Reset prediction
		rm(train.and.test.model.predict)
	}

cat("\nPreparing plot of summed votes... \n")

	# --------------------- # Output # --------------------- #
	# Make plot of summed votes 
#	win.graph(width=48,height=70)
	plot(bioclim_nz[,1],bioclim_nz[,2],cex=0.2,pch=16,col="lightcyan1", 
		xlab="Longitude", ylab="Latitude") # 0 species count

	for (i in 1:nb.species) 
		points(bioclim_nz[votes==i,1],bioclim_nz[votes==i,2],cex=0.2,
			pch=16,col=col.vec[i]) # 1 or greater species count
	if(nb.species<10) legend(168,-35, legend=0:nb.species,col=c("lightcyan1",col.vec),pch=15,bg="white")
	if(nb.species>10) legend(168,-34, legend=seq(0,nb.species,2),col=c("lightcyan1",col.vec[seq(0,nb.species,2)+1]),pch=15,bg="white")

	# Saving output in R workspace
	if (save.out)  {
		setwd("R:\\Intelligent systems for biosecurity\\INVASIVE_SPP\\_prediction_maps\\_overlay")
		info.matrix <- cbind(species,methods,p.lab)
		save(votes,info.matrix,file=paste("NZ_Votes_",sp.lab,"_",thresh,".Rdata",sep=""))
	}

	# Saving figure
	if (save.fig)  {
		setwd("R:\\Intelligent systems for biosecurity\\INVASIVE_SPP\\_prediction_maps\\_overlay")
		savePlot(paste("NZ_Votes_",sp.lab,"_thresh_",thresh,"_",col.set,".png",sep=""), type="png")
	}	

	# Saving ArcGIS file
	if (save.file) {	
		setwd("R:\\Intelligent systems for biosecurity\\INVASIVE_SPP\\_prediction_maps\\_overlay")
		dat.vote <- cbind(bioclim_nz[,1:2],votes)
		colnames(dat.vote) <- c("xcoord","ycoord","count")
		write.table(dat.vote, paste("NZ_Votes_",sp.lab,"_thresh_",thresh,"_lonlat_counts.txt",sep=""), row.names=F)
	}

	# Calculate total process time
	ptm2 <- proc.time()
	ptm <- (ptm2 - ptm1)[3]

	# Print process time
	cat("\nTime:\n")
		print(ptm)

	return(votes)
}