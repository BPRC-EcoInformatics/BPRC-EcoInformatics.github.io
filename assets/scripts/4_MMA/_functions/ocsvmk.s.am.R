ocsvmk.s.am <- function(sp.name, color.choice, nb=10, save.fig, file.out,do.kmeans=T,thresh=0,preva=50) {


# ------------------------ #--# ------------------------ #
# ------------------------ #--# ------------------------ #
# ------------------ Do ALL of OCSVMk ------------------ #
#   gather votes, plot,get best absences, make dataset   #
# ------------------------ #--# ------------------------ #
# ------------------------ #--# ------------------------ #

# INPUT:
# sp.name = character of species name
# color.choice = color of plot, "brown","heat","rainbow","venette","maxent"
# nb = 10, number of spaces. 
# save.fig = logical, to save figures
# file.out = logical, to save output in files
# do.kmeans = logical, T: (default) do kmeans
#		F: no kmeans calculation, just ocsvm plots
# thresh = 0, absences combined from predictions with probability = 0 (default).
#             if thresh = t (100>t>0), it will take marginal probabilities (percentage) < t
# preva = 50 (default), from 0 to 100

# OUTPUT:
# a) maps
# - OCSVM environmental suitability map based on 100 SVMs
# - above, with presences
# - just presences
# - random sample of relevant samples
# b) files
# - models.in
# - models.out
# - final.votes
# - k.means absences
# - presence/absence dataset
# ------------------------ #--# ------------------------ #


# ------------------------------ #
# Read 'votes.txt' of each space #
# ------------------------------ #

for(i in 1:nb){
	setwd(paste("C:\\MMA\\",species.name,"\\_OCSVM_results",sep=""))
	assign(paste("vote",i,sep=""),read.table("votes.txt",header=T))
	assign(paste("model.in.",i,sep=""),read.table("model_in.txt",header=T))
	assign(paste("model.out.",i,sep=""),read.table("model_out.txt",header=T))
	cat(paste("dimensions vote",i,sep=""),":\t",dim(get(paste("vote",i,sep=""))))
	cat("\n")
	}

# --------- #
# Sum votes # 
# --------- #

setwd(paste("C:\\MMA\\",species.name,"\\_OCSVM_results",sep=""))
for(i in 1:nb){
	if(i == 1) {
		final.vote <- get(paste("vote",i,sep=""))[,1]
		model.in <- get(paste("model.in.",i,sep=""))
		model.out <- get(paste("model.out.",i,sep=""))
		}
	else{
		final.vote <- final.vote + get(paste("vote",i,sep=""))[,1]
		model.in <- rbind(model.in, get(paste("model.in.",i,sep="")))
		model.out <- rbind(model.out, get(paste("model.out.",i,sep="")))
		}
	}

# ----------------------- #
# Prepare output and save #
# ----------------------- #

res <- cbind(final.vote, final.vote/100)
colnames(res) <- c("Vote","%")
write.table(res, paste(sp.name,"_final_votes.txt",sep=""), row.names=F)
cat("Saved output step1 \n")

# ---------- #
# World Plot #
# ---------- #

# DATA
setwd(paste("C:\\MMA\\",species.name,"\\_dist_data",sep=""))
species <- read.table(paste(sp.name,"_worldclim.txt",sep=""),header=T)

ind <- species[,2:3]
#inds <- read.table(paste(sp.name,"_data_with_gbif.txt",sep=""),header=T)
ind[,3] <- 1
colnames(ind)[3] <- "pres"
species <- species[which(ind[,3]!=0),]



# FILTER

#f <- function(x) sum(is.na(x))==0
#chosen.v <- (1:19)+39
#bool0 <- apply(species[,chosen.v],1,f)

# ------------------- #
# Save models.in/.out
# ------------------- #

setwd(paste("C:\\MMA\\",species.name,"\\_dist_data_OCSVM",sep=""))
write.table(model.in, paste(sp.name,"_models_in.txt",sep=""), row.names=F)
write.table(model.out,paste(sp.name,"_models_out.txt",sep=""), row.names=F)
cat("Saved output step2 \n")

# -------------- #
# Split and plot #
# -------------- #

spl <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
ind1 <- res[,2]<spl[1]
ind2 <- res[,2]>spl[1] & res[,2]<=spl[2]
ind3 <- res[,2]>spl[2] & res[,2]<=spl[3]
ind4 <- res[,2]>spl[3] & res[,2]<=spl[4]
ind5 <- res[,2]>spl[4] & res[,2]<=spl[5]
ind6 <- res[,2]>spl[5] & res[,2]<=spl[6]
ind7 <- res[,2]>spl[6] & res[,2]<=spl[7]
ind8 <- res[,2]>spl[7] & res[,2]<=spl[8]
ind9 <- res[,2]>spl[8] & res[,2]<=spl[9]
ind10 <- res[,2]>spl[9]

col <- rep("white", dim(res)[1])

if (color.choice == "brown") {
#	col[ind1] <- rgb(207, 207, 207, maxColorValue=255)
	col[ind1] <- rgb(224, 255, 255, maxColorValue=255) # <0.1
	col[ind2] <- rgb(255, 215, 0, maxColorValue=255)
	col[ind3] <- rgb(255, 193, 37, maxColorValue=255)
	col[ind4] <- rgb(255, 180, 34, maxColorValue=255)
	col[ind5] <- rgb(255, 185, 15, maxColorValue=255)
	col[ind6] <- rgb(255, 165, 0, maxColorValue=255)
	col[ind7] <- rgb(205, 104, 57, maxColorValue=255)
	col[ind8] <- rgb(139, 76, 57, maxColorValue=255)
	col[ind9] <- rgb(139, 54, 38, maxColorValue=255)
	col[ind10] <- rgb(139, 37, 0, maxColorValue=255)
	col.vec <- c(rgb(224, 255, 255, maxColorValue=255), 
			rgb(255, 215, 0, maxColorValue=255),
			rgb(255, 193, 37, maxColorValue=255),
			rgb(255, 180, 34, maxColorValue=255),
			rgb(255, 185, 15, maxColorValue=255),
			rgb(255, 165, 0, maxColorValue=255),
			rgb(205, 104, 57, maxColorValue=255),
			rgb(139, 76, 57, maxColorValue=255),
			rgb(139, 54, 38, maxColorValue=255),
			rgb(139, 37, 0, maxColorValue=255))
	}
	
	
if (color.choice == "rainbow") {
#	col[ind1] <- rainbow(13)[10]
	col[ind1] <- rgb(224, 255, 255, maxColorValue=255) # <0.1
	col[ind2] <- rainbow(13)[9]
	col[ind3] <- rainbow(13)[8]
	col[ind4] <- rainbow(13)[7]
	col[ind5] <- rainbow(13)[6]
	col[ind6] <- rainbow(13)[5]
	col[ind7] <- rainbow(13)[4]
	col[ind8] <- rainbow(13)[3]
	col[ind9] <- rainbow(13)[2]
	col[ind10] <- rainbow(13)[1]
	col.vec <- c(rgb(224, 255, 255, maxColorValue=255),rev(rainbow(13)[1:9]))
	}

if (color.choice == "heat") {
#	col[ind1] <- heat.colors(10)[10]
	col[ind1] <- rgb(224, 255, 255, maxColorValue=255) # <0.1
	col[ind2] <- heat.colors(10)[9]
	col[ind3] <- heat.colors(10)[8]
	col[ind4] <- heat.colors(10)[7]
	col[ind5] <- heat.colors(10)[6]
	col[ind6] <- heat.colors(10)[5]
	col[ind7] <- heat.colors(10)[4]
	col[ind8] <- heat.colors(10)[3]
	col[ind9] <- heat.colors(10)[2]
	col[ind10] <- heat.colors(10)[1]
	col.vec <- c(rgb(224, 255, 255, maxColorValue=255),rev(heat.colors(10)))
	}

if (color.choice == "venette") {
	col[ind1] <- rgb(224, 255, 255, maxColorValue=255) # <0.1
	col[ind2] <- rgb( 39,  64, 139, maxColorValue=255)
	col[ind3] <- rgb( 69, 139,   0, maxColorValue=255)
	col[ind4] <- rgb(192, 255,  62, maxColorValue=255)
	col[ind5] <- rgb(255, 255,   0, maxColorValue=255)
	col[ind6] <- rgb(255, 215,   0, maxColorValue=255)
	col[ind7] <- rgb(255, 165,   0, maxColorValue=255)
	col[ind8] <- rgb(255, 127,   0, maxColorValue=255)
	col[ind9] <- rgb(255,   0,   0, maxColorValue=255)
	col[ind10]<- rgb(205,   0,   0, maxColorValue=255)
	col.vec <- c(rgb(224, 255, 255, maxColorValue=255), 	# <0.1
			rgb(39, 64, 139, maxColorValue=255), # royalblue4 0.1-0.2
			rgb(69, 139, 0, maxColorValue=255), # chartreuse4 0.2-0.3
			rgb(192, 255, 62, maxColorValue=255), # olivedrab1 0.3-0.4
			rgb(255, 255, 0, maxColorValue=255), # yellow 0.4-0.5
			rgb(255, 215, 0, maxColorValue=255), # gold1 0.5-0.6
			rgb(255, 165, 0, maxColorValue=255), # orange 0.6-0.7
			rgb(255, 127, 0, maxColorValue=255),# darkorange1 0.7-0.8
			rgb(255, 0, 0, maxColorValue=255), # red 0.8-0.9
			rgb(205, 0, 0, maxColorValue=255)) # red3 0.9-1
	}
# honeydew1		240, 255, 240
# honeydew2 	224, 238, 224
# lightcyan1 	224, 255, 255
# lightcyan2 	209, 238, 238
# lightsteelblue1 202, 225, 255
# powderblue 	176, 224, 230
# slategrey1	198, 226, 255
# paleturquoise	187, 255, 255
# azure2		224, 238, 238
legend.loc <- c(-90,-20)
cex1 <- 0.2
pch1 <- 16
	
win.graph(width=55,height=70)
plot(x[,1],x[,2],cex=cex1,pch=pch1,col=col, xlab="Longitude", ylab="Latitude")
legend(legend.loc[1],legend.loc[2], legend=c(paste("<",spl[1]), 
					paste(paste(spl[1],"-"),spl[2]),
					paste(paste(spl[2],"-"),spl[3]),
					paste(paste(spl[3],"-"),spl[4]),
					paste(paste(spl[4],"-"),spl[5]),
					paste(paste(spl[5],"-"),spl[6]),
					paste(paste(spl[6],"-"),spl[7]),
					paste(paste(spl[7],"-"),spl[8]),
					paste(paste(spl[8],"-"),spl[9]),
					paste(">",spl[9])), 
					col=col.vec,
					pch=15,bg="white", cex=0.7)
#points(border.no.nz[,2], border.no.nz[,3], cex=0.2, xlab="Longitude", ylab="Latitude")
#title(paste("Prediction for World",sp.name,sep=","))
if(save.fig) savePlot(paste(sp.name,"_World_ensemble_oneSVM_",color.choice,"_no_pres.png",sep=""), type="png")

# -------------------------------- #
# Plot presences on top of ocsvm
# -------------------------------- #

points(species[,2], species[,3], col="magenta",pch=20, cex=0.5)
if(save.fig) savePlot(paste(sp.name,"_World_ensemble_oneSVM_",color.choice,"_pres.png",sep=""), type="png")

# ------------------- #
# Plot just presences 
#-------------------- #

plot(border[,2], border[,3], cex=0.2, xlab="Longitude", ylab="Latitude")
#plot(worldclim[!bool2,2],worldclim[!bool2,3],cex=0.01, lwd=0.05,xlab="Longitude", ylab="Latitude")
points(species[,2], species[,3], col="magenta",pch=20, cex=0.5)
if(save.fig) savePlot(paste(sp.name,"_World_presences.png",sep=""), type="png")

# ------------------------ #
# Make database for ArcGIS
# ------------------------ #

#if(file.out) {
#	dat <- cbind(x[,1],x[,2],res)
#	colnames(dat) <- c("xcoord", "ycoord","Vote","Prob")
#	write.table(dat, paste(sp.name,"_data_ArcGIS.txt",sep=""), row.names=F)
#}

# ---------------------------------- #
# Make database of Relevant Absences
# ---------------------------------- #

if (thresh == 0) bool.abs <- res[,2] == 0
if (thresh > 0)  bool.abs <- res[,2] < thresh/100

dat.abs <- cbind(x[bool.abs,1],x[bool.abs,2],res[bool.abs,], worldclim2[bool.abs,])
#if(sp.name=="ks") dat.abs <- dat.abs[1:450000,]
colnames(dat.abs)[1:4] <- c("xcoord", "ycoord","Vote","Prob")
if(file.out) write.table(dat.abs, paste(sp.name,"_data_absences_threshold_",thresh,".txt",sep=""), row.names=F)


#----------------- ------------------------ #
# Random undersampling of relevant absences
# ----------------------------------------- #

	n.abs <- round(nrow(species)*(100-preva)/preva)
	index <- sample(nrow(dat.abs), n.abs, replace=F)
	dat.abs.sample <- dat.abs[index,]
	plot(border.no.nz[,2], border.no.nz[,3], cex=0.2, xlab="Longitude", ylab="Latitude")
	points(dat.abs.sample[,1], dat.abs.sample[,2], col="blue",pch=20, cex=0.5)
	if(save.fig) savePlot(paste(sp.name,"_World_random_absences_threshold_",thresh,"_prev_",preva,".png",sep=""), type="png")
	if(file.out) write.table(dat.abs.sample, paste(sp.name,"_random_abs_threshold_",thresh,"_prev_",preva,".txt",sep=""), row.names=F)
	graphics.off()

cat("All plots done \n")


# ---------------------------------------------------- #
# K-means based undersampling, create presence/absence dataset
# ---------------------------------------------------- #

if (do.kmeans) {
	n.abs <- round(nrow(species)*(100-preva)/preva)
	data.abs.kmeans <- kmeans(dat.abs[,-(1:4)], n.abs, iter.max=100)
	pres <- cbind(species[,chosen.v], rep(1,nrow(species)))
	abs <- cbind(data.abs.kmeans$centers, rep(0,n.abs))
	colnames(pres)[ncol(pres)] <- colnames(abs)[ncol(abs)] <- "pres"
	data.pres.abs <- rbind(pres,abs)
	if(file.out) {
		write.table(data.pres.abs, paste(sp.name,"_data_pres_abs_threshold_",thresh,"_prev_",preva,".txt",sep=""), row.names=F)
		save(data.abs.kmeans, file=paste(sp.name,"_kmeans_undersampling_object_threshold_",thresh,"_prev_",preva,".RData",sep=""))
	}

	cat("K-means done \n")


}

}
