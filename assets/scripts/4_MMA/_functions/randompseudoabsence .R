randompseudoabsence <- function(sp.name, nb=10, save.fig, file.out, thresh=0,preva=50,random.nb) {


# ------------------------ #--# ------------------------ #
# ------------------------ #--# ------------------------ #
# ------------------ Do ALL of OCSVMk ------------------ #
#   gather votes, plot,get best absences, make dataset   #
# ------------------------ #--# ------------------------ #
# ------------------------ #--# ------------------------ #

# INPUT:
# sp.name = character of species name
# color.choice = color of plot, "brown","heat","rainbow","venette"
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
	setwd(paste("C:\\MMA\\",species.name,"\\_OCSVM_results\\",sp.name,sep=""))
	assign(paste("vote",i,sep=""),read.table("votes.txt",header=T))
	assign(paste("model.in.",i,sep=""),read.table("model_in.txt",header=T))
	assign(paste("model.out.",i,sep=""),read.table("model_out.txt",header=T))
	cat(paste("dimensions vote",i,sep=""),":\t",dim(get(paste("vote",i,sep=""))))
	cat("\n")
	}

# --------- #
# Sum votes # 
# --------- #

setwd(paste("C:\\MMA\\",species.name,"\\_OCSVM_results\\",sp.name,sep=""))
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

legend.loc <- c(-120,-10)
cex1 <- 0.2
pch1 <- 16
	
win.graph(width=70,height=40)


# ---------------------------------- #
# Make database of Relevant Absences
# ---------------------------------- #

if (thresh == 0) bool.abs <- res[,2] == 0
if (thresh > 0)  bool.abs <- res[,2] < thresh/100

dat.abs <- cbind(x[bool.abs,1],x[bool.abs,2],res[bool.abs,], worldclim2[bool.abs,])
#if(sp.name=="ks") dat.abs <- dat.abs[1:450000,]
colnames(dat.abs)[1:4] <- c("xcoord", "ycoord","Vote","Prob")
#if(file.out) write.table(dat.abs, paste(sp.name,"_data_absences_threshold_",thresh,".txt",sep=""), row.names=F)


#----------------- ------------------------ #
# Random undersampling of relevant absences
# ----------------------------------------- #

	n.abs <- nrow(species)*(100-preva)/preva
	index <- sample(nrow(dat.abs), n.abs, replace=F)
	dat.abs.sample <- dat.abs[index,]
	plot(border.no.nz[,2], border.no.nz[,3], cex=0.2, xlab="Longitude", ylab="Latitude")
	points(dat.abs.sample[,1], dat.abs.sample[,2], col="blue",pch=20, cex=0.5)
	if(save.fig) savePlot(paste(sp.name,"_World_random_absences_threshold_",thresh,"_prev_",preva,".png",sep=""), type="png")
	if(file.out) write.table(dat.abs.sample, paste(sp.name,"_random_abs_set",random.nb,"_threshold_",thresh,"_prev_",preva,".txt",sep=""), row.names=F)
	graphics.off()

}
