randompseudoabsence <- function(sp.name, thresh=0,preva=50,random.nb,file.out=T,save.fig=T) {


# INPUT:
# sp.name = character of species name
# thresh = 0, absences combined from predictions with probability = 0 (default).
#             if thresh = t (100>t>0), it will take marginal probabilities (percentage) < t
# preva = 50 (default), from 0 to 100
# random.nb = a random number to name output

# OUTPUT:
# a) maps
# - random absences with presences
# b) files
# - presence/absence dataset
# ------------------------ #--# ------------------------ #

# DATA
setwd(paste("C:\\MMA\\",species.name,"\\_dist_data",sep=""))
species <- read.table(paste(sp.name,"_worldclim.txt",sep=""),header=T)

ind <- species[,2:3]
ind[,3] <- 1
colnames(ind)[3] <- "pres"
species <- species[which(ind[,3]!=0),]

# ---------------------------------- #
# read in relevant absences from ocsvm (ALREADY DONE)
# ---------------------------------- #
#	setwd("R:\\Intelligent systems for biosecurity\\INVASIVE_SPP\\_dist_data_OCSVM")
#	all.absence.data.ocsvm <- read.table(paste(sp.name,"_data_absences_threshold_",thresh,".txt",sep=""))


#------------------------------------------ #
# Random undersampling of relevant absences
# ----------------------------------------- #

	setwd(paste("C:\\MMA\\",species.name,"\\_dist_data_OCSVM",sep=""))
	n.abs <- nrow(species)*(100-preva)/preva
	index <- sample(nrow(all.absence.data.ocsvm), n.abs, replace=F)
	dat.abs.sample <- all.absence.data.ocsvm[index,]
#	if(file.out) write.table(dat.abs.sample, paste(sp.name,"_random_abs_set",random.nb,"_threshold_",thresh,"_prev_",preva,".txt",sep=""), row.names=F)

#	win.graph(width=70,height=40)
#	plot(border.no.nz[,2], border.no.nz[,3], cex=0.2, xlab="Longitude", ylab="Latitude")
#	points(as.numeric(as.character(dat.abs.sample[,1])), as.numeric(as.character(dat.abs.sample[,2])), col="blue",pch=20, cex=0.5)
#	if(save.fig) savePlot(paste(sp.name,"_World_random_abs_set",random.nb,"_threshold_",thresh,"_prev_",preva,".png",sep=""), type="png")
#	points(species[,2], species[,3], col="red",pch=20, cex=0.5)
#	legend(-150,-30,col=c("red","blue"),pch=20,legend=c("Presence","Absence"))
#	if(save.fig) savePlot(paste(sp.name,"_World_pres_random_abs_set",random.nb,"_threshold_",thresh,"_prev_",preva,".png",sep=""), type="png")
#	graphics.off()

#------------------------------------------ #
# Output pres and abs dataset
# ----------------------------------------- #

pres <- cbind(species[,chosen.v], rep(1,nrow(species)))
abs <- cbind(dat.abs.sample, rep(0,n.abs))
abs <- abs[,-(1:4)]
colnames(pres)[ncol(pres)] <- colnames(abs)[ncol(abs)] <- "pres"
colnames(abs) <- colnames(pres)

dat.pres.abs <- rbind(pres,abs)
if(file.out) write.table(dat.pres.abs, paste(sp.name,"_data_pres_random_abs_set",random.nb,"_threshold_",thresh,"_prev_",preva,".txt",sep=""), row.names=F)
	

}
