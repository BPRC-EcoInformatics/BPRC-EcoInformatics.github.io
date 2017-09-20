plot.split.ocsvm <- function(species.name,color.choice="nuria",save.fig=F) {

#######################################
#######################################
# MAKE SURE YOU LOAD WORLDCLIM FIRST !!
#######################################
#######################################

# res = 2 column matrix of Prediction and Probability from model result
# x = 2 column matrix with coordinates of predictions to be colored into map
# bool = boolean for edges
# legend.loc = length 2 vector with coordinates for legend
# color.choice = character, "nuria", "rainbow", "heat"
# save.fig = logical

# Load data
#setwd("R:\\Intelligent systems for biosecurity\\INVASIVE_SPP\\_dist_data&gbif")
#data1 <- read.table(paste(species.name,"_data_pres_abs.txt",sep=""), header=TRUE)
#setwd(paste("R:\\Intelligent systems for biosecurity\\INVASIVE_SPP\\_OCSVM\\",species.name,sep=""))
#res <- cbind(final.vote, final.vote/100)
#colnames(res) <- c("Vote","%")


## FILTER
#f <- function(x) sum(is.na(x))==0
#bool <- apply(worldclim[,colnames(data1)[1:(ncol(data1)-1)]],1,f)
#worldclim2 <- worldclim[bool,colnames(data1)[1:(ncol(data1)-1)]]
#bool2 <- apply(worldclim.orig[,c(59,60)],1,f)
#x <- worldclim[bool,c(2,3)]


# Splits
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

	col <- rep("grey", dim(res)[1])
	if (color.choice == "nuria") {
		col[ind1] <- rgb(207, 207, 207, maxColorValue=255)
		col[ind2] <- rgb(255, 215, 0, maxColorValue=255)
		col[ind3] <- rgb(255, 193, 37, maxColorValue=255)
		col[ind4] <- rgb(255, 180, 34, maxColorValue=255)
		col[ind5] <- rgb(255, 185, 15, maxColorValue=255)
		col[ind6] <- rgb(255, 165, 0, maxColorValue=255)
		col[ind7] <- rgb(205, 104, 57, maxColorValue=255)
		col[ind8] <- rgb(139, 76, 57, maxColorValue=255)
		col[ind9] <- rgb(139, 54, 38, maxColorValue=255)
		col[ind10] <- rgb(139, 37, 0, maxColorValue=255)
		col.vec <- c(rgb(250, 250, 210, maxColorValue=255), 
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
		col[ind1] <- rainbow(12)[10]
		col[ind2] <- rainbow(12)[9]
		col[ind3] <- rainbow(12)[8]
		col[ind4] <- rainbow(12)[7]
		col[ind5] <- rainbow(12)[6]
		col[ind6] <- rainbow(12)[5]
		col[ind7] <- rainbow(12)[4]
		col[ind8] <- rainbow(12)[3]
		col[ind9] <- rainbow(12)[2]
		col[ind10] <- rainbow(12)[1]
		col.vec <- rev(rainbow(12)[1:10])
	}

	if (color.choice == "heat") {
		col[ind1] <- heat.colors(10)[10]
		col[ind2] <- heat.colors(10)[9]
		col[ind3] <- heat.colors(10)[8]
		col[ind4] <- heat.colors(10)[7]
		col[ind5] <- heat.colors(10)[6]
		col[ind6] <- heat.colors(10)[5]
		col[ind7] <- heat.colors(10)[4]
		col[ind8] <- heat.colors(10)[3]
		col[ind9] <- heat.colors(10)[2]
		col[ind10] <- heat.colors(10)[1]
		col.vec <- rev(heat.colors(10))
	}
	legend.loc <- c(-120,-10)
	cex1 <- 0.2
	pch1 <- 16
	
	win.graph(width=70,height=40)
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
					pch=15,bg="white", cex=0.8)


	if(plot.fig == "world") {
		points(border[,2], border[,3], cex=0.2, xlab="Longitude", ylab="Latitude")
		title(paste("Prediction for World",species.name,method.name,sep=","))
		if(save.fig) savePlot(paste(species.name,"_plotWorld_prediction_",method.name,".png",sep=""), type="png")
	}
	if(plot.fig == "nz") {
		points(border[,2], border[,3], cex=0.2, xlab="Longitude", ylab="Latitude")
#		points(worldclim.orig[!bool2,2],worldclim.orig[!bool2,3],cex=1,pch=15,col="grey")
		title(paste("Prediction for NZ",species.name,method.name,sep=","))
		if(save.fig) savePlot(paste(species.name,"_plotNZ_prediction_",method.name,".png",sep=""), type="png")
	}




}
