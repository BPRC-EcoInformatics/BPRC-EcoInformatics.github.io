plot.split <- function(sp.name,method.name,color.choice="venette",save.fig=F,exp=F,xlim=NULL,ylim=NULL,pdf=F) {

#######################################
#######################################
# MAKE SURE YOU LOAD WORLDCLIM FIRST !!
#######################################
#######################################

# res = 2 column matrix of Prediction and Probability from model result
# x = 2 column matrix with coordinates of predictions to be colored into map
# bool = boolean for edges
# legend.loc = length 2 vector with coordinates for legend
# color.choice = character, "nuria", "rainbow", "heat", "maxent"
# save.fig = logical
# exp = F, plotted at 0.1 intervals. If T, higher end is finer
# xlim, ylim = to fit plotted map

# Splits
	spl <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
	if (exp) spl <- c(0.1,0.3,0.5,0.7,0.9,0.925,0.95,0.975,0.99)
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


	if (color.choice == "maxent") {
		spl <- c(0,0.08,0.15,0.23,0.31,0.38,0.46,0.54,0.62,0.69,0.77,0.85,0.92)
#		if (exp) spl <- c(0.1,0.3,0.5,0.7,0.9,0.925,0.95,0.975,0.99)
		ind1 <- res[,2]==spl[1]
		ind2 <- res[,2]>spl[1] & res[,2]<=spl[2]
		ind3 <- res[,2]>spl[2] & res[,2]<=spl[3]
		ind4 <- res[,2]>spl[3] & res[,2]<=spl[4]
		ind5 <- res[,2]>spl[4] & res[,2]<=spl[5]
		ind6 <- res[,2]>spl[5] & res[,2]<=spl[6]
		ind7 <- res[,2]>spl[6] & res[,2]<=spl[7]
		ind8 <- res[,2]>spl[7] & res[,2]<=spl[8]
		ind9 <- res[,2]>spl[8] & res[,2]<=spl[9]
		ind10 <- res[,2]>spl[9] & res[,2]<=spl[10]
		ind11 <- res[,2]>spl[10] & res[,2]<=spl[11]
		ind12 <- res[,2]>spl[11] & res[,2]<=spl[12]
		ind13 <- res[,2]>spl[12] & res[,2]<=spl[13]
		ind14 <- res[,2]>spl[13]

		col <- rep("white", dim(res)[1])

		col[ind1] <- rgb(0, 0, 255, maxColorValue=255) # <0.1
		col[ind2] <- rgb(26,66,255, maxColorValue=255)
		col[ind3] <- rgb(15,237,255, maxColorValue=255)
		col[ind4] <- rgb(15,237,255, maxColorValue=255)
		col[ind5] <- rgb(17,253,200, maxColorValue=255)
		col[ind6] <- rgb(50,252,121, maxColorValue=255)
		col[ind7] <- rgb(0,255,0, maxColorValue=255)
		col[ind8] <- rgb(0,255,0, maxColorValue=255)
		col[ind9] <- rgb(128,255,0, maxColorValue=255)
		col[ind10] <- rgb(200,255,18, maxColorValue=255)
		col[ind11] <- rgb(255,236,0, maxColorValue=255)
		col[ind12] <- rgb(255,152,13, maxColorValue=255)
		col[ind13] <- rgb(255,80,0, maxColorValue=255)
		col[ind14] <- rgb(255,0,0, maxColorValue=255)

		col.vec <- c(rgb(0, 0, 255, maxColorValue=255), 
				rgb(26,66,255, maxColorValue=255),
				rgb(15,237,255, maxColorValue=255),
				rgb(15,237,255, maxColorValue=255),
				rgb(17,253,200, maxColorValue=255),
				rgb(50,252,121, maxColorValue=255),
				rgb(0,255,0, maxColorValue=255),
				rgb(0,255,0, maxColorValue=255),
				rgb(128,255,0, maxColorValue=255),
				rgb(200,255,18, maxColorValue=255),
				rgb(255,236,0, maxColorValue=255),
				rgb(255,152,13, maxColorValue=255),
				rgb(255,80,0, maxColorValue=255),
				rgb(255,0,0, maxColorValue=255))


	}

	if (color.choice == "brown") {
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

	legend.loc <- c(-150,20)
	cex1 <- 0.2
	pch1 <- 16
	
	if (is.null(xlim)) X11(width=15,height=10)
	if (!is.null(xlim)) { 
#		win.graph() 
		cex1 <- 0.8
		pch1 <- 15
	}

	plot(x[,1],x[,2],cex=cex1,pch=pch1,col="lightcyan1", xlab="Longitude", ylab="Latitude",xlim=xlim, ylim=ylim)
	points(x[,1],x[,2],cex=cex1,pch=pch1,col=col)
# if plotting with Core2duo:
#	legend(legend.loc[1]-10,legend.loc[2]+10, legend=c(paste("<",spl[1]), 
if (color.choice!="maxent") {
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
}

if (color.choice=="maxent") {
	legend(legend.loc[1],legend.loc[2], legend=c(paste("<",spl[1]), 
					paste(paste(spl[1],"-"),spl[2]),
					paste(paste(spl[2],"-"),spl[3]),
					paste(paste(spl[3],"-"),spl[4]),
					paste(paste(spl[4],"-"),spl[5]),
					paste(paste(spl[5],"-"),spl[6]),
					paste(paste(spl[6],"-"),spl[7]),
					paste(paste(spl[7],"-"),spl[8]),
					paste(paste(spl[8],"-"),spl[9]),
					paste(paste(spl[9],"-"),spl[10]),
					paste(paste(spl[10],"-"),spl[11]),
					paste(paste(spl[11],"-"),spl[12]),
					paste(paste(spl[12],"-"),spl[13]),
					paste(">",spl[13])), 
					col=col.vec,
					pch=15,bg="white", cex=0.8)
}


#	title(paste("Prediction for World",sp.name,method.name,sep=", "))
	if(save.fig) savePlot(paste("World_prediction_",sp.name,"_",method.name,"_",p.lab,"_",bg,nb.var,"_",reduction,".png",sep=""), type="png")


}
