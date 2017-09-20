plot.prediction.color <- function(result) {

# does the splitting and plotting

	spl <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
	ind1 <- result[,2]<spl[1]
	ind2 <- result[,2]>spl[1] & result[,2]<=spl[2]
	ind3 <- result[,2]>spl[2] & result[,2]<=spl[3]
	ind4 <- result[,2]>spl[3] & result[,2]<=spl[4]
	ind5 <- result[,2]>spl[4] & result[,2]<=spl[5]
	ind6 <- result[,2]>spl[5] & result[,2]<=spl[6]
	ind7 <- result[,2]>spl[6] & result[,2]<=spl[7]
	ind8 <- result[,2]>spl[7] & result[,2]<=spl[8]
	ind9 <- result[,2]>spl[8] & result[,2]<=spl[9]
	ind10 <- result[,2]>spl[9]

col <- rep("grey", dim(result)[1])
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

	legend.loc <- c(168,-35)
	cex1 <- 0.4
	pch1 <- 16

plot(x.nz[,1],x.nz[,2],cex=cex1,pch=pch1,col="lightcyan1", xlab="Longitude", ylab="Latitude")
points(x[,1],x[,2],cex=0.2,pch=16,col=col, xlab="Longitude", ylab="Latitude")
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