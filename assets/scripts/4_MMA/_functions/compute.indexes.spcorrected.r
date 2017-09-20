
compute.indexes.spcorrected <- function(prob.test,test_coord,coord.spind_order){
##### Predictions

# Extract predictions
predictions.spind <- prob.test
# Extract actuals
actuals.spind <- test_coord$pred
# Create data frame with actuals and predicted values
data.spind <- as.data.frame(cbind(actuals.spind,prob.test))
data.spind[,1] <- as.numeric(data.spind[,1])
data.spind[,2] <- as.numeric(data.spind[,2])
#cbind(c(actuals.spind[actuals.spind==1],actuals.spind[actuals.spind==0]),c(prob.test[actuals.spind==1],prob.test[actuals.spind==0]))

# Extract coordinates from presence/absence points in the training dataset
##### For plotting to check if the conversion to x-y coordinates as integer worked!
#x.world <- colFromX(raster_border,data.frame(rasterToPoints(raster_border))[,1])
#y.world <- rowFromY(raster_border,data.frame(rasterToPoints(raster_border))[,2])

#plot(range(2250,0), range(1125,0), type="n", ylim= rev(range(1125,0)) ) 
#abline(v=(seq(0,2250,1)), col="lightgray", lty="dotted", ylim= rev(range(1125,0)))
#abline(h=(seq(0,1125,1)), col="lightgray", lty="dotted", ylim= rev(range(1125,0)))
#points(x.world,y.world,col='black', cex=0.75,ylim= rev(range(1125,0)))
#points(coord.spind[1:(length(actuals.spind[actuals.spind==1])),1],coord.spind[1:(length(actuals.spind[actuals.spind==1])),2],col='red',pch=18, cex=0.75,ylim= rev(range(1125,0)))
#points(coord.spind[(length(actuals.spind[actuals.spind==1])+1):(dim(coord.spind)[1]),1],coord.spind[(length(actuals.spind[actuals.spind==1])+1):(dim(coord.spind)[1]),2],col='forestgreen',pch=18, cex=0.75,ylim= rev(range(1125,0)))


## Calculation of sptially corrected indexes
#coord.spind <- coord.spind_order[,c(1,2)]
coord.spind <- as.data.frame(coord.spind_order[,c(1,2)])
#coord.spind[,1] <- as.numeric(coord.spind[,1])
#coord.spind[,2] <- as.numeric(coord.spind[,2])

print('Print dimension - compute.indexes.spcorrected')
print(dim(coord.spind))
print(length(data.spind[,2]))
#print(class(coord.spind))
#print(coord.spind[,1])
print(coord.spind[,2])
print(as.numeric(data.spind[,2]))
#print(nrow(coord.spind))
#print(ncol(coord.spind))
#print(dim(coord.spind))
#print(head(data.spind))


# spatial autocorrelation
ac<- NA
si1<- NULL
si2<- NULL


if(length(unique(data.spind[,2]))!=1){
print('Calculating spatial autocorrelation')
ac<-acfft(coord.spind,as.numeric(data.spind[,2]))
print('Index threshold dependent')
si1 <-th.dep(data.spind,coord.spind,spatial=TRUE)
print('Index threshold independent')
si2 <-th.indep(data.spind,coord.spind,spatial=TRUE)
}

out <- t(as.matrix(c(ac,si1$sensitivity,si1$kappa,si1$specificity,si2$TSS,si2$AUC),1,6))
colnames(out) <- c("sp_autocorrelation","Recall","Kappa","Specificity","TSS","AUC")
rownames(out) <- "Value"
return(out)

}
