plot.uncertainty <- function(sp.name,method.name,save.fig,just.pred=F,data.out=F) {

# ------ #
# Output
# Plots map of prediction, weights, weighted prediction

# Input
# sp.name = species name
# method.name = method's name
# save.fig = logical
# just.pred = logical, if you want to just do predictions
# data.out = logical, to output weights data in a file (x,y,weight)

# ------ #

# ----------------- #
# DATA of prediction
# ----------------- #

wd3 <- paste("C:\\MMA\\",species.name,"\\_prediction_maps\\_predict_NZ",sep="")
setwd(wd3)
load(paste("NZ_predictions_",sp.name,"_",method.name,"_",p.lab,".Rdata",sep=""))
res <- cbind(train.and.test.model.predict$pred.test,train.and.test.model.predict$prob.test)
colnames(res) <- c("Prediction","Probability")

plot.prediction.color(res)
title(paste("Prediction for New Zealand",sp.name,method.name,sep=", "))

if (save.fig) {
	setwd(paste("C:\\MMA\\",species.name,"\\_prediction_maps\\_uncertainty_NZ",sep=""))
	savePlot(paste("NZ_Prediction_",sp.name,"_",method.name,"_",p.lab,".png",sep=""), type="png")
	}

# ----------------- #
# Calculate weights
# ----------------- #
if (!just.pred) {

# Change to temporary output directory
wd <- paste("C:\\MMA\\",species.name,"\\_prediction_maps\\_uncertainty_NZ",sep="")
setwd(wd)

# Create boolean that gives uncertain predictions
world.bool2 <- (species.read.model.test$mean+species.read.model.test$sd)>0.5 & (species.read.model.test$mean-species.read.model.test$sd)<0.5

# Only do weights if there is uncertainty
if (sum(world.bool2)>0) {
	weight <- rep(1,dim(worldclim2)[1])
	world2 <- species.read.model.test[world.bool2]
	world.bool3 <- (world2$mean+world2$sd)>0.5 & (world2$mean-world2$sd) <0.5
	weight.f <- function(x){abs(x-0.5)}
	weight <- apply(as.matrix(world2$mean),1,weight.f)/world2$sd
	weight <- ifelse(weight>1,1,weight) 

# Plot uncertainty map
	res2 <- cbind(train.and.test.model.predict$pred.test,weight)
	plot.prediction.color(res2)
	title(paste("Weights",sp.name,method.name,sep=", "))

	if (save.fig) {
		setwd(paste("C:\\MMA\\",species.name,"\\_prediction_maps\\_uncertainty_NZ",sep=""))
		savePlot(paste("NZ_Weights_",sp.name,"_",method.name,"_",p.lab,"_nrep",n.rep,".png",sep=""), type="png")
		}


# Plot prediction with weight
	res3 <- cbind(train.and.test.model.predict$pred.test,train.and.test.model.predict$prob.test*weight)
	plot.prediction.color(res3)
	title(paste("Weighted Prediction",sp.name,method.name,sep=", "))

	if (save.fig) {
		setwd(paste("C:\\MMA\\",species.name,"\\_prediction_maps\\_uncertainty_NZ",sep=""))
		savePlot(paste("NZ_Weighted_Prediction_",sp.name,"_",method.name,"_",p.lab,"_nrep",n.rep,".png",sep=""), type="png")
		}

	if (data.out) {
		setwd(paste("C:\\MMA\\",species.name,"\\_prediction_maps\\_uncertainty_NZ",sep=""))
		dat <- cbind(x[,1],x[,2],weight)
		colnames(dat) <- c("xcoord","ycoord","weight")
		write.table(dat,paste("NZ_Weights_",sp.name,"_",method.name,"_",p.lab,"_nrep",n.rep,".txt",sep=""),row.names=F)
		}
	
	}

# if there is no uncertainty
	else #print("no uncertainty") 
		write.table("Hooray! No uncertainty!",paste("NZ_Prediction_",sp.name,"_",method.name,"_",p.lab,"no_uncertainty.txt",sep=""),
			row.names=F,col.names=F)


}

}