#  Author: Tak Ikeda
# Date: 10 July 2013
# Modified by Audrey Lustig - January 2017
# Sensitivity and robustness analysis


##################################################################################################################################   
###################################################     Create variables and directories 
# Name species, universal throughout folder names and output
#start the timer
# Set main directories
root_wd <- "C:/Users/Ursula Torres/Documents/Ecoinformatics_workshop"
main_wd <- "C:/Users/Ursula Torres/Documents/Ecoinformatics_workshop/4_MMA" 
data_wd <- "C:/Users/Ursula Torres/Documents/Ecoinformatics_workshop/1_CleaningData"
setwd(main_wd)

species.name='Ccapitata'
s
##################################################################################################################################  
###################################################     Load packages and functions

# load packages for mapping, statistical analysis
require(sp)
require(kernlab)
require(MASS)
require(class)
require(e1071)
require(randomForest)
require(varSelRF)
require(ROCR)
require(rpart)
require(party)
require(nnet)
require(grDevices)
require(klaR)
require(raster)
require(spind)
require(maptools) 

# load MMA functions
func.dir <- paste(main_wd,"_functions",sep="/")
sourceDir <- function(path, trace = TRUE, ...) {
    for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
       if(trace) cat(nm," ")           
       source(file.path(path, nm), ...)
       if(trace) cat('\n')
    }
}
sourceDir(func.dir)

##################################################################################################################################  
###################################################     Load species occurences
# Load the table of latitue and longitude coordinates
all_data <- read.table(paste(data_wd,'/',species.name,'/',species.name,'_clim_range_withoutRep.txt',sep=''),header=T,sep='\t')
all_data <- all_data[c("xcoord","ycoord", names(all_data)[1:35])] 



##################################################################################################################################  
###################################################     Load world climate data
# Choose the appropriate dataset
# P2DR1 - consists of 35 variables derived from temperature, precipitation, radiation and water-balance soil moisture variables of the CLIMOND dataset 
rdata.name <- "P2DR1" 

# Load worldclim dataset: table of latitutde and longitude and associated worldclim data that define regions of the world
load(paste(main_wd,"/_database/",rdata.name,".Rdata",sep=""))

# Transform dataframe to spatialpoints
rdata.coord<-get(rdata.name)
coordinates(rdata.coord) <- c("xcoord", "ycoord")
climate.points <- SpatialPoints(rdata.coord)
proj4string(climate.points) <- "+proj=longlat +datum=WGS84"

# Extract world administrative country border

border <-readShapePoly(paste(main_wd,"/_database/_data_border/gadm28_adm0.shp",sep=""))

proj4string(border) <- "+proj=longlat +datum=WGS84" 





##################################################################################################################################  
###################################################     Read presence/absence data and associated climond varibales.
# Create directory for named after species.name
setwd(paste(main_wd, sep="/"))
dir.create(species.name,showWarnings=F)
setwd(species.name)

print('Create output folders')
# Create 9 folders to order outputs
dir.create("_condor",showWarnings=F) 
dir.create("_dist_data",showWarnings=F)
dir.create("_dist_data_OCSVM",showWarnings=F)
dir.create("_OCSVM_results",showWarnings=F)
dir.create("_results_compar.models",showWarnings=F)
dir.create("_prediction_maps",showWarnings=F)
dir.create("_data",showWarnings=F)
setwd(paste(main_wd,species.name,"/_prediction_maps",sep="/"))
#dir.create("_predict_NZ",showWarnings=F)
dir.create("_predict_world",showWarnings=F)
setwd(paste(main_wd,species.name,sep="/"))




##################################################################################################################################  
###################################################     Choose background data
## Now you have your location points, climate data in these locations and ouput folders defined. But you also need some points that define regions where the studied species aren’t found. 
# One option is to define a background? region to sample either at random or using the three step analysis, which (we hope) captures environmental conditions that the studied species could disperse to, but haven't.

# choose background (random or 3step)
bg<-"3step"

#for(bg in c('random','3step')){
print(paste('set background:',species.name, '-',bg,sep=' '))

# Choose number of variables 
nb.var=35
print(paste('set number of variables:', species.name, '-',bg,'-',nb.var,sep=' '))


if(bg =='3step') {
	# Add the background files from 3 step analysis
	bgfile1<-paste(root_wd,"/3_3stepanalysis_firststep/",species.name,"/",species.name,"_opt_dist.csv",sep="")
	bndpool=read.csv(file =bgfile1, header =T,sep=" ") #make sure column names match with worldclim data
	colnames(bndpool) <- colnames(P2DR1)
	threestep_bg <- as.data.frame(bndpool)
	save(threestep_bg,file=paste(main_wd,'/',species.name,"/_data/",species.name,"_3stepcustom_bg.Rdata",sep=""))


	custom_bg <- threestep_bg
	dir.create(paste(main_wd,species.name,"_OCSVM_results",bg,sep="/"))
	dir.create(paste(main_wd,species.name,"_OCSVM_results",bg,nb.var,sep="/"))
	dir.create(paste(main_wd,species.name,"_dist_data",bg,sep="/"))
	dir.create(paste(main_wd,species.name,"_dist_data",bg,nb.var,sep="/"))
	
	# Check if BG bounds presence points 
	setwd(paste(main_wd,species.name,"_dist_data",bg,nb.var,sep="/")) 
	spdata.check <-all_data
	require(maps)	
	X11(width=15,height=10)
	map('world', fill = TRUE, col = "white")
	points(custom_bg$xcoord,custom_bg$ycoord,col='blue',pch=18, cex=0.5)
	points(spdata.check$xcoord,spdata.check$ycoord,col='red',pch=18, cex=0.5)
	legend(-150,-40,col=c('blue','red'),legend=c("Background","Presences"),pch=16)
	savePlot(paste(species.name,"_presence_background_check_",bg,nb.var,".png",sep=""), type="png")
	cat("background check done \n")
	dim(custom_bg)
	graphics.off()

	# assign BG data for pseudo-absence selection 
	worldclim <- custom_bg
	chosen.v <- c(3:(nb.var+2))
	custom_bg <- custom_bg[,c(1,2,chosen.v)]
	worldclim2 <- custom_bg


	###################################################################################################################
	# Environmental profiling
	setwd(paste(main_wd,"/_database/_spaces",sep=""))
	for (i in 1:10) 
  		assign(paste("space",i,sep=""),
	   		read.table(paste("OneSVM_space_parameter_",i,".txt",sep=""), header=T))

	#OCSVM
	print(paste('Environmental profiling: ', species.name, '-',bg,'-',nb.var,sep=' '))
	setwd(paste(main_wd,species.name,"_dist_data",bg,nb.var,sep="/")) 
	species <- all_data#[,c(1,2,chosen.v)]
	names(species)<-names(custom_bg)
	# output OCSVM files 
	setwd(paste(main_wd,species.name,"_OCSVM_results",bg,nb.var,sep="/"))

	cv.k <- 10
	if (nrow(species)<100 & nrow(species)>=10) cv.k <- 5
	if (nrow(species)<10) {
		cv.k <- NA
		write.table("dataset less than 10 rows","error.txt",row.names=F,col.names=F)
	}

	test <- OneSVM.param.cross.clim(data1.x=species[,chosen.v],
						data0.x=worldclim2[,chosen.v], 
						nb.models=10, k=cv.k,  			###############change to 100			
						prior=rbind(space1,space2,space3,space4,space5,
						space6,space7,space8,space9,space10),
						file=T, save.out=T,bg=bg)



	# Read in and join OCSVM results 
	setwd(paste(main_wd,species.name,"_OCSVM_results",bg,nb.var,sep="/"))
	load(paste("one_class_svm_out_",bg,nb.var,".RData",sep=''))
	outOneSVM<-OneSVM.param.cross.clim.read(out$models.in, path=NULL, file=T, bg=bg)

	# OCSVMk function 
	th <- 0
	pr <- 50
	ocsvmk(species.name,color.choice="venette",save.fig=T,file.out=T,do.kmeans=T, thresh= th,prev=pr,bg=bg,nb.var=nb.var) 

	# Plot representative/relevant absences from kmeans clusters 
	kmeans.plot(species.name, bg, thresh= th, save.fig=T, file.out=T,prev=pr,nb.var)

	# Create a file of pres/abs with their respective coordinates
	setwd(paste(main_wd,species.name,"_dist_data_OCSVM",sep="/"))
	abs_coord<-read.table(paste(species.name,"_kmeans_abs_thresh_",th,"_prev_",pr,"_",bg,nb.var,".txt",sep=""),h=T)
	abs_coord<-cbind(abs_coord[,-(3:4)],rep(0,nrow(abs_coord)))
	pres_coord<-cbind(species[,c(1:2,chosen.v)], rep(1,nrow(species)))
	colnames(pres_coord)[ncol(pres_coord)] <- colnames(abs_coord)[ncol(abs_coord)] <- "pres"
	data.pres.abs_coord <- rbind(pres_coord,abs_coord)

	graphics.off()
	write.table(data.pres.abs_coord, paste(species.name,"_data_pres_abs_thresh_",th,"_prev_",pr,"_coord_",bg,nb.var,".txt",sep=""), row.names=F)
}

if(bg =='random') {
	#custom_bg <- threestep_bg
	dir.create(paste(main_wd,species.name,"_OCSVM_results",bg,sep="/"))
	dir.create(paste(main_wd,species.name,"_OCSVM_results",bg,nb.var,sep="/"))
	dir.create(paste(main_wd,species.name,"_dist_data",bg,sep="/"))
	dir.create(paste(main_wd,species.name,"_dist_data",bg,nb.var,sep="/"))

	require(dismo)
	require(raster)

	# create spatial point in the object all_data_spp
	all_data_spp <- all_data
	coordinates(all_data_spp) <- ~xcoord+ycoord
	projection(all_data_spp) <- CRS('+proj=longlat +datum=WGS84')

	# circles with a radius of 50 km
	xcircle <- circles(all_data_spp, d=50000, lonlat=TRUE)
	# transform to polygon
	pol<-xcircle@polygons
	# select only buffer area that overlap with terrestrial area
	spborder<-as(border,"SpatialPolygons")
	cpol_border<-raster::intersect(pol,spborder)
	projection(cpol_border) <- projection(all_data_spp)
	# sample randomly from all buffer/circles (random location all circle)
	bgcircle = spsample(cpol_border, dim(all_data_spp)[1], type='stratified')
	xy <- as.data.frame(bgcircle)
	coordtest=xy
	coordinates(coordtest)=c("x1","x2") #transform numbers in spatial data (coordinates)
	ind_zerotest=zerodist(coordtest,zero=0.17,unique.ID=T) # identified redundant coordinate within a resolution of 10 deg
	ind_tokeep=which(duplicated(ind_zerotest)==F)
	xy2=as.data.frame(coordtest[ind_tokeep,])


	setwd(paste(main_wd,species.name,"_dist_data",bg,nb.var,sep="/")) 
	require(maps)		
	X11(width=15,height=10)
	map('world', fill = TRUE, col = "white")
	lines(cpol_border,col='blue')
	points(all_data_spp$xcoord,all_data_spp$ycoord, cex=0.75, pch=20, col='red')
	points(xy2, cex=0.75, pch=20, col='green')
	legend(-150,-40,col=c('red','green'),legend=c("Presence","Absence"),pch=16)
	savePlot(paste(species.name,"_presence_background_check_",bg,nb.var,".png",sep=""), type="png")
	cat("background check done \n")
	graphics.off()
	
	dir<-Sys.glob(paste(root_wd,"/0_Data/Climond_CM10_1975H_Bio_ESRI_V1.2/CM10_1975H_Bio_V1.2/bio*/hdr.adf",sep=''))
	# Select the climatic variables to use
	dir<-dir[1:nb.var]

	bio<-lapply(dir,raster) # transform data into raster
	climond<-stack(bio) # concatenates multiple vectors contained in bio into a single vector along with a factor indicating where each observation 
	dist_clim<-extract(climond,xy2,df=T)
	abs_coord<- cbind(xy2,dist_clim[,2:(nb.var+1)], rep(0,nrow(xy2)))
	names(abs_coord) <- c("xcoord", "ycoord", names(rdata.coord)[1:(nb.var)],"pres")
	dim(abs_coord)

	th <- 0
	pr <- 50
	chosen.v <- c(3:(nb.var+2))
	pres_coord<-cbind(all_data[,c(1:2,chosen.v)], rep(1,nrow(all_data)))
	colnames(pres_coord)<- colnames(abs_coord)

	data.pres.abs_coord <- rbind(pres_coord,abs_coord)
	data.pres.abs_coord <- data.pres.abs_coord[complete.cases(data.pres.abs_coord),] 
	setwd(paste(main_wd,species.name,"_dist_data_OCSVM",sep="/"))
	write.table(data.pres.abs_coord, paste(species.name,"_data_pres_abs_thresh_",th,"_prev_",pr,"_coord_",bg,nb.var,".txt",sep=""), row.names=F)
	write.table(data.pres.abs_coord[,3:ncol(data.pres.abs_coord)], paste(species.name,"_data_pres_abs_thresh_",th,"_prev_",pr,"_",bg,nb.var,".txt",sep=""), row.names=F)
}

reduction<-"randomForest"
#for (reduction in c('randomForest','PCA')){
print(paste('Set reduction technique:',species.name, '-',bg,'-',nb.var,'-',reduction,sep=' '))

if(reduction == 'randomForest'){
	setwd(paste(main_wd,species.name,"_dist_data_OCSVM",sep="/"))
	dir.create(paste(main_wd,species.name,"_dist_data_OCSVM",reduction,sep="/"))
	# Random Forest
	print('variable selection - ocsvm')
	vs <- varSel.species(species.name,thresh=th, preva=pr,bg,nb.var) 
	write.table(vs$rf_initial, paste(main_wd,"/",species.name,"/_dist_data_OCSVM/",reduction,"/",species.name,"_varSel_rf_initial_thresh_",th,"_prev_",pr,"_",bg,nb.var,".txt",sep=""),row.names=F)
	write.table(vs$step_initial, paste(main_wd,"/",species.name,"/_dist_data_OCSVM/",reduction,"/",species.name,"_varSel_step_initial_thresh_",th,"_prev_",pr,"_",bg,nb.var,".txt",sep=""),row.names=F)
	c(length(vs$rf_initial), length(vs$step_initial))

	## Which Variables were selected?
	vs_best <- read.table(paste(main_wd,"/",species.name,"/_dist_data_OCSVM/",reduction,"/",species.name,"_varSel_rf_initial_thresh_",th,"_prev_",pr,"_",bg,nb.var,".txt",sep=""),header=T)
	print("number of variable selected:", length(as.character(c(vs_best)$x)))
	dat <- read.table(paste(main_wd,"/",species.name,"/_dist_data_OCSVM/",species.name,"_data_pres_abs_thresh_",th,"_prev_",pr,"_",bg,nb.var,".txt",sep=""),header=T)
	dat_best <- dat[,c(as.character(c(vs_best)$x),"pres")]
	write.table(dat_best, paste(main_wd,"/",species.name,"/_dist_data_OCSVM/",reduction,"/",species.name,"_data_pres_abs_best_",bg,nb.var,".txt",sep=""),row.names=F)

	##Subset variables for dataframe with coordinates
	dat_coord <- read.table(paste(main_wd,"/",species.name,"/_dist_data_OCSVM/",species.name,"_data_pres_abs_thresh_",th,"_prev_",pr,"_coord_",bg,nb.var,".txt",sep=""),header=T)
	dat_coord_best <- dat_coord[,c("xcoord","ycoord",as.character(c(vs_best)$x),"pres")]
	write.table(dat_coord_best, paste(main_wd,"/",species.name,"/_dist_data_OCSVM/",reduction,"/",species.name,"_data_pres_abs_best_coord_",bg,nb.var,".txt",sep=""),row.names=F)
}

if(reduction == 'PCA'){
	setwd(paste(main_wd,species.name,"_dist_data_OCSVM",sep="/"))
	dir.create(paste(main_wd,species.name,"_dist_data_OCSVM",reduction,sep="/"))
	library(ade4)
	print('variable selection - pca')
	vs.pca <- varSel.species.pca(species.name,thresh=th, preva=pr,bg,nb.var) 
	write.table(vs.pca, paste(main_wd,"/",species.name,"/_dist_data_OCSVM/",reduction,"/",species.name,"_varSel_pca_initial_thresh_",th,"_prev_",pr,"_",bg,nb.var,".txt",sep=""),row.names=F,quote=FALSE)

	## Which Variables were selected?
	vs_best <- vs.pca[,1]
	dat <- read.table(paste(main_wd,"/",species.name,"/_dist_data_OCSVM/",species.name,"_data_pres_abs_thresh_",th,"_prev_",pr,"_",bg,nb.var,".txt",sep=""),header=T)
	dat_best <- dat[,c(as.character(vs_best),"pres")]
	print(paste("number of variable selected:",  length(c(as.character(vs_best)))))
	write.table(dat_best, paste(main_wd,"/",species.name,"/_dist_data_OCSVM/",reduction,"/",species.name,"_data_pres_abs_best_",bg,nb.var,".txt",sep=""),row.names=F,quote=FALSE)

	##Subset variables for dataframe with coordinates
	dat_coord <- read.table(paste(main_wd,"/",species.name,"/_dist_data_OCSVM/",species.name,"_data_pres_abs_thresh_",th,"_prev_",pr,"_coord_",bg,nb.var,".txt",sep=""),header=T)
	dat_coord_best <- dat_coord[,c("xcoord","ycoord",as.character(vs_best),"pres")]
	write.table(dat_coord_best, paste(main_wd,"/",species.name,"/_dist_data_OCSVM/",reduction,"/",species.name,"_data_pres_abs_best_coord_",bg,nb.var,".txt",sep=""),row.names=F)
}


# ------------ #
# Import data
# ------------ #
	setwd(paste(main_wd,species.name,"_dist_data_OCSVM",sep="/"))
	print('read best variable')	
	data <- read.table(paste(main_wd,"/",species.name,"/_dist_data_OCSVM/",reduction,"/",species.name,"_data_pres_abs_best_",bg,nb.var,".txt",sep=""), header=TRUE)	
	data <- as.data.frame(data[sample(nrow(data)),])
	data.x <- scale(data[,1:(ncol(data)-1)], center=TRUE, scale=TRUE)
	data.x <- as.data.frame(data.x)
	data.y <- data[,ncol(data)]
	# Prevalence
	prevalence <- sum(data.y)

#-------------------#
#  PARAMETRIZATION  |
#-------------------#
# Set output directory
	setwd(paste(main_wd,species.name,"_results_compar.models",sep="/"))###UT
	cv.k <- 10
	if (prevalence < 50) cv.k <- 5
	knrep <- 200/cv.k
      knrep<-2
  
	###resampling method
	resampling_method<-"cv"#choose the resampling method that you used to extrac parameters

print(paste('Parametrization:',species.name, '-',bg,'-',nb.var,'-',reduction,sep=' '))

# Parametrization of the KNN by bootstrapping/CV
param.knn(x=data.x, y=data.y, knn.rep=knrep, resampling=resampling_method, kfold=cv.k, nrep=knrep)
graphics.off()
cat("done knn cv \n")

# Parametrization of the Support Vector Machine by bootstrapping/CV
thresh_SVM<-0
prev_SVM<-NA
param.svm(x=data.x, y=data.y, kernL="rbfdot", resampling=resampling_method, kfold=cv.k, nrep=knrep,thresh=thresh_SVM,prev=prev_SVM)
graphics.off()
cat("done param.svm, resample = cv, threshold =",0, "\n")

param.nnet(x=data.x , y=data.y , nnet.rep =knrep, resampling = resampling_method, nrep=knrep, file=T)

cat("done nnet cv \n")


#----------------------#
#  Parameter settings  |
#----------------------#
###knn
load(paste(species.name,"_param_knn_",resampling_method,"_output_object_",bg,nb.var,"_",reduction,".Rdata",sep=""))
out_knn<-as.data.frame(out$k.distribution)
ind_max<-which.max(out_knn)
best_k<-as.numeric(colnames(out_knn[ind_max]))
if(resampling_method=="cv")c.k1<-best_k
if(resampling_method=="boot")b.k1<-best_k

###svm
thresh_SVM<-0
prev_SVM<-NA
out_svm<-read.table(paste(species.name,"_param_svm_",resampling_method,"_best_models_thresh_",thresh_SVM,"_","prev_",prev_SVM,"_",bg,nb.var,"_",reduction,".txt",sep=""),h=T,row.names=NULL)

if(resampling_method=="cv")c.s1<-out_svm[1,"C"]
if(resampling_method=="boot")b.s1<-out_svm[1,"C"]

if(resampling_method=="cv")c.s2<-out_svm[1,"gamma"]
if(resampling_method=="boot")b.s2<-out_svm[1,"gamma"]


###nnet
out_nnet<-read.table(paste(species.name,"_param_nnet_",resampling_method,"_best_",bg,nb.var,"_",reduction,".txt",sep=""),h=T)

if(resampling_method=="cv")c.n1<-out_nnet[1,"Size"]
if(resampling_method=="boot")b.n1<-out_nnet[1,"Size"]

if(resampling_method=="cv")c.n2<-out_nnet[1,"Maxiter"]
if(resampling_method=="boot")b.n2<-out_nnet[1,"Maxiter"]

if(resampling_method=="cv")c.n3<-out_nnet[1,"Decay"]
if(resampling_method=="boot")b.n3<-out_nnet[1,"Decay"]



# ------------------------------------------------------------ #
# Import data- again to include the full set of pres-abs data
# ------------------------------------------------------------ #
	setwd(paste(main_wd,species.name,"_dist_data_OCSVM",reduction,sep="/"))
	data_coord <- read.table(paste(species.name,"_data_pres_abs_best_coord_",bg,nb.var,".txt",sep=""), header=TRUE)#include coordinates for chapter3
	
	# extract presence/absence coordinates separatly to match data.frame with new coordinate systems
	# need this step to avoid the issue of having a presence and absence in the same cell that will mess up the order function
	data_coord_pres <- data_coord[data_coord$pres==1,]
	data_coord_pres_order <- data_coord_pres[with(data_coord_pres, order(xcoord,ycoord)),]
	data_coord_abs <- data_coord[data_coord$pres==0,]	
	data_coord_abs_order <- data_coord_abs[with(data_coord_abs, order(xcoord,ycoord)),]
	data_coord_order <- rbind(data_coord_pres_order,data_coord_abs_order)

	# sample rows arbitrarely to mix the coordinate again
	sample_row <- sample(nrow(data_coord))
	data_coord_order <- as.data.frame(data_coord_order[sample_row,])


	data <- data_coord_order[,-c(1,2)] # data best without coordinates
	data.x <- scale(data[,1:(ncol(data)-1)], center=TRUE, scale=TRUE)
	data.x <- as.data.frame(data.x)
	data.y <- data[,ncol(data)]
	# Prevalence
	prevalence <- sum(data.y)

	data.x_coord <- scale(data_coord_order[,3:(ncol(data_coord_order)-1)], center=TRUE, scale=TRUE)
	data.x_coord <- as.data.frame(data.x_coord)
	data.x_coord <-cbind(data_coord_order[,1:2],data.x_coord)
	data.y_coord <- data_coord_order[,ncol(data_coord_order)]

#-------------------------#
#  Fix Weights for nnet   |
#-------------------------#
# given data.x, data.y, nnet.par
# Set output directory
print(paste('Compare results:',species.name, '-',bg,'-',nb.var,'-',reduction,sep=' '))
setwd(paste(main_wd,species.name,"_results_compar.models",sep="/"))


nnet.Wts.cv <- best.nnet.weights(species.name,nreps=10,  ######################################################## CARFULL TO PUT BACK 1000!!!
					data.x,data.y,
					nnet.par=c(c.n1,c.n2,c.n3),
					save.out=T)


#--------------#
#  COMPARISON  |
#--------------#
# CROSS-VALIDATION
setwd(paste(main_wd,species.name,"_results_compar.models",sep="/"))
print(paste('Cross validation:',species.name, '-',bg,'-',nb.var,'-',reduction,sep=' '))
#sourceDir(func.dir)
result.compar.cv <- compar.models.saved2(x=data.x_coord,
                               y=data.y_coord, 
						 methods=c("LDA","QDA","LOG","NB","CART","CTREE","KNN","SVM","NNET"),
                               knn.par = c.k1, 
                               svm.par = c("rbfdot",c.s1,c.s2),                                
                               nnet.par = c(c.n1,c.n2,c.n3),
                               Wts = nnet.Wts.cv,
                               resampling="cv",
                               nrep=10,   ######################################################## CARFULL TO PUT BACK 200!!!
						 nbtree   = 1000,                              
						 train.frac=0.75,
                               k=cv.k,
                               file = T, 
                               plots=T, 
                               plots.all=T, 
                               boxplots=T, 
                               bagging=NULL,
                               boosting=NULL)
						 #sp=TRUE,
						 #coord.spind_order=coord.spind_order)

# --------------------------------------------------- #
# Best models & ranked criteria tables saved into files
# --------------------------------------------------- #
print(paste('Best model criteria:',species.name, '-',bg,'-',nb.var,'-',reduction,sep=' '))
best.model.criteria(species.name,resampling_method)

#####################
### PREDICT World ###
#####################
setwd(paste(main_wd,species.name,"_dist_data_OCSVM",sep="/"))####UT  /

#-------- Train and test with best model -------
setwd(paste(main_wd,species.name,"_results_compar.models",sep="/"))
print(paste('Train and test:',species.name, '-',bg,'-',nb.var,'-',reduction,sep=' '))
nnet.Wts.cv <- read.table(paste(species.name,"_best_weights_",c.n1,"_",c.n2,"_",c.n3,"_",bg,nb.var,"_",reduction,".txt",sep=""),header=T)
nnet.Wts.cv <- c(nnet.Wts.cv$x)

# --------------------------------------------------------------------- #
# Plotting and saving output according to name of method and resampling --QDA cv
# --------------------------------------------------------------------- #
worldclim <- get(rdata.name)
f <- function(x) sum(is.na(x))==0
bool <- apply(worldclim[,colnames(data)[1:(ncol(data)-1)]],1,f)
worldclim2 <- worldclim[bool,colnames(data)[1:(ncol(data)-1)]]
x <- worldclim[bool,c(1,2)]

for (method.name in c("LDA","QDA","LOG","NB","CART","CTREE","KNN","SVM","NNET")){ # Capital letters!
setwd(paste(main_wd,species.name,"_results_compar.models",sep="/"))
print(paste('Algo:',species.name, '-',bg,'-',nb.var,'-',reduction,'-',method.name,sep=' '))
resampling <-  resampling_method

if (method.name == "KNN") {
	if (resampling == "boot") p1 <- b.k1  
	if (resampling == "cv")   p1 <- c.k1 
	p2 <- 0; p3 <- 0; p.lab <- paste(p1,sep="_")
}
if (method.name == "SVM") {
	if (resampling == "boot") {p1 <- b.s1; p2 <- b.s2}
	if (resampling == "cv")   {p1 <- c.s1; p2 <- c.s2}
	p3 <- 0; p.lab <- paste(p1,p2,sep="_")
}
if (method.name == "NNET") {
	if (resampling == "boot") {p1 <- b.n1; p2 <- b.n2; p3 <- b.n3; nnet.Wts <- nnet.Wts.boot}
	if (resampling == "cv")   {p1 <- c.n1; p2 <- c.n2; p3 <- c.n3; nnet.Wts <- nnet.Wts.cv}
	p.lab <- paste(p1,p2,p3,sep="_")
}
if (method.name != "KNN" & method.name != "SVM" & method.name != "NNET") {
	p1 <- p2 <- p3 <- 0
	p.lab <- ""
}


#sourceDir(func.dir)
train.and.test.model.predict <- train.and.test.model(train=data_coord_order[,-c(1,2)],
                                 method=method.name,
                                 test=worldclim2,
                                 test2=NULL,
                                 test3=NULL,
                                 resp.var.name="pres", 
                                 knn.par = p1, 
                                 svm.par = c("rbfdot",p1,p2), 
                                 nnet.par = c(p1,p2,p3),
                                 Wts = nnet.Wts,
                                 nbtree=1000, #put back to 10000
                                 priors=NULL,
                                 prob=T)


setwd(paste(main_wd,species.name,"_prediction_maps/_predict_world",sep="/"))
save(train.and.test.model.predict,file=paste("World_predictions_",species.name,"_",method.name,"_",p.lab,"_",resampling,"_",bg,nb.var,"_",reduction,".Rdata",sep=""))

# Results
res <- cbind(train.and.test.model.predict$pred.test,train.and.test.model.predict$prob.test)
colnames(res) <- c("Prediction","Probability")

plot.split(species.name,method.name,save.fig=F)
savePlot(paste("World_prediction_",species.name,"_",method.name,"_",p.lab,"_",resampling,"_",bg,nb.var,"_",reduction,".png",sep=""),type="png")

# ------------------------ #
# Make database for ArcGIS
# ------------------------ #
# Save things greater than threshold
thresh <- 0
dat.prob <- cbind(x[res[,2]>=thresh,1],x[res[,2]>=thresh,2],res[res[,2]>=thresh,2])
site=paste("site",1:length(res[,1]),sep="")
dat.prob <-data.frame(dat.prob,site)
colnames(dat.prob)[1:3] <- c("xcoord", "ycoord","Prob")write.table(dat.prob, paste("World_",species.name,"_",method.name,"_",p.lab,"_threshold_",thresh,"_",resampling,"_data_prob_",bg,nb.var,"_",reduction,".txt",sep=""), row.names=F)

dat.vote <- data.frame(x,res[,1],site=paste("site",1:length(res[,1]),sep=""))
colnames(dat.vote)[1:3] <- c("xcoord", "ycoord","Vote")
write.table(dat.vote, paste("World_",species.name,"_",method.name,"_",p.lab,"_",resampling,"_data_vote1_",bg,nb.var,"_",reduction,".txt",sep=""), row.names=F)
graphics.off()
}

###################################################################################################################################
setwd(paste(main_wd,species.name,"_dist_data_OCSVM",sep="/")) 
graphics.off()
gc()
#}# reduction loop
#}# bg loop



