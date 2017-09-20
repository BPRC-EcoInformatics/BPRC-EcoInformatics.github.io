library(ecospat)
library(ggplot2)

species_name<-"Ccapitata"

main_wd<-"C:/Users/Ursula Torres/Documents/Ecoinformatics_workshop"
setwd(paste(main_wd,"/6_Extrapolation",sep=""))
dir.create(species_name)
setwd(species_name)

###################################################################################################
###########Compute exdet index between calibration and projection dataset            ##############
########### Exdet index quantifies extrapolation. Tool developed by Mesgaran 2014    ##############
###################################################################################################

#prepare calibration data
calibration_data<-read.table(paste(main_wd, "/1_CleaningData/",species_name,"/",species_name,"_clim_range_withoutRep.txt",sep=""),h=T,sep='\t')
choose.var<-c(1,10,25)#use the variables selected in the MMA framework
calibration_data<-calibration_data[,choose.var]

#prepare projection data
load(paste(main_wd,"/4_MMA/_database/P2DR1.Rdata",sep=""))
worldclim<-get("P2DR1")
projection_data<-worldclim[,choose.var]#choose same variables as the calibration data

#compute extrapolation metric
exdet<-ecospat.exdet(calibration_data,projection_data)

#join coordinates of the projection data
exdet<-cbind(exdet,worldclim[,c("xcoord","ycoord")])

###################################################################################################
###########Plot extrapolation map            ######################################################
###########                                 ######################################################
###################################################################################################
breaks_scale<-c(0.01,0.1,0.33,0.66,0.9,0.99)
scalerange <- range(exdet$exdet)
gradientends<-c(scalerange[1],scalerange[1]/2,scalerange[2]/2,scalerange[2])
colorends<-c("blue","green","red")




png(paste(species_name,"_extrapolation.png",sep=""),width=936,height= 455)
p1<-ggplot(data=exdet, aes(y=ycoord, x=xcoord)) +
  geom_raster(aes(fill=exdet))+
  scale_fill_gradientn(name="exdet",colours=c("#0000FFFF","#FFFFFFFF","#666600"))+
  labs(title="Exdet map",x= "Longitude", y = "Latitude")+
  theme(panel.background = element_blank())
p1
dev.off()


png(paste(species_name,"_extrapolation.png",sep=""),width=936,height= 455)
p1<-ggplot(data=exdet, aes(y=ycoord, x=xcoord)) +
  geom_raster(aes(fill=exdet))+
  scale_fill_gradientn(name="exdet",colours=colorends,values=gradientends)+
  labs(title="Exdet map",x= "Longitude", y = "Latitude")+
  theme(panel.background = element_blank())
p1
dev.off()