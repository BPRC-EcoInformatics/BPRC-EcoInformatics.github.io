library("ade4", lib.loc="C:/Users/ytakeuchi/R")
library("sp", lib.loc="C:/Users/ytakeuchi/R")
library("maptools", lib.loc="C:/Users/ytakeuchi/R")
library("rgdal", lib.loc="C:/Users/ytakeuchi/R")
library("raster", lib.loc="C:/Users/ytakeuchi/R")
library("rgeos", lib.loc="C:/Users/ytakeuchi/R")
library("RColorBrewer", lib.loc="C:/Users/ytakeuchi/R")
library("plyr", lib.loc="C:/Users/ytakeuchi/R")
library("reshape2", lib.loc="C:/Users/ytakeuchi/R")
library("ggplot2", lib.loc="C:/Users/ytakeuchi/R")
library("bindrcpp", lib.loc="C:/Users/ytakeuchi/R")
library("dplyr", lib.loc="C:/Users/ytakeuchi/R")



main_wd <- "C:/Users/ytakeuchi/APHIS/Ecoinformatics_workshop" 
#read the database
rdata.name <- "P2DR1"
load(paste(main_wd,"/4_MMA/_database/",rdata.name,".Rdata",sep=""))###### UT added / before database
species_name<-"Ccapitata"

#read the pca output from the niche analysis

load("C:/Users/ytakeuchi/APHIS/Ecoinformatics_workshop/2_NicheAnalysis/Ccapitata/continent/temp_precipitation_moisture_elevation/pca_Ccapitata.rda")

#read niche grid from the niche analysis
load("C:/Users/ytakeuchi/APHIS/Ecoinformatics_workshop/2_NicheAnalysis/Ccapitata/continent/temp_precipitation_moisture_elevation/grid_z1_Ccapitata_100_0.025_none_.rda")
load("C:/Users/ytakeuchi/APHIS/Ecoinformatics_workshop/2_NicheAnalysis/Ccapitata/continent/temp_precipitation_moisture_elevation/grid_z2_Ccapitata_100_0.025_none_.rda")

#read niche grid from the niche analysis
load("C:/Users/ytakeuchi/APHIS/Ecoinformatics_workshop/2_NicheAnalysis/Ccapitata/continent/temp_precipitation_moisture_elevation/scores_Ccapitata.rda")

scores.sp1<-scores$scores.sp1
scores.sp2<-scores$scores.sp2
scores.clim12<-scores$scores.clim12
scores.clim1<-scores$scores.clim1
scores.clim2<-scores$scores.clim2



#main_wd <-"/home/audrey/Desktop/USDA"
shape_file<-readShapePoly(paste(main_wd, "/0_Data/biomes/test_intersect.shp",sep=''))
coord_data<-P2DR1[,1:2]
coordinates(coord_data)= c("xcoord", "ycoord")
coord_data1<-coord_data


points = SpatialPointsDataFrame(coord_data,P2DR1[,1:2])
o = over(points,shape_file[,c("COUNTRY","CONTINENT")])
points$country <- o$"COUNTRY"

clim.aus <- P2DR1[which(as.character(points$country)=='Australia'),]
clim.nz <- P2DR1[which(as.character(points$country)=='New Zealand'),]
clim.usat <- P2DR1[which(as.character(points$country)=='United States'),]
clim.usaf <- P2DR1[which(as.character(points$country)=='Canada'),]
clim.usa <- rbind(clim.usat,clim.usaf)
clim.fr <- P2DR1[which(as.character(points$country)=='France'),] 


# Projet the climate of studied area into pca space
temp.usa<-suprow(pca.cal, scale(clim.usa[3:dim(clim.usa)[2]], center=pca.cal$cent, scale=pca.cal$norm))
scores.usa<-cbind(clim.usa,temp.usa$lisup)

temp.aus<-suprow(pca.cal,scale(clim.aus[3:dim(clim.aus)[2]], center=pca.cal$cent, scale=pca.cal$norm))
scores.aus<-cbind(clim.aus,temp.aus$lisup)

temp.nz<-suprow(pca.cal,scale(clim.nz[3:dim(clim.nz)[2]], center=pca.cal$cent, scale=pca.cal$norm))
scores.nz<-cbind(clim.nz,temp.nz$lisup)

temp.fr<-suprow(pca.cal,scale(clim.fr[3:dim(clim.fr)[2]], pca.cal$cent, pca.cal$norm))
scores.fr<-cbind(clim.fr,temp.fr$lisup)

temp.usa<-suprow(pca.cal, clim.usa[3:dim(clim.usa)[2]],center = T, scale = T)
scores.usa<-cbind(clim.usa,temp.usa$lisup)

temp.aus<-suprow(pca.cal,clim.aus[3:dim(clim.aus)[2]],center = T, scale = T)
scores.aus<-cbind(clim.aus,temp.aus$lisup)

temp.nz<-suprow(pca.cal,clim.nz[3:dim(clim.nz)[2]],center = T, scale = T)
scores.nz<-cbind(clim.nz,temp.nz$lisup)

temp.fr<-suprow(pca.cal,clim.fr[3:dim(clim.fr)[2]],center = T, scale = T)
scores.fr<-cbind(clim.fr,temp.fr$lisup)

temp.usa<-suprow(pca.cal, clim.usa[3:dim(clim.usa)[2]])
scores.usa<-cbind(clim.usa,temp.usa$lisup)

temp.aus<-suprow(pca.cal,clim.aus[3:dim(clim.aus)[2]])
scores.aus<-cbind(clim.aus,temp.aus$lisup)

temp.nz<-suprow(pca.cal,clim.nz[3:dim(clim.nz)[2]])
scores.nz<-cbind(clim.nz,temp.nz$lisup)

temp.fr<-suprow(pca.cal,clim.fr[3:dim(clim.fr)[2]])
scores.fr<-cbind(clim.fr,temp.fr$lisup)

#extracting the niche categories for the country climate
Z <- t(as.matrix(z1$w + 2 * z2$w))[,nrow(as.matrix(z1$z.uncor)):1]

m <- matrix(c(1:4), ncol = 4, byrow = TRUE)
layout (m)
colz1 = "#00FF0050"; colz2 = "#FF000050"; colinter = "#0000FF50"; colZ1 = "green3"; colZ2 = "red3";

image(x=z2$x,y=z2$y,z=Z, col = c("#FFFFFF00", colz1, colz2, colinter))
points(scores.sp1[,c("Axis1")],scores.sp1[,c("Axis2")],col="green",pch=16)#scores occ native range
points(scores.sp2[,c("Axis1")],scores.sp2[,c("Axis2")],col="red",pch=16)#scores occ invasive range
points(scores.usa[,c("Axis1")],scores.usa[,c("Axis2")],col="black",pch=16)#scores env nz

image(x=z2$x,y=z2$y,z=Z, col = c("#FFFFFF00", colz1, colz2, colinter))
points(scores.sp1[,c("Axis1")],scores.sp1[,c("Axis2")],col="green",pch=16)#scores occ native range
points(scores.sp2[,c("Axis1")],scores.sp2[,c("Axis2")],col="red",pch=16)#scores occ invasive range
points(scores.aus[,c("Axis1")],scores.aus[,c("Axis2")],col="black",pch=16)#scores env nz

image(x=z2$x,y=z2$y,z=Z, col = c("#FFFFFF00", colz1, colz2, colinter))
points(scores.sp1[,c("Axis1")],scores.sp1[,c("Axis2")],col="green",pch=16)#scores occ native range
points(scores.sp2[,c("Axis1")],scores.sp2[,c("Axis2")],col="red",pch=16)#scores occ invasive range
points(scores.nz[,c("Axis1")],scores.nz[,c("Axis2")],col="black",pch=16)#scores env nz

image(x=z2$x,y=z2$y,z=Z, col = c("#FFFFFF00", colz1, colz2, colinter))
points(scores.sp1[,c("Axis1")],scores.sp1[,c("Axis2")],col="green",pch=16)#scores occ native range
points(scores.sp2[,c("Axis1")],scores.sp2[,c("Axis2")],col="red",pch=16)#scores occ invasive range
points(scores.fr[,c("Axis1")],scores.fr[,c("Axis2")],col="black",pch=16)#scores env nz



#generate matrix with categories (stability, unfilling and expansion)
library(raster)

# Create niche categories
Z <- t(as.matrix(z1$w + 2 * z2$w))[,nrow(as.matrix(z1$z.uncor)):1]

Xusa <- sapply(scores.usa$Axis1, findInterval, z2$x)
Yusa <- sapply(scores.usa$Axis2, findInterval, z2$y)
niche.usa <- matrix(0,ncol=3,nrow=length(Xusa))
for (i in seq(1,length(Xusa))){
	if(Xusa[i]!=0 & Yusa[i]!=0) niche.usa[i,]<-c(clim.usa[i,1],clim.usa[i,2],Z[Xusa[i],Yusa[i]])
	if(Xusa[i]==100){niche.usa[i,]<-c(clim.usa[i,1],clim.usa[i,2],0)}
	if(Yusa[i]==100){niche.usa[i,]<-c(clim.usa[i,1],clim.usa[i,2],0)}
}
scores.usa$niche_class<-as.factor(niche.usa[,3])
scores.usa[is.na(scores.usa$niche_class),c("niche_class")]<-0
scores.usa = within(scores.usa, {
  stability = ifelse(niche_class == 3, 1, 0)
  expansion = ifelse(niche_class == 2, 1, 0)
  unfilling = ifelse(niche_class == 1, 1, 0)
  })


Xaus <- sapply(scores.aus$Axis1, findInterval, z2$x)
Yaus <- sapply(scores.aus$Axis2, findInterval, z2$y)
niche.aus <- matrix(0,ncol=3,nrow=length(Xaus))
for (i in seq(1,length(Xaus))){
	if(Xaus[i]!=0 & Yaus[i]!=0) niche.aus[i,]<-c(clim.aus[i,1],clim.aus[i,2],Z[Xaus[i],Yaus[i]])
	if(Xaus[i]==100){niche.aus[i,]<-c(clim.aus[i,1],clim.aus[i,2],0)}
	if(Yaus[i]==100){niche.aus[i,]<-c(clim.aus[i,1],clim.aus[i,2],0)}
}
scores.aus$niche_class<-as.factor(niche.aus[,3])
scores.aus[is.na(scores.aus$niche_class),c("niche_class")]<-0
scores.aus = within(scores.aus, {
  stability = ifelse(niche_class == 3, 1, 0)
  expansion = ifelse(niche_class == 2, 1, 0)
  unfilling = ifelse(niche_class == 1, 1, 0)
  })

Xnz <- sapply(scores.nz$Axis1, findInterval, z2$x)
Ynz <- sapply(scores.nz$Axis2, findInterval, z2$y)
niche.nz <- matrix(0,ncol=3,nrow=length(Xnz))
for (i in seq(1,length(Xnz))){
	if(Xnz[i]!=0 & Ynz[i]!=0) niche.nz[i,]<-c(clim.nz[i,1],clim.nz[i,2],Z[Xnz[i],Ynz[i]])
	if(Xnz[i]==100){niche.nz[i,]<-c(clim.nz[i,1],clim.nz[i,2],0)}
	if(Ynz[i]==100){niche.nz[i,]<-c(clim.nz[i,1],clim.nz[i,2],0)}
}
scores.nz$niche_class<-as.factor(niche.nz[,3])
scores.nz[is.na(scores.nz$niche_class),c("niche_class")]<-0
scores.nz = within(scores.nz, {
  stability = ifelse(niche_class == 3, 1, 0)
  expansion = ifelse(niche_class == 2, 1, 0)
  unfilling = ifelse(niche_class == 1, 1, 0)
  })

Xfr <- sapply(scores.fr$Axis1, findInterval, z2$x)
Yfr <- sapply(scores.fr$Axis2, findInterval, z2$y)
niche.fr <- matrix(0,ncol=3,nrow=length(Xfr))
for (i in seq(1,length(Xfr))){
	if(Xfr[i]!=0 & Yfr[i]!=0) niche.fr[i,]<-c(clim.fr[i,1],clim.fr[i,2],Z[Xfr[i],Yfr[i]])
	if(Xfr[i]==100){niche.fr[i,]<-c(clim.fr[i,1],clim.fr[i,2],0)}
	if(Yfr[i]==100){niche.fr[i,]<-c(clim.fr[i,1],clim.fr[i,2],0)}
}
scores.fr$niche_class<-as.factor(niche.fr[,3])
scores.fr[is.na(scores.fr$niche_class),c("niche_class")]<-0
scores.fr = within(scores.fr, {
  stability = ifelse(niche_class == 3, 1, 0)
  expansion = ifelse(niche_class == 2, 1, 0)
  unfilling = ifelse(niche_class == 1, 1, 0)
  })


#######################Plot raster 
#Create colours
library(RColorBrewer)
library(ggplot2)
library(dplyr)
mycolours<-brewer.pal(4,"Paired")
#native_colours<-mycolours[2]
native_colours<-"#009E73"
invasive_colours<-c("gold3")
inter<-mycolours[2]
colours_niche<-c("#CCCCCC",native_colours,invasive_colours,inter)
names(colours_niche)<-c("0","1","2","3")

label_names<-c("Outside of the known niche", "Unfilling", "Expansion","Stability")
legend_names<-data.frame(label_names,cat_code=c("0","1","2","3"))

scores.usa<-droplevels(scores.usa)
cat_code.usa=data.frame(cat_code=as.factor(levels(scores.usa$niche_class)))
legend_names_sp.usa<-right_join(legend_names,cat_code.usa,by="cat_code")

scores.aus<-droplevels(scores.aus)
cat_code.aus=data.frame(cat_code=as.factor(levels(scores.aus$niche_class)))
legend_names_sp.aus<-right_join(legend_names,cat_code.aus,by="cat_code")

scores.nz<-droplevels(scores.nz)
cat_code.nz=data.frame(cat_code=as.factor(levels(scores.nz$niche_class)))
legend_names_sp.nz<-right_join(legend_names,cat_code.nz,by="cat_code")

scores.fr<-droplevels(scores.fr)
cat_code.fr=data.frame(cat_code=as.factor(levels(scores.fr$niche_class)))
legend_names_sp.fr<-right_join(legend_names,cat_code.fr,by="cat_code")


#Plot
library("gridExtra", lib.loc="C:/Users/ytakeuchi/R")
library("labeling", lib.loc="C:/Users/ytakeuchi/R")
library("digest", lib.loc="C:/Users/ytakeuchi/R")


p1<-ggplot(data=scores.usa, aes(y=ycoord, x=xcoord)) +
  geom_raster(aes(fill=as.factor(niche_class)))+
  scale_fill_manual(name="Niche categories",values=c(colours_niche),labels=legend_names_sp.usa[,1])+
  labs(x= "Longitude", y = "Latitude")+
  theme(panel.background = element_blank())

p2<-ggplot(data=scores.aus, aes(y=ycoord, x=xcoord)) +
  geom_raster(aes(fill=as.factor(niche_class)))+
  scale_fill_manual(name="Niche categories",values=c(colours_niche),labels=legend_names_sp.aus[,1])+
  labs(x= "Longitude", y = "Latitude")+
  theme(panel.background = element_blank())

p3<-ggplot(data=scores.nz, aes(y=ycoord, x=xcoord)) +
  geom_raster(aes(fill=as.factor(niche_class)))+
  scale_fill_manual(name="Niche categories",values=c(colours_niche),labels=legend_names_sp.nz[,1])+
  labs(x= "Longitude", y = "Latitude")+
  theme(panel.background = element_blank())

p4<-ggplot(data=scores.fr, aes(y=ycoord, x=xcoord)) +
  geom_raster(aes(fill=as.factor(niche_class)))+
  scale_fill_manual(name="Niche categories",values=c(colours_niche),labels=legend_names_sp.fr[,1])+
  labs(x= "Longitude", y = "Latitude")+
  theme(panel.background = element_blank())
x11()
png(filename = paste('C:/Users/ytakeuchi/APHIS/Ecoinformatics_workshop/2_NicheAnalysis/',species_name,"nichedynamics.png",sep="_"), width = 1600, height = 1000,res=130)
grid.arrange(p1,p4,p2,p3,ncol=2)
dev.off()

