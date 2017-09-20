niche_analysis_function <- function(main_wd, species_name,native,invasive,type_analysis,number_variable,PROJ,iterations,R){


#######################################################################################################################################
# Type of analysis - choose the extent on which the environemntal background will calculated 
# it can be countries, continent or ecoregion
cat('Setting up directories'); cat('\n')

setwd(paste(main_wd, "/2_NicheAnalysis/",species_name,sep=''))
dir.create(type_analysis)
setwd(type_analysis)
if(number_variable==29){var_name='temp_precipitation'}
if(number_variable==37){var_name='temp_precipitation_moisture_elevation'}
dir.create(var_name)
setwd(var_name)


Xvar<-c(3:number_variable)


# selection of the type of analysis.
# If PROJ =F, the models are calibrated on both ranges.
# If PROJ =T, the models are calibrated on species 1 range only and projected to range 2. 
# Analyses where both ranges are needed (ex: LDA) are not done

cat('Load environmental variables'); cat('\n')
# load environemental background for all sites in the native range (column names should be x,y,X1,X2,...,Xn)
clim1<-read.table(paste(main_wd, "/1_CleaningData/",species_name,"/mask_",species_name,"_", type_analysis ,"_native_df.txt",sep=""),h=T,sep=' ')
colnames(clim1)<-c("x","y", paste("cbio",1:35,sep=""))

# load environemental background for all sites in the non-native range (column names should be x,y,X1,X2,...,Xn)
clim2<-read.table(paste(main_wd, "/1_CleaningData/",species_name,"/mask_",species_name,"_", type_analysis ,"_invasive_df.txt",sep=""),h=T,sep=' ')
colnames(clim2)<-c("x","y", paste("cbio",1:35,sep=""))

# Global climate for both ranges
clim12<-rbind(clim1,clim2)

# Occurence data for native and invasive range 
occ.sp1<-native[,1:37]
colnames(occ.sp1)[1:2]<-c("x","y")
occ.sp2<-invasive[,1:37]
colnames(occ.sp2)[1:2]<-c("x","y")


# if PROJ = F
row.w.1.occ<-1-(nrow(occ.sp1)/nrow(rbind(occ.sp1,occ.sp2))) # prevalence of occ1 (weight by the number of occurences recorded)
row.w.2.occ<-1-(nrow(occ.sp2)/nrow(rbind(occ.sp1,occ.sp2))) # prevalence of occ2
row.w.occ<-c(rep(0, nrow(clim1)),rep(0, nrow(clim2)),rep(row.w.1.occ, nrow(occ.sp1)),rep(row.w.2.occ, nrow(occ.sp2)))

row.w.1.env<-1-(nrow(clim1)/nrow(clim12))  # prevalence of clim1
row.w.2.env<-1-(nrow(clim2)/nrow(clim12))  # prevalence of clim2
row.w.env<-c(rep(row.w.1.env, nrow(clim1)),rep(row.w.2.env, nrow(clim2)),rep(0, nrow(occ.sp1)),rep(0, nrow(occ.sp2)))

fac<-as.factor(c(rep(1, nrow(clim1)),rep(2, nrow(clim2)),rep(1, nrow(occ.sp1)),rep(2, nrow(occ.sp2))))

# if PROJ = T
row.w.occ.PROJT<-c(rep(0, nrow(clim1)),rep(0, nrow(clim2)),rep(1, nrow(occ.sp1)),rep(0, nrow(occ.sp2)))
row.w.env.PROJT<-c(rep(1, nrow(clim1)),rep(0, nrow(clim2)),rep(0, nrow(occ.sp1)),rep(0, nrow(occ.sp2)))

# global dataset for the analysis and rows for each sub dataset
clim1=as.data.frame(sapply(clim1,as.numeric))
clim2=as.data.frame(sapply(clim2,as.numeric))

data.env.occ<-rbind(clim1[,1:37],clim2[,1:37],occ.sp1,occ.sp2)[Xvar]
row.clim1<-1:nrow(clim1)
row.clim2<-(nrow(clim1)+1):(nrow(clim1)+nrow(clim2))
row.clim12<-1:(nrow(clim1)+nrow(clim2))
row.sp1<-(nrow(clim1)+nrow(clim2)+1):(nrow(clim1)+nrow(clim2)+nrow(occ.sp1))
row.sp2<-(nrow(clim1)+nrow(clim2)+nrow(occ.sp1)+1):(nrow(clim1)+nrow(clim2)+nrow(occ.sp1)+nrow(occ.sp2))



######################################################################################################################################
# measures niche overlap along the two first axes of a PCA calibrated on all the pixels of the study areas
cat('Principal component analysis'); cat('\n')
if(PROJ == F){  #fit of the analyse using occurences from both ranges		
  pca.cal <-dudi.pca(data.env.occ,row.w = row.w.env, center = T, scale = T, scannf = F, nf = 2)
}
if(PROJ == T){	#fit of the analyse using occurences from range 1		
  pca.cal <-dudi.pca(data.env.occ,row.w = row.w.env.PROJT, center = T, scale = T, scannf = F, nf = 2)
}

ecospat.plot.contrib(pca.cal$co,pca.cal$eig)
pca.cal$co[order(-abs(pca.cal$co[,1])),]

#remove files no longer needed
rm(list=c("data.env.occ"))

save(pca.cal, file=paste("pca_",species_name,".rda",sep=""))

# predict the scores on the axes
scores.clim12<- pca.cal$li[row.clim12,]
scores.clim1<- pca.cal$li[row.clim1,]
scores.clim2<- pca.cal$li[row.clim2,]
scores.sp1<- pca.cal$li[row.sp1,]
scores.sp2<- pca.cal$li[row.sp2,]

#save the scores
scores<-list(scores.clim12,scores.clim1,scores.clim2,scores.sp1,scores.sp2)
names(scores)<-c("scores.clim12","scores.clim1","scores.clim2","scores.sp1","scores.sp2")
save(scores,file=paste("scores_",species_name,".rda",sep=""))

rm(list=c("scores"))

######################################################################################################################################
# calculation of occurence density (kernel in the PCA space)
cat('Calculating occurence density '); cat('\n')
th.sp<-0.025 # Quantile to use of the envrionemntal density used to delimit marginal climates
th.env<-"none"
z1<-ecospat.grid.clim.dyn(scores.clim12,scores.clim1,scores.sp1,R,th.sp=th.sp)
z2<-ecospat.grid.clim.dyn(scores.clim12,scores.clim2,scores.sp2,R,th.sp,th.sp)

save(z1,file=paste("grid_z1",species_name,R,th.sp,th.env,".rda",sep="_"))
save(z2,file=paste("grid_z2",species_name,R,th.sp,th.env,".rda",sep="_"))

cat('Calculating index expansion/stability/unfilling'); cat('\n')
# # calculating the index expansion/stability/unfilling and classical D metric and Hellingers metric
# for different uantile used for delimiting marginal climates when calculating the envrionemntal density 
index_0<-ecospat.niche.dyn.index(z1,z2,intersection=0)
index_5th<-ecospat.niche.dyn.index(z1,z2,intersection=0.05)
index_10th<-ecospat.niche.dyn.index(z1,z2,intersection=0.10)
index_15th<-ecospat.niche.dyn.index(z1,z2,intersection=0.15)
index_20th<-ecospat.niche.dyn.index(z1,z2,intersection=0.20)
index_25th<-ecospat.niche.dyn.index(z1,z2,intersection=0.25)
index_30th<-ecospat.niche.dyn.index(z1,z2,intersection=0.30)
index_35th<-ecospat.niche.dyn.index(z1,z2,intersection=0.35)
index_40th<-ecospat.niche.dyn.index(z1,z2,intersection=0.40)
index_45th<-ecospat.niche.dyn.index(z1,z2,intersection=0.45)
index_50th<-ecospat.niche.dyn.index(z1,z2,intersection=0.50)
index_55th<-ecospat.niche.dyn.index(z1,z2,intersection=0.55)
index_60th<-ecospat.niche.dyn.index(z1,z2,intersection=0.60)
index_65th<-ecospat.niche.dyn.index(z1,z2,intersection=0.65)
index_70th<-ecospat.niche.dyn.index(z1,z2,intersection=0.70)
index_75th<-ecospat.niche.dyn.index(z1,z2,intersection=0.75)
index_80th<-ecospat.niche.dyn.index(z1,z2,intersection=0.80)
index_85th<-ecospat.niche.dyn.index(z1,z2,intersection=0.85)
index_90th<-ecospat.niche.dyn.index(z1,z2,intersection=0.90)
index_95th<-ecospat.niche.dyn.index(z1,z2,intersection=0.95)
index_whole<-ecospat.niche.dyn.index(z1,z2,intersection=NA)

overlap<-ecospat.niche.overlap(z1,z2,cor=T)

#testing significance for D metric and Hellinger's metric
eq<-ecospat.niche.equivalency.test(z1,z2,iterations)
sim1_2<-ecospat.niche.similarity.test(z1,z2,iterations)#,one.sided=F)
sim2_1<-ecospat.niche.similarity.test(z2,z1,iterations)#,one.sided=F)

save(eq,file=paste("eq",species_name,"R",R,"it",iterations,".rda",sep="_"))
save(sim1_2,file=paste("sim1_2",species_name,"R",R,"it",iterations,".rda",sep="_"))
save(sim2_1,file=paste("sim2_1",species_name,"R",R,"it",iterations,".rda",sep="_"))

new_index<-list(index_0,index_5th,index_10th,index_15th,index_20th,index_25th,index_30th,index_35th,index_40th,index_45th,index_50th,index_55th,index_60th,index_65th,index_70th,index_75th,index_80th
,index_85th,index_90th,index_95th,index_whole)
names(new_index)<-c("index_0","index_5th","index_10th","index_15th","index_20th","index_25th","index_30th","index_35th","index_40th","index_45th","index_50th","index_55th","index_60th","index_65th","index_70th","index_75th","index_80th"
,"index_85th","index_90th","index_95th","index_whole")
results_index<-do.call("rbind",lapply(new_index,"[[",2))

D<-c(eq$obs$D,eq$p.D)
names(D)<-c("index","pvalue")
I<-c(eq$obs$I,eq$p.I)
names(I)<-c("index","pvalue")

#D<-overlap$D
#I<-overlap$I

results<-list(results_index,D,I)
names(results)<-c("dynamic index","D","I")

cat('Save results'); cat('\n')
save(results,file=paste("index",species_name,"R",R,"it",iterations,".rda",sep="_"))

######################################################################################################################################
# Plot
#create colours
cat('Plotting niche space'); cat('\n')
add.alpha <- function(col, alpha=1){
if(missing(col))
stop("Please provide a vector of colours.")
apply(sapply(col, col2rgb)/255, 2,
function(x)
rgb(x[1], x[2], x[3], alpha=alpha))
}

new_brown<-"burlywood4"
new_brown_alpha<- add.alpha(new_brown, alpha=0.4)

new_green<-add.alpha("chartreuse4",alpha=0.4)
new_red<-add.alpha("firebrick4",alpha=0.4)

# density grids
#pdf(paste(species_name,"densgrid.pdf",sep="_"),width = 13, height = 9 )
#m <- matrix(c(1:2), ncol = 2, byrow = TRUE)
#layout (m)

#ecospat.plot.niche (z1, title=paste(species_name,type_analysis,sep=" "), name.axis1="PC1", name.axis2="PC2", cor=T)
#ecospat.plot.niche (z2, title=paste(species_name,type_analysis,sep=" "), name.axis1="PC1", name.axis2="PC2", cor=T)
#dev.off()


#dynamic index
#pdf(paste(species_name,"nichegrid_quantile25.pdf",sep="_"),width = 13, height = 9 )	
png(filename = paste(main_wd,'/2_NicheAnalysis/',species_name,"/",species_name,"_nichegrid_quantile25.png",sep=""), width = 1600, height = 1000,res=130)

m <- matrix(c(1:2), ncol = 2, byrow = TRUE)
layout (m)

# quant = 75th quantile - Quantile to use of the environemntal density used to delimit marginal climates
# if interest=1 plot native density
# if interest=2 p[lot non-native density
# colz1: unfilling area, colz2: expension area, colinter: overlap, colZ1 native extent, colZ2: non-native extent


if (species_name=='Ldispar'){name_title='L. dispar'}
if (species_name=='Hvitripennis'){name_title='H. vitripennis'}
if (species_name=='Ccapitata'){name_title='C. capitata'}



ecospat.plot.niche.dyn(z1,z2,0.25,title='',interest=1,colz1=new_green,colz2=new_red,colinter=new_brown_alpha,colZ1=new_green,colZ2=new_red)
#legend("topright" ,inset=c(-0.2,0),  fill = c(new_green,new_red,new_brown_alpha,0,0,0,0), lty = c(0,0,0,1,2,1,2), c('Unfilling area','Expansion area', 'Overlap area', 'Native extent', 'Native extent without marginal envrionments','Invaded extent','Invaded extent without marginal envrionments'),  col = c('white','white','white',new_green, new_green,new_red, new_red), merge = T, border=rep('white',7))
legend(min(pca.cal$li[,1])+1, max(pca.cal$li[,2])-1 ,  fill = c(new_green,new_red,new_brown_alpha,0,0,0,0), lty = c(0,0,0,1,2,1,2), c('Unfilling area','Expansion area', 'Overlap area', 'Native extent', 'Native extent without marginal envrionments','Invaded extent','Invaded extent without marginal envrionments'),  col = c('white','white','white',new_green, new_green,new_red, new_red), merge = T, border=rep('white',7))

title(sub = paste('Density kernel: Native range - Environmental background: ',type_analysis,sep=''), cex=0.1, font.main= 6)
ecospat.plot.niche.dyn(z1,z2,0.25,title='',interest=2,
colz1=new_green,colz2=new_red,colinter=new_brown_alpha,colZ1=new_green,colZ2=new_red)
title(sub = paste('Density kernel: Invaded range - Environmental background: ',type_analysis,sep=''), cex=0.1, font.main= 6)
title(main = paste(name_title,sep=''), outer=TRUE, line=-1.5)
title(main = paste('Niche analysis: unfilling - ', round(index_25th$dynamic.index.w[3],digit=3), ', expansion - ', round(index_25th$dynamic.index.w[1],digit=3),', stability - ', round(index_25th$dynamic.index.w[2],digit=3),sep=''), cex=0.1,outer=TRUE, line=-3, font.main= 6)
dev.off()

pdf(paste(species_name,"nichegrid_noquantile.pdf",sep="_"),width = 13, height = 9 )	
m <- matrix(c(1:2), ncol = 2, byrow = TRUE)
layout (m)

#whole niche
ecospat.plot.niche.dyn(z1,z2,1,title='',interest=1,
colz1=new_green,colz2=new_red,colinter=new_brown_alpha,colZ1=new_green,colZ2=new_red)
legend(min(pca.cal$li[,1])+1, max(pca.cal$li[,2])-1 ,  fill = c(new_green,new_red,new_brown_alpha,0,0), lty = c(0,0,0,1,1), c('Unfilling area','Expansion area', 'Overlap area', 'Native extent', 'Invaded extent'),  col = c('white','white','white',new_green, new_red), merge = T, border=rep('white',5))
title(sub = paste('Density kernel: Native range - Envrionmental background: ',type_analysis,sep=''), cex=0.1, font.main= 6)

ecospat.plot.niche.dyn(z1,z2,1,title='',interest=2,
colz1=new_green,colz2=new_red,colinter=new_brown_alpha,colZ1=new_green,colZ2=new_red)
#title(main = paste('Envrionemntal background: ',type_analysis,sep=''), line=0.1, cex=0.1, font.main= 6)
#title(main = paste('Density kernel: Invasive range',sep=''), line=1.1, cex=0.1, font.main= 6)
title(sub = paste('Density kernel: Invaded range - Envrionmental background: ',type_analysis,sep=''), cex=0.1, font.main= 6)
title(main = paste(name_title,sep=''), outer=TRUE, line=-1.5)
title(main = paste('Niche analysis: unfilling - ', round(index_whole$dynamic.index.w[3],digit=3), ', expansion - ', round(index_whole$dynamic.index.w[1],digit=3),', stability - ', round(index_whole$dynamic.index.w[2],digit=3),sep=''), cex=0.1,outer=TRUE, line=-3, font.main= 6)
dev.off()

#correlation circle and explained variance
pdf(paste(species_name,"correlation_circle_PCAenv.pdf",sep="_"))	
ecospat.plot.contrib(pca.cal$co,pca.cal$eig)
dev.off()

rm(list=c("pca.cal","z1","z2","scores.clim12","scores.clim1","scores.clim2"
,"scores.sp1","scores.sp2","index_0","index_5th","index_10th","index_15th","index_20th","index_25th","index_30th","index_35th","index_40th","index_45th","index_50th","index_55th","index_60th","index_65th","index_70th","index_75th","index_80th"
,"index_85th","index_90th","index_95th","index_whole","eq","sim1_2","sim2_1"))

#rm(list=c("pca.cal","z1","z2","scores.clim12","scores.clim1","scores.clim2"
#,"scores.sp1","scores.sp2","index_0","index_5th","index_10th","index_15th","index_20th","index_25th","index_30th","index_35th","index_40th","index_45th","index_50th","index_55th","index_60th","index_65th","index_70th","index_75th","index_80th"
#,"index_85th","index_90th","index_95th","index_whole"))

setwd(paste(main_wd, "/2_NicheAnalysis/",species_name,sep=''))

}



