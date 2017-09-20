##############################################
##########SPECIES ENSEMBLE MAPS ##############
##############################################
library(tidyr)
library(ade4)
library(reshape2)
library(plyr)
library(dplyr)
library(Hmisc)
library(ggplot2)
library(raster)
library(sp)
library(RColorBrewer)

#change the name code of the species of interest
species_name<-"Ccapitata"

#select pseudo-absence generation method used
#either "3step" or "random'
bg<-"3step"  

#select variable selection method used
#either "PCA" or "randomForest"
reduction<-"randomForest"



#############################Computing ensemble  metrics ########################################

main_wd<-"E:/PHD/R_analysis/Ecoinformatics_workshop"
setwd(paste(main_wd,"/5_Ensemble",sep=""))
dir.create(species_name)
setwd(species_name)

#----------------------by vote------------------------------#
#probability > 0.5 is considered as presence

files_vote<-Sys.glob(paste(main_wd,"/4_MMA/",species_name,"/_prediction_maps/_predict_world/World_",species_name,"_*_data_vote1_",bg,"*",reduction,".txt",sep=""))
vote<-lapply(files_vote,function(x)read.table(x,h=T))
names(vote)<-c("CART","CTREE","KNN","LDA","LOG","NB","NNET","QDA","SVM")

vote_melt<-melt(vote,id.vars = c("xcoord", "ycoord","Vote","site"))
names(vote_melt)[5]<-"model"

sum_site<-aggregate(Vote~site,data=vote_melt,sum)

sum_site_coordinates<-full_join(sum_site,vote[[1]][,c(1:2,4)],by=c("site"))


#----------------------by weighted averaging +SD ------------------------------#

#read files with the scores for all the models

scores<-read.table(Sys.glob(paste(main_wd,"/4_MMA/",species_name,"/_results_compar.models/ranked_criteria_all_cv_",bg,"*",reduction,".txt",sep="")),h=T)


#compute somer's D using auc (Breiner et al 2015)
somers_function<-function(x){S<-2*(x$auc-0.5);as.data.frame(S)}
somers<-somers_function(scores)
somers<-cbind(somers,model=rownames(scores))



#check if there is a model that has a somer's value <0 (i.e AUC<0.5) 
#remove it if that is true

bad_model<-somers[somers$value<0,]
bad_model

#read file of probabilities
dir<-Sys.glob(paste(main_wd,"/4_MMA/",species_name,"/_prediction_maps/_predict_world/","World_*_*_threshold_0_cv_data_prob_",bg,"*",reduction,".txt",sep=""))
species_res<-lapply(dir,read.table,h=T)
names(species_res)<-c("CART","CTREE","KNN","LDA","LOG","NB","NNET","QDA","SVM")




#join all dataframes into one and add a column with model label
species_res_melt<-melt(species_res,id.vars = c("xcoord", "ycoord","Prob","site"))
names(species_res_melt)[5]<-"model"

#join somer's value
species_full<-full_join(species_res_melt,somers,by="model")
#calculate weight for averaging
weigth<-sum(somers$S)

#caculate weighted average/site
probxS<-mutate(species_full,weight_prob=(Prob*S))
sum_site<-aggregate(weight_prob~site,data=probxS,sum)
avg_site<-mutate(sum_site,average_prob=(weight_prob/weigth))
avg_site<-avg_site[,c("site","average_prob")]


#calculate weighted standard deviation/site
sd_wtd<-function(x){
            y<-subset(x,select=c("Prob","S"))
            res<-sqrt(wtd.var(y[,c("Prob")],y[,c("S")]))
            res                                                         
            }
  
sd_site<-ddply(species_full,.(site),sd_wtd)

#----------------------by median  ------------------------------#

median_site<-ddply(species_full,.(site),function(x) {y<-x[,c("Prob")];res<-median(y);res})


#----------------------by mean ------------------------------#

mean_site<-ddply(species_full,.(site),function(x) {y<-x[,c("Prob")];res<-mean(y);res})


#----------------------by best  ------------------------------#

best_model_auc<-rownames(scores[which.max(scores$auc),])

best<-species_full[species_full$model%in%best_model_auc,c("site","Prob")]


#----------------------by median PCA ------------------------------#


#prepare the data for pca
pca_data<-spread(species_res_melt,model,Prob)
pca_data<-pca_data[,4:12]

#perform a pca and get the top models (top half)
pca_res<-dudi.pca(pca_data,center=T,scannf=F,nf=2)
var_contribution<-pca_res$co^2##get variable contribution

#select the top half of the models that contribute the most to the axis 1
top_models<-names(tail(sort(var_contribution[,"Comp1"]),n=ceiling(dim(var_contribution)[1]/2)))           
top_pca<-species_full[species_full$model%in%top_models,]

#median between top models from pca

median_pca<-ddply(top_pca,.(site),function(x) {y<-x[,c("Prob")];res<-median(y);res})

#----------------------min predicted probability/site ------------------------------#
min_site<-ddply(species_full,.(site),function(x) {y<-x[,c("Prob")];res<-min(y);res})


#----------------------max predicted probability /site ------------------------------#
max_site<-ddply(species_full,.(site),function(x) {y<-x[,c("Prob")];res<-max(y);res})


#----------------------join all consensus approaches to coordinates ------------------------------#


final_res<-full_join(species_full[species_full$model=="CART",c("xcoord","ycoord","site")],avg_site,by="site")
final_res<-full_join(final_res,sd_site,by="site")
names(final_res)[5]<-"sd"

final_res<-full_join(final_res,median_site,by="site")
names(final_res)[6]<-"median"

final_res<-full_join(final_res,mean_site,by="site")
names(final_res)[7]<-"mean"

final_res<-full_join(final_res,best,by="site")
names(final_res)[8]<-"best"


final_res<-full_join(final_res,median_pca,by="site")
names(final_res)[9]<-"median_pca"


final_res<-full_join(final_res,sum_site_coordinates[,1:2],by="site")
names(final_res)[10]<-"vote"


final_res<-full_join(final_res,min_site,by="site")
names(final_res)[11]<-"min"


final_res<-full_join(final_res,max_site,by="site")
names(final_res)[12]<-"max"




#----------------------by zeta metric------------------------------------------------#
  #zeta is a metric developed by Golay 2013, it's the stadard deviation corrected by the mean probability 
#(coefficient of variation adapted for bouded data)


zeta_function<-function(sd,mean,max,min,nb_models){sd/(sqrt(nb_models/(nb_models-1))*sqrt(mean*(max-min)-mean^2))}
final_res<-mutate(final_res,zeta=zeta_function(sd,mean,max,min,9))


#write final table with all consensus results

write.table(final_res[,c(1:2,4:13)],paste("World_",species_name,"_consensus_results.txt",sep=""),row.names=F)


#----------------------hyperconsensus ------------------------------#


binary_final_res<-mutate(final_res,bin_average_prob=ifelse(average_prob>0.5,1,0),bin_median=ifelse(median>0.5,1,0),
                         bin_mean=ifelse(mean>0.5,1,0),bin_best=ifelse(best>0.5,1,0),bin_median_pca=ifelse(median_pca>0.5,1,0),
                         bin_vote=ifelse(vote>=5,1,0))

write.table(binary_final_res,paste("World_",species_name,"_binary_consensus_results.txt",sep=""),row.names=F)

hyper_consensus4<-apply(binary_final_res[,c("bin_average_prob","bin_median_pca","bin_best","bin_vote")],1,sum)
hyper_consensus4_final<-data.frame(binary_final_res[,c("xcoord","ycoord")],hyper_consensus4)

write.table(hyper_consensus4_final,paste("World_",species_name,"_hyper_consensus4_results.txt",sep=""),row.names=F)


############################# Mapping ########################################

final_res<-read.table(paste("NZ_",species_name,"_consensus_results.txt",sep=""),h=T)
hyper_consensus4_final<-read.table(paste("NZ_",species_name,"_hyper_consensus4_results.txt",sep=""),h=T)


##defining extent and creating an empty raster with that extent
xmin<-min(final_res[,c("xcoord")])
xmax<-max(final_res[,c("xcoord")])
ymin<-min(final_res[,c("ycoord")])
ymax<-max(final_res[,c("ycoord")])



##create raster with the extent defined previously
r <- raster(xmn=xmin,xmx=xmax,ymn=ymin,ymx=ymax)
##rasterize
raster_avg_prob<- rasterize(final_res[, 1:2], r, final_res[,c("average_prob")], fun=mean)
raster_sd<- rasterize(final_res[, 1:2], r, final_res[,c("sd")], fun=mean)
raster_median<- rasterize(final_res[, 1:2], r, final_res[,c("median")], fun=mean)
raster_vote<- rasterize(final_res[, 1:2], r, final_res[,c("vote")], fun=max)
raster_mean<- rasterize(final_res[, 1:2], r, final_res[,c("mean")], fun=mean)
raster_pca_median<- rasterize(final_res[, 1:2], r, final_res[,c("median_pca")], fun=mean)
raster_best<- rasterize(final_res[, 1:2], r, final_res[,c("best")], fun=mean)
raster_min<- rasterize(final_res[, 1:2], r, final_res[,c("min")], fun=mean)
raster_max<- rasterize(final_res[, 1:2], r, final_res[,c("max")], fun=mean)
raster_zeta<- rasterize(final_res[, 1:2], r, final_res[,c("zeta")], fun=mean)

##export rasters


writeRaster(raster_avg_prob,paste(species_name,"avg_prob.rst",sep=""),format="IDRISI",overwrite=T)
writeRaster(raster_median,paste(species_name,"median.rst",sep=""),format="IDRISI",overwrite=T)
writeRaster(raster_vote,paste(species_name,"vote.rst",sep=""),format="IDRISI",overwrite=T)
writeRaster(raster_pca_median,paste(species_name,"pca_median.rst",sep=""),format="IDRISI",overwrite=T)
writeRaster(raster_best,paste(species_name,"best.rst",sep=""),format="IDRISI",overwrite=T)
writeRaster(raster_mean,paste(species_name,"mean.rst",sep=""),format="IDRISI",overwrite=T)


##plot with ggplot, colours suitable for colour blind, numeric scale

#create break scale for legend
breaks_scale<-c(0.01,0.1,0.33,0.66,0.9,0.99)
breaks_scale_uniform<-seq(0,1,0.2)

#create reverse palette for vote plot
custom_YlGnBu=brewer.pal(name="YlGnBu", n=nlevels(as.factor(final_res$vote)))
if(nlevels(as.factor(final_res$vote))>9) custom_YlGnBu<-c(custom_YlGnBu,"#000000")
names(custom_YlGnBu)=rev(levels(as.factor(final_res$vote)))




png(paste(species_name,"_weighted_avg.png",sep=""),width=936,height= 455)
p1<-ggplot(data=final_res, aes(y=ycoord, x=xcoord)) +
  geom_raster(aes(fill=average_prob))+
  scale_fill_gradientn(name="Probability",colours=c("#0000FFFF","#FFFFFFFF","#666600"),limits=c(0,1),breaks=breaks_scale_uniform)+
  labs(title="Weighted averaged probability",x= "Longitude", y = "Latitude")+
  theme(panel.background = element_blank())
p1
dev.off()


png(paste(species_name,"_sd.png",sep=""),width=936,height= 455)
p2<-ggplot(data=final_res, aes(y=ycoord, x=xcoord)) +
  geom_raster(aes(fill=sd))+
  scale_fill_gradientn(name="Probability",colours=c("#CCCCCC","#999999","#000000"))+
  labs(title="Weighted standard deviation",x= "Longitude", y = "Latitude")+
  theme(panel.background = element_blank())
p2
dev.off()

png(paste(species_name,"_median_pca.png",sep=""),width=936,height= 455)
p3<-ggplot(data=final_res, aes(y=ycoord, x=xcoord)) +
  geom_raster(aes(fill=median_pca))+
  scale_fill_gradientn(name="Probability",colours=c("#0000FFFF","#FFFFFFFF","#666600"),limits=c(0,1),breaks=breaks_scale_uniform)+
  labs(title="Median-PCA probability",x= "Longitude", y = "Latitude")+
  theme(panel.background = element_blank())
p3
dev.off()

png(paste(species_name,"_vote.png",sep=""),width=936,height= 455)
p4<-ggplot(data=final_res, aes(y=ycoord, x=xcoord)) +
  geom_raster(aes(fill=as.factor(vote)))+
  scale_fill_manual(name="Vote",values=custom_YlGnBu)+
  labs(title="Number of model votes",x= "Longitude", y = "Latitude")+
  theme(panel.background = element_blank())
p4
dev.off()

png(paste(species_name,"_best.png",sep=""),width=936,height= 455)
p5<-ggplot(data=final_res, aes(y=ycoord, x=xcoord)) +
  geom_raster(aes(fill=best))+
  scale_fill_gradientn(name="Probability",colours=c("#0000FFFF","#FFFFFFFF","#666600"),limits=c(0,1),breaks=breaks_scale_uniform)+
  labs(title="Best model probability",x= "Longitude", y = "Latitude")+
  theme(panel.background = element_blank())
p5
dev.off()

png(paste(species_name,"_hyperconsensus.png",sep=""),width=936,height= 455)
p6<-ggplot(data=hyper_consensus4_final, aes(y=ycoord, x=xcoord)) +
  geom_raster(aes(fill=as.factor(hyper_consensus4)))+
  scale_fill_manual(name="Vote",values=custom_YlGnBu)+
  labs(title="Overlap of consensus approaches",x= "Longitude", y = "Latitude")+
  theme(panel.background = element_blank())
p6
dev.off()

png(paste(species_name,"_zeta.png",sep=""),width=936,height= 455)
p7<-ggplot(data=final_res, aes(y=ycoord, x=xcoord)) +
  geom_raster(aes(fill=zeta))+
  scale_fill_gradientn(name="Zeta",colours=c("#CCCCCC","#999999","#000000"))+
  labs(title="Zeta",x= "Longitude", y = "Latitude")+
  theme(panel.background = element_blank())
p7
dev.off()



png(paste(species_name,"_min.png",sep=""),width=936,height= 455)
p8<-ggplot(data=final_res, aes(y=ycoord, x=xcoord)) +
  geom_raster(aes(fill=min))+
  scale_fill_gradientn(name="Probability",colours=c("#0000FFFF","#FFFFFFFF","#666600"),limits=c(0,1),breaks=breaks_scale_uniform)+
  labs(title="Min",x= "Longitude", y = "Latitude")+
  theme(panel.background = element_blank())
p8
dev.off()


png(paste(species_name,"max.png",sep=""),width=936,height= 455)
p8<-ggplot(data=final_res, aes(y=ycoord, x=xcoord)) +
  geom_raster(aes(fill=max))+
  scale_fill_gradientn(name="Probability",colours=c("#0000FFFF","#FFFFFFFF","#666600"),limits=c(0,1),breaks=breaks_scale_uniform)+
  labs(title="Max",x= "Longitude", y = "Latitude")+
  theme(panel.background = element_blank())
p8
dev.off()



##plot for colour blind people and for exercice verbal scale
breaks_scale_verbal<-c(0.05,0.215,0.495,0.78,0.99)
labels_verbal<-c("Very unlikely","Unlikely","About as likely as not","Likely","Very likely")

png(paste(species_name,"_weighted_avg_verbalscale.png",sep=""),width=936,height= 455)
p9<-ggplot(data=final_res, aes(y=ycoord, x=xcoord)) +
  geom_raster(aes(fill=average_prob))+
  scale_fill_gradientn(name="Probability",colours=c("#0000FFFF","#FFFFFFFF","#666600"),limits=c(0,1),breaks=breaks_scale_verbal,labels=labels_verbal)+
  labs(title="Weighted averaged probability",x= "Longitude", y = "Latitude")+
  theme(panel.background = element_blank())
p9
dev.off()


