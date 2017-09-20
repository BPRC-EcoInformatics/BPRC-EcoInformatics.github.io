niche_analysis_function2<- function(species_name,native,invasive,type_analysis,number_variable, clim1, clim2){

PROJ=F

colnames(clim1)<-c("x","y", paste("cbio",1:35,sep=""))
colnames(clim2)<-c("x","y", paste("cbio",1:35,sep=""))
clim12<-rbind(clim1,clim2)

# Occurence data for native and invasive range 
occ.sp1<-native[,1:37]
colnames(occ.sp1)[1:2]<-c("x","y")
occ.sp2<-invasive[,1:37]
colnames(occ.sp2)[1:2]<-c("x","y")

row.w.1.occ<-1-(nrow(occ.sp1)/nrow(rbind(occ.sp1,occ.sp2))) # prevalence of occ1 (weight by the number of occurences recorded)
row.w.2.occ<-1-(nrow(occ.sp2)/nrow(rbind(occ.sp1,occ.sp2))) # prevalence of occ2
row.w.occ<-c(rep(0, nrow(clim1)),rep(0, nrow(clim2)),rep(row.w.1.occ, nrow(occ.sp1)),rep(row.w.2.occ, nrow(occ.sp2)))

row.w.1.env<-1-(nrow(clim1)/nrow(clim12))  # prevalence of clim1
row.w.2.env<-1-(nrow(clim2)/nrow(clim12))  # prevalence of clim2
row.w.env<-c(rep(row.w.1.env, nrow(clim1)),rep(row.w.2.env, nrow(clim2)),rep(0, nrow(occ.sp1)),rep(0, nrow(occ.sp2)))

fac<-as.factor(c(rep(1, nrow(clim1)),rep(2, nrow(clim2)),rep(1, nrow(occ.sp1)),rep(2, nrow(occ.sp2))))

# global dataset for the analysis and rows for each sub dataset
clim1=as.data.frame(sapply(clim1,as.numeric))
clim2=as.data.frame(sapply(clim2,as.numeric))

data.env.occ<-rbind(clim1[,1:37],clim2[,1:37],occ.sp1,occ.sp2)[Xvar]
row.clim1<-1:nrow(clim1)
row.clim2<-(nrow(clim1)+1):(nrow(clim1)+nrow(clim2))
row.clim12<-1:(nrow(clim1)+nrow(clim2))
row.sp1<-(nrow(clim1)+nrow(clim2)+1):(nrow(clim1)+nrow(clim2)+nrow(occ.sp1))
row.sp2<-(nrow(clim1)+nrow(clim2)+nrow(occ.sp1)+1):(nrow(clim1)+nrow(clim2)+nrow(occ.sp1)+nrow(occ.sp2))

pca.cal <-dudi.pca(data.env.occ,row.w = row.w.env, center = T, scale = T, scannf = F, nf = 2)

th.sp<-0.025 # Quantile to use of the envrionemntal density used to delimit marginal climates
th.env<-"none"
z1<-ecospat.grid.clim.dyn(scores.clim12,scores.clim1,scores.sp1,R,th.sp=th.sp)
z2<-ecospat.grid.clim.dyn(scores.clim12,scores.clim2,scores.sp2,R,th.sp,th.sp)

index_25th<-ecospat.niche.dyn.index(z1,z2,intersection=0.25)

overlap<-ecospat.niche.overlap(z1,z2,cor=T)

D<-overlap$D
I<-overlap$I

result <- list(z1,z2,index_25th,D,I)

return{result}
}
