# Program to calculate performance measures for desired number of random 
# sets of pseudo-absences. Exports the performance measures in separate files.

perf.meas.rand.pseudo <- function(nrep.range, 
				species.name=species.name, 
				thresh=thresh, 
				prev=prev,
				b.k,b.s,b.n,
				c.k,c.s,c.n) {

time1 <- proc.time()

for (i in nrep.range) {
cat("rep for random pseudo-absences, i =",i,"\n")
	randompseudoabsence(sp.name=species.name,thresh=0,preva=prev,random.nb=i)
	dat.random.abs <- read.table(paste(species.name,"_data_pres_random_abs_set",i,"_threshold_",thresh,"_prev_",prev,".txt",sep=""),header=T)
	dat_best <- dat.random.abs[,c(as.character(c(vs_best)$x),"pres")]
	write.table(dat_best,paste(species.name,"_data_pres_random_abs_best_set",i,"_threshold_",thresh,"_prev_",prev,".txt",sep=""),row.names=F)


# Import data
	setwd(paste("C:\\MMA\\",species.name,"\\_dist_data_OCSVM",sep=""))
	data <- read.table(paste(species.name,"_data_pres_random_abs_best_set",i,"_threshold_",thresh,"_prev_",prev,".txt",sep=""), header=TRUE)
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
	setwd(paste("C:\\MMA\\",species.name,"\\_results_compar.models",sep=""))
	cv.k <- 10
	if (prevalence < 50) cv.k <- 5
	knrep <- 200/cv.k



#----------------------#
#  Same Parameter settings  |
#----------------------#

# boot:
setwd(paste("C:\\MMA\\",species.name,"\\_results_compar.models",sep=""))
b.k1 <- b.k
b.s1 <- b.s[1]
b.s2 <- b.s[2]
b.n1 <- b.n[1]
b.n2 <- b.n[2]
b.n3 <- b.n[3]

# cv:
c.k1 <- c.k
c.s1 <- c.s[1]
c.s2 <- c.s[2]
c.n1 <- c.n[1]
c.n2 <- c.n[2]
c.n3 <- c.n[3]

#-------------------------#
#  Fix Weights for nnet   |
#-------------------------#


setwd(paste("C:\\MMA\\",species.name,"\\_results_compar.models",sep=""))
nnet.Wts.boot <- read.table(paste(species.name,"_best_weights_",b.n1,"_",b.n2,"_",b.n3,".txt",sep=""),header=T)
nnet.Wts.cv <- nnet.Wts.boot
nnet.Wts.boot <- c(nnet.Wts.boot$x)
nnet.Wts.cv <- c(nnet.Wts.cv$x)




#--------------#
#  COMPARISON  | SOMETHING NOT WORKING WITH 0.632+ERROR !!!err.mat dimension... and 3 NAs...
#--------------#

## Comparison of the performance of supervised models using different measures

setwd(paste("C:\\MMA\\",species.name,"\\_results_compar.models",sep=""))

# CROSS-VALIDATION
win.graph()
result.compar.cv <- compar.models(x=data.x,
                               y=data.y, 
                               methods=c("LDA","QDA","LOG","NB","CART","CTREE","KNN","SVM","NNET"),
                               knn.par = b.k1, 
                               svm.par = c("rbfdot",b.s1,b.s2),                                
                               nnet.par = c(b.n1,b.n2,b.n3),
                               Wts = nnet.Wts.boot,
                               resampling="cv",
                               nrep=knrep,
                               train.frac=0.75,
                               k=cv.k,
                               file = T, 
                               plots=T, 
                               plots.all=T, 
                               boxplots=T, 
                               bagging=NULL,
                               boosting=NULL)





# --------------------------------------------------- #
# Script for best model selection based on multiple criteria #
# --------------------------------------------------- #

source("C:\\MMA\\_functions\\best.model.criteria.r")

# --------------------------------------------------- #
# Best models & ranked criteria tables saved into files
# --------------------------------------------------- #

best.model.criteria(species.name,"cv",thresh=NA,prev=60)



#spm1 <- read.table("ranked_criteria_all_cv.txt",skip=1,header=F)[8,c(18,14,12,20,22)]
#svm.performance.measure[i,] <- c(t(spm1))
spm <- read.table("ranked_criteria_all_cv.txt",skip=1,header=F)[,c(2*(1:10),22)]
lda.performance.measure.all[i,]   <- c(t(spm[1,]))
qda.performance.measure.all[i,]   <- c(t(spm[2,]))
log.performance.measure.all[i,]   <- c(t(spm[3,]))
nb.performance.measure.all[i,]    <- c(t(spm[4,]))
cart.performance.measure.all[i,]  <- c(t(spm[5,]))
ctree.performance.measure.all[i,] <- c(t(spm[6,]))
knn.performance.measure.all[i,]   <- c(t(spm[7,]))
svm.performance.measure.all[i,]   <- c(t(spm[8,]))
nnet.performance.measure.all[i,]  <- c(t(spm[9,]))




graphics.off()

# do for other models too
write.table(lda.performance.measure.all,paste(species.name,"_performance_measure_all_lda_threshold_",thresh,"_prev_",prev,".txt",sep=""),row.names=F)
write.table(qda.performance.measure.all,paste(species.name,"_performance_measure_all_qda_threshold_",thresh,"_prev_",prev,".txt",sep=""),row.names=F)
write.table(log.performance.measure.all,paste(species.name,"_performance_measure_all_log_threshold_",thresh,"_prev_",prev,".txt",sep=""),row.names=F)
write.table(nb.performance.measure.all,paste(species.name,"_performance_measure_all_nb_threshold_",thresh,"_prev_",prev,".txt",sep=""),row.names=F)
write.table(knn.performance.measure.all,paste(species.name,"_performance_measure_all_knn_threshold_",thresh,"_prev_",prev,".txt",sep=""),row.names=F)
write.table(cart.performance.measure.all,paste(species.name,"_performance_measure_all_cart_threshold_",thresh,"_prev_",prev,".txt",sep=""),row.names=F)
write.table(ctree.performance.measure.all,paste(species.name,"_performance_measure_all_ctree_threshold_",thresh,"_prev_",prev,".txt",sep=""),row.names=F)
write.table(svm.performance.measure.all,paste(species.name,"_performance_measure_all_svm_threshold_",thresh,"_prev_",prev,".txt",sep=""),row.names=F)
write.table(nnet.performance.measure.all,paste(species.name,"_performance_measure_all_nnet_threshold_",thresh,"_prev_",prev,".txt",sep=""),row.names=F)


}

out <- list()
out$lda <- lda.performance.measure.all
out$qda <- qda.performance.measure.all
out$log <- log.performance.measure.all
out$nb <- nb.performance.measure.all
out$knn <- knn.performance.measure.all
out$cart <- cart.performance.measure.all
out$ctree <- ctree.performance.measure.all
out$svm <- svm.performance.measure.all
out$nnet <- nnet.performance.measure.all

return(out)
save(out, file = paste(species.name,"_performance_measure_all_threshold_",thresh,"_prev_",prev,".Rdata",sep=""))
cat("for", length(nrep), "reps,", (proc.time() - time1)[3]/60, "minutes \n")

}

