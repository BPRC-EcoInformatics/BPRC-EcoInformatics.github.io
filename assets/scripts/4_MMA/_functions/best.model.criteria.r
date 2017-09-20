# Criteria for best model selection

best.model.criteria <- function(sp.name,method,save.file=T,object.assign=F,thresh=NA,prev=NA,sp=FALSE) {

### -------------------------------------------------- ###
### -------------------------------------------------- ###
# inputs:
# sp.name = character of species.name
# method = "boot", "cv", "holdout"
# save.file = T, to save in appropriate directory
# object.assign = F, do not return anything
#			if T, can assign to an object 

#outputs:
# Table for ranked criteria and best model based on all
# Table for ranked criteria and best model based on subset
### -------------------------------------------------- ###
### -------------------------------------------------- ###


# ---------------------------------------------------------------- #
# Set directory for species and load criteria from model comparison
# ---------------------------------------------------------------- #
 
if (is.na(thresh)) setwd(paste(main_wd,species.name,"_results_compar.models",sep="/"))  ####modified
if (!is.na(thresh)) setwd(paste(main_wd,species.name,"_results_compar.models,",thresh,prev,sep="/"))
ranks.table <- read.table(paste(sp.name,"_compar_models_",method,"_general_ranks_table_",bg,nb.var,"_",reduction,".txt",sep=""),header=T)
load(paste(sp.name,"_compar_models_",method,"_output_object_",bg,nb.var,"_",reduction,".Rdata",sep=""))

# ----------------------- #
# Extract criteria values
# ----------------------- #

indexes <- read.table(paste(sp.name,"_compar_models_",method,"_general_indexes_",bg,nb.var,"_",reduction,".txt",sep=""),header=T)
rank1 <- rank(-indexes[,1])
rank2 <- rank(-indexes[,2])
rank3 <- rank(-indexes[,3])
rank4 <- rank(-indexes[,4])
rank5 <- rank(-indexes[,5])
rank6 <- rank(-indexes[,6])
rank7 <- rank(-indexes[,7])
uncert <- read.table(paste(sp.name,"_compar_models_",method,"_general_uncertainty_",bg,nb.var,"_",reduction,".txt",sep=""),header=T)
rank8 <- rank(uncert[,1])
error <- read.table(paste(sp.name,"_compar_models_",method,"_general_error_",bg,nb.var,"_",reduction,".txt",sep=""),header=T)
rank9 <- rank(error[,1])
auc <- out$aucROC[rownames(indexes),]
rank10 <- rank(-auc)

ranks.table2 <- cbind(indexes[,1],rank1, 
				indexes[,2],rank2, 
				indexes[,3],rank3,
				indexes[,4],rank4, 
				indexes[,5],rank5, 
				indexes[,6],rank6,
				indexes[,7],rank7,
				uncert, rank8,
				error, rank9,
				auc, rank10)


ranks.table2 <- cbind(ranks.table2,apply(ranks.table2[,seq(2,ncol(ranks.table2),2)],1,sum))
rownames(ranks.table2) <- rownames(ranks.table)
colnames(ranks.table2)[seq(1,14,2)] <- colnames(indexes)
colnames(ranks.table2)[ncol(ranks.table2)] <- "score"


best.model <- c(rownames(ranks.table2)[ranks.table2$score==min(ranks.table2$score)])

# Use ALL criteria, but if there is more than one chosen model, 
# refer to subset and take a vote:
# subset = uncert, error, tss, roc auc

#if (length(best.model)>1) {
uncert <- read.table(paste(sp.name,"_compar_models_",method,"_general_uncertainty_",bg,nb.var,"_",reduction,".txt",sep=""),header=T)
rank8 <- rank(uncert[,1])
error <- read.table(paste(sp.name,"_compar_models_",method,"_general_error_",bg,nb.var,"_",reduction,".txt",sep=""),header=T)
rank9 <- rank(error[,1])
auc <- out$aucROC[rownames(indexes),]
rank10 <- rank(-auc)

	ranks.table3 <- cbind(indexes[,3],rank3,
					indexes[,5],rank5,
					indexes[,7],rank7,
					indexes[,6],rank6,
					uncert, rank8,
					error, rank9,
					auc, rank10)
	colnames(ranks.table3)[1] <- colnames(indexes)[3]
	colnames(ranks.table3)[3] <- colnames(indexes)[5]
	colnames(ranks.table3)[5] <- colnames(indexes)[7]
	colnames(ranks.table3)[7] <- colnames(indexes)[6]
	ranks.table3 <- cbind(ranks.table3,apply(ranks.table3[,seq(2,ncol(ranks.table3),2)],1,sum))
	rownames(ranks.table3) <- rownames(ranks.table)
	#colnames(ranks.table3)[1] <- c("TSS")
	colnames(ranks.table3)[ncol(ranks.table3)] <- "score"

	best.model2 <- c(rownames(ranks.table3)[ranks.table3$score==min(ranks.table3$score)])



if(sp){
indexes_sp <- read.table(paste(sp.name,"_compar_models_",method,"_general_indexes_sp_",bg,nb.var,"_",reduction,".txt",sep=""),header=T)
rank11 <- rank(-indexes[,1])
rank12 <- rank(-indexes[,2])
rank13 <- rank(-indexes[,3])
rank14 <- rank(-indexes[,4])
rank15 <- rank(-indexes[,5])
rank16 <- rank(-indexes[,6])


	ranks.table4 <- cbind(indexes_sp[,2],rank11,
					indexes[,3],rank13,
					indexes[,5],rank15,
					indexes[,4],rank14,
					uncert, rank8,
					error, rank9,
					indexes[,6],rank16)
	colnames(ranks.table4)[1] <- colnames(indexes_sp)[2]
	colnames(ranks.table4)[3] <- colnames(indexes_sp)[3]
	colnames(ranks.table4)[5] <- colnames(indexes_sp)[5]
	colnames(ranks.table4)[7] <- colnames(indexes_sp)[4]
	colnames(ranks.table4)[13] <- colnames(indexes_sp)[6]
	ranks.table4 <- cbind(ranks.table4,apply(ranks.table4[,seq(2,ncol(ranks.table4),2)],1,sum))
	rownames(ranks.table4) <- rownames(ranks.table)
	#colnames(ranks.table3)[1] <- c("TSS")
	colnames(ranks.table4)[ncol(ranks.table4)] <- "score"

	best.model3 <- c(rownames(ranks.table4)[ranks.table4$score==min(ranks.table4$score)])

}


if(sp){models.list <- c(best.model,best.model2,best.model3)}
else{models.list <- c(best.model,best.model2)}
unique.models <- unique(models.list)
vote <- NULL
for (i in 1:length(unique.models)) vote[i] <- sum(unique.models[i]==models.list)
best.model <- unique.models[which.max(vote)]

write.table(ranks.table3,paste("ranked_criteria_subset_",method,"_",bg,nb.var,"_",reduction,".txt",sep=""),row.names=T)
if(sp){write.table(ranks.table4,paste("ranked_criteria_subset_sp_",method,"_",bg,nb.var,"_",reduction,".txt",sep=""),row.names=T)}

#}

# Save in file
write.table(ranks.table2,paste("ranked_criteria_all_",method,"_",bg,nb.var,"_",reduction,".txt",sep=""),row.names=T)
write.table(best.model,paste("ranked_highest_model_",method,"_",bg,nb.var,"_",reduction,".txt",sep=""),row.names=F,col.names=F)

# if assign to object
if (object.assign) return(best.model)

}


