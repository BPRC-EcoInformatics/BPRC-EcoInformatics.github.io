#### Script written by Ursula Torres
#### 13/10/15
#### Part of this script comes from Broennimann et al. 2012 (for niche overlap functions and occ.prep.functions)
#### Part of this script comes from Vicent Q. Vu. (2011) (for ggbiplot.R)


#######################################################################################################################################
# set species/directory

main_wd <-"E:/PHD/R_analysis/Ecoinformatics_workshop"

#select species
species_name='Ccapitata'

dir.create(paste(main_wd, "/2_NicheAnalysis/",species_name,sep=''))
setwd(paste(main_wd, "/2_NicheAnalysis/",species_name,sep=''))

#create a file where error messages will be recorded
zz <- file("all.Rout", open="wt")
sink(zz, type="message")

# Load functions and packages
#source("E:/PHD/R_analysis/Chapter1_niche_shift_MMA/niche_paper/_functions/plot.contrib.R")#Broennimann function
library(ecospat)
library(ggplot2)


#######################################################################################################################################
# Reorganized data set

# Load occurance, climate data and species range 
data_complete=read.table(paste(main_wd, "/1_CleaningData/",species_name,"/",species_name,"_clim_range_withoutRep.txt",sep=""),h=T,sep='\t')

# Select only complete data (which should be the case already)
data_complete<-data_complete[,c("xcoord","ycoord",paste("cbio",1:35,sep=""),"range")]

# Separate the native data from the non-native one
native=subset(data_complete,range=="Native",select=xcoord:range)
invasive=subset(data_complete,range=="Invasive",select=xcoord:range)

# Save the data into a text file
write.table(native,paste("native_",species_name,".txt",sep=""),row.names=F)
write.table(invasive,paste("invasive_",species_name,".txt",sep=""),row.names=F)


#######################################################################################################################################
# Select environmental background 
# choose the extent on which the environemntal background will calculated 
# it can be countries, continent or ecoregion
type_analysis='continent'

# If PROJ =F, the models are calibrated on both ranges.
# If PROJ =T, the models are calibrated on species 1 range only and projected to range 2. 
# Analyses where both ranges are needed (ex: LDA) are not done

PROJ = F

# Select climatic parameters for the niche analysis
#number_variable <- c(29,37)

#number of interation for the tests of equivalency and similarity
iterations<-1

#resolution of the gridding of the climate space
R<-100

#variables to be used (35 variables: column 37, 27 variables: column 29)
column_variable<-37




#######################################################################################################################################
source(paste(main_wd, "/2_NicheAnalysis/niche_analysis_function.R",sep=''))


niche_analysis_function(main_wd, species_name,native,invasive,type_analysis,column_variable,PROJ,iterations,R)
	


#close the file where error messages are recorded##############
sink(type="message") 
close(zz)

