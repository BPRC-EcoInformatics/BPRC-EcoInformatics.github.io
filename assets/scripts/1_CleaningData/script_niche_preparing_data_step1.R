library(maptools)
library(rgdal)
library(sp)
library(raster)
library(rgeos)
library(RColorBrewer)
library(plyr)
library(reshape2)

#######################################################################################################################################
# set directory
species_name<-"Ccapitata"
main_wd <-"C:/Users/Ursula Torres/Documents/Ecoinformatics_workshop"
dir.create(paste(main_wd, "/1_CleaningData/",species_name,sep=''))
setwd(paste(main_wd, "/1_CleaningData/",species_name,sep=''))

######### Reading shapefile with biomes/ecoregions/realm/countries/continents
shape_file<-readShapePoly(paste(main_wd, "/0_Data/biomes/test_intersect.shp",sep=''))


#######################################################################################################################################
######## REMOVE SPATIAL REDUNDANCY
# Creating spatial points from coordinates and removing redundant coordinates
#coord_data1<-read.table(paste(main_wd,"/0_Data/SpeciesOccurences/Ldispar_GBIF_OrginalData_September2016.csv",sep=''),header=T,sep='\t')
coord_data1<-read.csv(paste(main_wd,"/0_Data/SpeciesOccurences/Ccapitata_occurencesGBIF_AndLitterature_September2016.csv",sep=''),header=T,sep='\t')
#coord_data1<-read.table(paste(main_wd,"/0_Data/SpeciesOccurences/Hvitripennis_occurenceDataGBIF_ISSG_CABI_October2016.csv",sep=''),header=T,sep='\t')
init_data <- coord_data1
dim(init_data)

if (species_name=='Ccapitata'){
coord_data1[,1] <- as.numeric(as.character(coord_data1[,1]))
coord_data1[,2] <-as.numeric(as.character(coord_data1[,2]))
}

############################## RESTART FROM HERE WHEN INDICATED 
# Remove occurences not geolocalised 
coord_data1 <- coord_data1[complete.cases(coord_data1[,c(1,2)]),]
dim(coord_data1)

#Search for records with geographic coordinates with same latitude, longitude
which(coord_data1[,1]==coord_data1[,2])
#Search for records with geographic coordinates 0,0
which(coord_data1[,1]==0)

#remove duplicates
coord=coord_data1[,1:2]
coordinates(coord)=c("xcoord","ycoord") #transform numbers in spatial data (coordinates)
ind_zero=zerodist(coord,zero=0.17,unique.ID=T) # identified redundant coordinate within a resolution of 10 deg
ind_tokeep=which(duplicated(ind_zero)==F)
coord_data1=coord_data1[ind_tokeep,]

coord_data<-coord_data1
coordinates(coord_data)= c("xcoord", "ycoord")
points = SpatialPointsDataFrame(coord_data,coord_data1)



#######################################################################################################################################
########  IDENTIFY ECOREGION, REALM, BIOME, COUNTRY AND CONTINENT FOR EACH OCCURENCE
# Overlay points and shape files 
o = over(points,shape_file[,c("ECO_NAME","WWF_REALM2","WWF_MHTNAM","COUNTRY","CONTINENT")])
points$ecoregion <- o$"ECO_NAME"
points$realm <- o$"WWF_REALM2"
points$biome <- o$"WWF_MHTNAM"
points$country <- o$"COUNTRY"
points$continent <- o$"CONTINENT"
species<-as.data.frame(points)

# Identifying points outside the terrestrial polygons  
species_NA1<-subset(species, is.na(species$country)|is.na(species$ecoregion))
species_NA<-subset(species, is.na(species$country)|is.na(species$ecoregion))
coordinates(species_NA)= c("xcoord", "ycoord")
points_NA = SpatialPointsDataFrame(species_NA,species_NA1) # coordinates of point not geolocalised within the terrestrial area
dim(species_NA1)


# To identified the ecoregion of coordinates identified as outside of terrestrial area:
# We search for the minimal distances with the closest terrestrial polygons in the units of the current projection
# and replace the ecoregion by the closest value
m_ecoregions = gDistance(points_NA, shape_file, byid=TRUE) 
# This function calculates the distances between the two geometries in the units of the current projection;
# It uses planar coordinates for euclidean distances;
# byid=TRUE calculate distance for every point
# Ignore warnings
col_ecoregions = apply(m_ecoregions, 2, function(x) which(x==min(x))) 
labels = unlist(shape_file[col_ecoregions,c("ECO_NAME","WWF_REALM2","WWF_MHTNAM","COUNTRY","CONTINENT")])
rows_to_change<-names(col_ecoregions)
labels<-as.data.frame(labels)
species[c(rows_to_change),c("ecoregion","realm","biome","country","continent")]<-labels

# document the coordinates used for definying the ecoregions
species[c(rows_to_change),c("xcoord.1","ycoord.1")]<-coordinates(shape_file[col_ecoregions,])
names(species)<-c("xcoord","ycoord","countrycode","ecoregion","realm","biome","country","continent","ecoregion_xcoord","ecoregion_ycoord")




# Check if countries of the geolocation as indicated in the original data set match the new geolocalisation
Countrycode_2Letter<-read.csv(paste(main_wd,"/0_Data/CountryCode_2Letter_3Letter.csv",sep=''),header=T,sep='\t')
index=NULL
for (i in 1:dim(species)[1]){ 
	if(is.na(species$countrycode[i])){ # if the name of the country is not documented check coordinates were the same
		if(species$xcoord[i] != species$ecoregion_xcoord[i] | species$ycoord[i] != species$ecoregion_ycoord[i]){index<-c(index,i)}
	}
	else{
		original_country <- as.character(species$countrycode[i])
		if(species_name=='Ldispar' | species_name=='Hvitripennis'){
     		# find the country name from it abbreviation
			countryIndex <- which(as.character(Countrycode_2Letter[,2])==as.character(species$countrycode[i]))
			original_country <-as.character(Countrycode_2Letter[countryIndex,1])
		}
		if(original_country!=as.character(species$country[i])){index<-c(index,i)}
	}
} 
species[index,c('xcoord','ycoord','countrycode','country','ecoregion_xcoord','ecoregion_ycoord')]

write.table(species[index,c('xcoord','ycoord','countrycode','country','ecoregion_xcoord','ecoregion_ycoord')],paste(species_name,'_checkIdentificationEcoregion.txt',sep=''),row.names=F, quote = FALSE, sep="\t")

# check location based on geoplaner http://www.geoplaner.com/
# For Ldispar 2 country were wrongly identify in the original data set -> use country variable and not countrycode!
# For the mdefly the data were sometimes wrongly labeled (change of country name, change of signs))



if(species_name=='Ccapitata'){
backup_species <- species
Corrected_coordinates<-read.csv('Ccapitata_checkIdentificationEcoregion_CorrectionByHand.csv',sep='\t',header=T)
Corrected_coordinates$countrycode <- as.character(Corrected_coordinates$countrycode)
Corrected_data <- species[,c('xcoord','ycoord','countrycode')]
Corrected_data$countrycode <- as.character(Corrected_data$countrycode)
for (i in seq(1,length(Corrected_coordinates$countrycode))){
	if (Corrected_coordinates$countrycode[i] == 'Eritea') Corrected_coordinates$countrycode[i]='Eritrea'
	if (Corrected_coordinates$countrycode[i] == 'Bostwana') Corrected_coordinates$countrycode[i]='Botswana'
	if (Corrected_coordinates$countrycode[i] == 'Lybia') Corrected_coordinates$countrycode[i]='Libya'
	if (Corrected_coordinates$countrycode[i] == 'Saudia Arabia') Corrected_coordinates$countrycode[i]='Saudi Arabia'
	if (Corrected_coordinates$countrycode[i] == 'Equatorial Guinea') Corrected_coordinates$countrycode[i]='Egypt' #mistake when dowloaded from GBIF!
	if (Corrected_coordinates$countrycode[i] == 'Yemen') Corrected_coordinates$countrycode[i]='B'
	if (Corrected_coordinates$countrycode[i] == 'Benin') Corrected_coordinates$countrycode[i]='Yemen'
	if (Corrected_coordinates$countrycode[i] == 'B') Corrected_coordinates$countrycode[i]='Benin'
	if (Corrected_coordinates$countrycode[i] == 'Congo DRC') Corrected_coordinates$countrycode[i]='Congo'
	if (Corrected_coordinates$countrycode[i] == 'Sao Tome') Corrected_coordinates$countrycode[i]='Sao Tome and Principe'
	
}
temp_data <-  Corrected_coordinates[,c("Corrected_xcoord","Corrected_ycoord","countrycode")]
names(temp_data) <- c('xcoord','ycoord','coutrycode')
Corrected_data[index,] <- temp_data
Corrected_data$countrycode <- as.factor(Corrected_data$countrycode)

#Start again with the new data set!
coord_data1 <- Corrected_data
}

 


#remove duplicates
#coord=coord_data1[,1:2]
#coordinates(coord)=c("xcoord","ycoord") #transform numbers in spatial data (coordinates)
#ind_zero=zerodist(coord,zero=0.17,unique.ID=T) # identified redundant coordinate within a resolution of 10 deg
#ind_tokeep=which(duplicated(ind_zero)==F)
#coord_data1=coord_data1[ind_tokeep,]

#dim(coord_data1)

if (species_name=='Hvitripennis'){
Locality <- coord_data1[,4]
coord_data1 <- coord_data1[,1:3]
}

#######################################################################################################################################
######## DOCUMENT THE RANGE OF THE SPECIES
range=NULL
for (i in 1:dim(species)[1]){
	if(species_name=='Ccapitata'){
		if(species$continent[i]=="Africa"){range[i]<-"Native"} # Africa for the Medfly and America for the gyspsy moth
		else {range[i]<-"Invasive"}
	} 
	if(species_name=='Ldispar'){
		if(species$continent[i]=="Europe"){range[i]<-"Native"} # Africa for the Medfly and America for the gyspsy moth
		else {range[i]<-"Invasive"}
	} 
	if(species_name=='Hvitripennis'){
		native_range <- c('Mexico','Alabama','Arkansas','Florida','Georgia','Louisiana', 'Mississippi','North Carolina','South Carolina','Texas')
		if(Locality[i] %in% native_range){range[i]<-"Native"}
		else {range[i]<-"Invasive"} 	
	}
 
}
species$range<-as.factor(range)

#cbind(as.character(species$range), as.character(Locality))

write.table(species,paste(species_name,'_range_ecoregionDocumented.txt',sep=''),row.names=F, quote = FALSE, sep="\t")
write.table(species[c(rows_to_change),],paste(species_name,'_range_ecoregionDocumented_IssueWithCoordinates.txt',sep=''),row.names=F, quote = FALSE, sep="\t")





#######################################################################################################################################
######## EXTRACT CLIMATIC VARIABLES FOR OCCURENCE POINTS AND THREE TYPES OF BACKGROUND

### Read climond database - deg 10
dir<-Sys.glob(paste(main_wd,"/0_Data/Climond_CM10_1975H_Bio_ESRI_V1.2/CM10_1975H_Bio_V1.2/bio*/hdr.adf",sep=''))
# Select the climatic variables to use
dir<-dir[1:35]

bio<-lapply(dir,raster) # transform data into raster
climond<-stack(bio) # concatenates multiple vectors contained in bio into a single vector along with a factor indicating where each observation originated.

### Extract climatic variables for presence points
dist_clim<-extract(climond,species[,c("xcoord","ycoord")],df=T)
dist_clim<-data.frame(dist_clim,species[,c("xcoord","ycoord","range","ecoregion","realm","biome","country","continent")])


#######################################################################################################################################
######## Dealing with  with data uncorrectly geolocalised

# There are several ways for dealing with data uncorrectly geolocalised
# 1. Remove the data from the analysis
# 2. Search the closest terrestrial polygons using a buffer search zone in the units of the current projection;

dist_clim_init <- dist_clim

###### 1 . Remove the data from the analysis
length(which(is.na(dist_clim[,2])))
# 21 points were found to not have their climate documented for L.dispar
# 30 points were found to not have their climate documented for C.Capitata
# WARNINGS: the border of the terrestrial shapefile and of climond points are slightly different
# Therefore points identified as not correctly localised here can be different than in the step before

dist_clim_na.omit<-dist_clim_init[complete.cases(dist_clim_init[,1:36]),]


###### 2. Search the closest terrestrial polygons using a buffer search zone in the units of the current projection;
# Identify NA values (point not properly geolocalised))
dist_clim_NA<-dist_clim[is.na(dist_clim$hdr.1),]
# Search for the closest terrestrial polygons using a buffer of 10km (resolution of one cell)
points_NA_dist<-extract(climond,dist_clim_NA[,c("xcoord","ycoord")],buffer=10000)
## arrange data
list_columns_points_NA_dist<-lapply(points_NA_dist,function(x)as.list(as.data.frame(x)))
points_NA_dist_merged<-melt(list_columns_points_NA_dist)
points_NA_dist_clim<-aggregate(value~L1+L2,data=points_NA_dist_merged,mean)
points_NA_dist_clim<-dcast(points_NA_dist_clim,L1~L2)
points_NA_dist_clim<-points_NA_dist_clim[,c("L1",names(dist_clim)[2:36])]

## indices to replace
ind_to_replace<-which(is.na(dist_clim$hdr.1))
ind_to_replace_order<-ind_to_replace[as.numeric(points_NA_dist_clim$L1)]

## replace NA values with mean values
dist_clim[ind_to_replace_order,c(2:36)]<-points_NA_dist_clim[,c(2:36)]
## remove NA values that even with buffer could not been solved
dist_clim<-na.omit(dist_clim)


##change names columns
colnames(dist_clim)[2:36]<-paste("cbio",01:35,sep="")
write.table(dist_clim[,c(2:dim(dist_clim)[2])],paste(species_name,"_clim_range_withoutRep.txt",sep=''),row.names=F, quote = FALSE, sep="\t")



#######################################################################################################################################
#####extract climate variables for background
#####list of polygons that will be used for cutting the data 

invasive_continent<-unlist(unique(subset(dist_clim,dist_clim$range=="Invasive",select=continent)))
native_continent<-unlist(unique(subset(dist_clim,dist_clim$range=="Native",select=continent)))

invasive_countries<-unlist(unique(subset(dist_clim,dist_clim$range=="Invasive",select=country)))
native_countries<-unlist(unique(subset(dist_clim,dist_clim$range=="Native",select=country)))

invasive_ecoregion<-unlist(unique(subset(dist_clim,dist_clim$range=="Invasive",select=ecoregion)))
native_ecoregion<-unlist(unique(subset(dist_clim,dist_clim$range=="Native",select=ecoregion)))


#subset the polygons using the list previously generated
shape_countries_invasive<-shape_file[shape_file$COUNTRY %in% invasive_countries & shape_file$CONTINENT %in% invasive_continent,]
shape_countries_native<-shape_file[shape_file$COUNTRY %in% native_countries & shape_file$CONTINENT %in% native_continent,]

shape_continent_invasive<-shape_file[shape_file$CONTINENT %in% invasive_continent,]
shape_continent_native<-shape_file[shape_file$CONTINENT %in% native_continent,]


shape_ecoregion_invasive<-shape_file[shape_file$ECO_NAME %in% invasive_ecoregion & shape_file$CONTINENT %in% invasive_continent,]
shape_ecoregion_native<-shape_file[shape_file$ECO_NAME %in% native_ecoregion&shape_file$CONTINENT %in% native_continent,]



#extract by mask raster values
mask_countries_invasive<-mask(climond,shape_countries_invasive)
mask_countries_invasive_df<-na.exclude(as.data.frame(mask_countries_invasive,xy=T))
mask_countries_native<-mask(climond,shape_countries_native)
mask_countries_native_df<-na.exclude(as.data.frame(mask_countries_native,xy=T))

mask_continent_invasive<-mask(climond,shape_continent_invasive)
mask_continent_invasive_df<-na.exclude(as.data.frame(mask_continent_invasive,xy=T))
mask_continent_native<-mask(climond,shape_continent_native)
mask_continent_native_df<-na.exclude(as.data.frame(mask_continent_native,xy=T))


mask_ecoregion_invasive<-mask(climond,shape_ecoregion_invasive)
mask_ecoregion_invasive_df<-na.exclude(as.data.frame(mask_ecoregion_invasive,xy=T))
mask_ecoregion_native<-mask(climond,shape_ecoregion_native)
mask_ecoregion_native_df<-na.exclude(as.data.frame(mask_ecoregion_native,xy=T))



write.table(mask_countries_invasive_df,paste("mask_",species_name,"_countries_invasive_df.txt",sep=""),row.names=F)
write.table(mask_countries_native_df,paste("mask_",species_name,"_countries_native_df.txt",sep=""),row.names=F)

write.table(mask_continent_invasive_df,paste("mask_",species_name,"_continent_invasive_df.txt",sep=""),row.names=F)
write.table(mask_continent_native_df,paste("mask_",species_name,"_continent_native_df.txt",sep=""),row.names=F)

write.table(mask_ecoregion_invasive_df,paste("mask_",species_name,"_ecoregion_invasive_df.txt",sep=""),row.names=F)
write.table(mask_ecoregion_native_df,paste("mask_",species_name,"_ecoregion_native_df.txt",sep=""),row.names=F)

