extract.range <- function(data, species.name,
			save.fig=FALSE,
			save.file=FALSE,
			color.pts="red",
			t1,t2,
			p=NULL)
{
####################################### 
#######################################
# MAKE SURE YOU LOAD WORLDCLIM FIRST !!
#######################################
#######################################

# extracts worldclim points from range of temperature, precipitation 
# for a species and returns subset, boolean vector, and figures.
# Input:
# - data = worldclim or pre-selected country
# - species.name = character, 2 or 3 letters, for species name, used for saving files
# - plot.fig, save.fig = logical
# - color.pts = color other than red
# - t1,2, p = c(min,max) of survival, vector of length 2. 
# Contents can be NULL to signify unknown upper/lower limit.

t_bool <- p_bool <- TRUE
#f <- function(x) sum(is.na(x))==0
#bool2 <- apply(worldclim[,c(59,60)],1,f)


#if (is.null(t1) & is.null(t2))   t_bool <- TRUE
#if (!is.null(t1) & is.null(t2))  t_bool <- data$cbio06 >= t1 #& !is.na(data$cbio06)
#if (is.null(t1) & !is.null(t2))  t_bool <- data$cbio05 <= t2 #& !is.na(data$cbio05)
#if (!is.null(t1) & !is.null(t2)) t_bool <- data$cbio06 >= t1 & data$cbio05 <= t2 #& !is.na(data$cbio05)
##if (!is.null(t1) & !is.null(t2)) t_bool <- data$cbio06 >= t1 & data$cbio05 <= t2 & data$cbio05 >= t1 & data$cbio06 <= t2 

#if (is.null(t1) & is.null(t2))   t_bool <- TRUE
#if (!is.null(t1) & is.null(t2))  t_bool <- data$cbio01 >= t1 #& !is.na(data$cbio01)
#if (is.null(t1) & !is.null(t2))  t_bool <- data$cbio01 <= t2 #& !is.na(data$cbio01)
#if (!is.null(t1) & !is.null(t2)) t_bool <- data$cbio01 >= t1 & data$cbio01 <= t2 #& !is.na(data$cbio05)

if (is.null(t1) & is.null(t2))   t_bool <- TRUE
if (!is.null(t1) & is.null(t2))  t_bool <- data$cbio01 >= t1 | data$cbio06 >= t1
if (is.null(t1) & !is.null(t2))  t_bool <- data$cbio01 <= t2 | data$cbio05 <= t2
if (!is.null(t1) & !is.null(t2)) t_bool <- data$cbio01 >= t1 & data$cbio01 <= t2 | data$cbio06 >= t1 & data$cbio05 <= t2


# if (is.null(p1) & is.null(p2))   p_bool <- TRUE
# if (!is.null(p1) & is.null(p2))  p_bool <- data$cbio14 > p1 #& !is.na(data$cbio14)
# if (is.null(p1) & !is.null(p2))  p_bool <- data$cbio13 < p2 #& !is.na(data$cbio13)
# if (!is.null(p1) & !is.null(p2)) p_bool <- data$cbio14 > p1 #& data$cbio13 < p2 & !is.na(data$cbio13)

subdata <- data[t_bool & p_bool,]

#if(save.file) write.table(subdata,paste(species.name,"_extract.txt",sep=""),row.names=F)
if(save.file) write.table(subdata,paste(species.name,"_extract_T",t1,"_",t2,".txt",sep=""),row.names=F)

plot(x.nz[,1],x.nz[,2], ylab="Latitude",xlab="Longitude",main=paste("Locations for Possible Survival, ",species.name,sep=""),col="honeydew", pch=15, cex=0.7)
#points(subdata$xcoord,subdata$ycoord,col=color.pts, pch=15, cex=0.7)
points(subdata$xcoord,subdata$ycoord,col=color.pts,pch=15,cex=0.2)

#points(worldclim[!bool2,2],worldclim[!bool2,3], col="black", cex=0.7, pch=15)
if(save.fig) savePlot(paste(species.name,"_extract_T",t1,"_",t2,".png",sep=""),type="png")

#return(subdata)
}

