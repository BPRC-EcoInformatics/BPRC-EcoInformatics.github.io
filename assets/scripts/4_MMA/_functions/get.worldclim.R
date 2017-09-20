# Assign best Worldclim database values to presence locations

# Author  :   Takayoshi Ikeda
# Updated :   June 2012

get.worldclim <- function(sp.name,wd,rdata="worldclim") {
 


# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# input:
# 	sp.name = species name, two or three characters
#	wd = work directory with the data
#	rdata = name of rdata file to use as background data
# 		= "worldclim" for original 19 bioclim + alt (default)
# output: 
#	- worldclim, alt, bio, tmin, tmax, prec data for presence locations of sp.name
#	- map of presence locations
#	- map of closest worldclim locations
#
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #
# - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - # - #


# ------------------------------------ #
# Coordinates of worldclim and species #
# ------------------------------------ #

	setwd(wd)
	dat <- read.table(paste(sp.name,"_data_no_replicates.txt",sep=""), header=T)
	cat("Read in species data \n")

# Initialize
	world <- get(rdata)
	xw <- world$xcoord
	yw <- world$ycoord
	xs <- dat[,1]
	ys <- dat[,2]

# --------------------- #
# Take regions of world #
# --------------------- #


# North America, 
	n.am <- c(-140,20,-55,70)

# South America, central America 
	s.am <- c(-110,-57,-35,20)

# North American islands
	n.am.i <- c(-180,0,-150,35)

# South American islands
	s.am.i <- c(-180,-50,-120,0)

# Europe
	eur <- c(-35,35,50,75)

# West Asia
	w.asia <- c(50,20,100,75)

# East Asia
	e.asia <- c(100,20,180,75)

# Southeast Asia
	se.asia <- c(60,-15,180,20)

# North Africa
	n.afr <- c(-35,0,50,35)

# South Africa
	s.afr <- c(-35,-40,60,0)

# Australasia
	aus <- c(110,-60,180,-15)

  
# Regions with worldclim and species data in them
	n.am.w <- xw>=n.am[1] & xw<=n.am[3] & yw>=n.am[2] & yw<=n.am[4]
	n.am.s <- xs>=n.am[1] & xs<=n.am[3] & ys>=n.am[2] & ys<=n.am[4]
	s.am.w <- xw>=s.am[1] & xw<=s.am[3] & yw>=s.am[2] & yw<=s.am[4]
	s.am.s <- xs>=s.am[1] & xs<=s.am[3] & ys>=s.am[2] & ys<=s.am[4]
	n.am.i.w <- xw>=n.am.i[1] & xw<=n.am.i[3] & yw>=n.am.i[2] & yw<=n.am.i[4]
	n.am.i.s <- xs>=n.am.i[1] & xs<=n.am.i[3] & ys>=n.am.i[2] & ys<=n.am.i[4]
	s.am.i.w <- xw>=s.am.i[1] & xw<=s.am.i[3] & yw>=s.am.i[2] & yw<=s.am.i[4]
	s.am.i.s <- xs>=s.am.i[1] & xs<=s.am.i[3] & ys>=s.am.i[2] & ys<=s.am.i[4]
	eur.w <- xw>=eur[1] & xw<=eur[3] & yw>=eur[2] & yw<=eur[4]
	eur.s <- xs>=eur[1] & xs<=eur[3] & ys>=eur[2] & ys<=eur[4]
	w.asia.w <- xw>=w.asia[1] & xw<=w.asia[3] & yw>=w.asia[2] & yw<=w.asia[4]
	w.asia.s <- xs>=w.asia[1] & xs<=w.asia[3] & ys>=w.asia[2] & ys<=w.asia[4]
	e.asia.w <- xw>=e.asia[1] & xw<=e.asia[3] & yw>=e.asia[2] & yw<=e.asia[4]
	e.asia.s <- xs>=e.asia[1] & xs<=e.asia[3] & ys>=e.asia[2] & ys<=e.asia[4]
	se.asia.w <- xw>=se.asia[1] & xw<=se.asia[3] & yw>=se.asia[2] & yw<=se.asia[4]
	se.asia.s <- xs>=se.asia[1] & xs<=se.asia[3] & ys>=se.asia[2] & ys<=se.asia[4]
	n.afr.w <- xw>=n.afr[1] & xw<=n.afr[3] & yw>=n.afr[2] & yw<=n.afr[4]
	n.afr.s <- xs>=n.afr[1] & xs<=n.afr[3] & ys>=n.afr[2] & ys<=n.afr[4]
	s.afr.w <- xw>=s.afr[1] & xw<=s.afr[3] & yw>=s.afr[2] & yw<=s.afr[4]
	s.afr.s <- xs>=s.afr[1] & xs<=s.afr[3] & ys>=s.afr[2] & ys<=s.afr[4]
	aus.w <- xw>=aus[1] & xw<=aus[3] & yw>=aus[2] & yw<=aus[4]
	aus.s <- xs>=aus[1] & xs<=aus[3] & ys>=aus[2] & ys<=aus[4]

  
# the rest of the world that was missed
  
  other.w <- !(n.am.w+s.am.w+n.am.i.w+s.am.i.w+eur.w+w.asia.w+e.asia.w+se.asia.w+n.afr.w+s.afr.w+aus.w)
  other.s <- !(n.am.s+s.am.s+n.am.i.s+s.am.i.s+eur.s+w.asia.s+e.asia.s+se.asia.s+n.afr.s+s.afr.s+aus.s)
  

  
cat("Regions set \n")

# ----------------------- #
# Closest worldclim point #
# ----------------------- #

# New worldclim datasets for regions
	world.n.am <- world[n.am.w,]
	world.s.am <- world[s.am.w,]
	world.n.am.i <- world[n.am.i.w,]
	world.s.am.i <- world[s.am.i.w,]
	world.eur <- world[eur.w,]
	world.w.asia <- world[w.asia.w,]
	world.e.asia <- world[e.asia.w,]
	world.se.asia <- world[se.asia.w,]
	world.n.afr <- world[n.afr.w,]
	world.s.afr <- world[s.afr.w,]
	world.aus <- world[aus.w,]
    world.other <- world[other.w,]
  
# Subsets of dat within regions
	dat.n.am <- dat[n.am.s,]
	dat.s.am <- dat[s.am.s,]
	dat.n.am.i <- dat[n.am.i.s,]
	dat.s.am.i <- dat[s.am.i.s,]
	dat.eur <- dat[eur.s,]
	dat.w.asia <- dat[w.asia.s,]
	dat.e.asia <- dat[e.asia.s,]
	dat.se.asia <- dat[se.asia.s,]
	dat.n.afr <- dat[n.afr.s,]
	dat.s.afr <- dat[s.afr.s,]
	dat.aus <- dat[aus.s,]
    dat.other <- dat[other.s,]

# Get closest by euclidean distance
	min.dist.n.am <- min.eucl.dist(dat.n.am,world.n.am)
	min.dist.s.am <- min.eucl.dist(dat.s.am,world.s.am)
	min.dist.n.am.i <- min.eucl.dist(dat.n.am.i,world.n.am.i)
	min.dist.s.am.i <- min.eucl.dist(dat.s.am.i,world.s.am.i)
	min.dist.eur <- min.eucl.dist(dat.eur,world.eur)
	min.dist.w.asia <- min.eucl.dist(dat.w.asia,world.w.asia)
	min.dist.e.asia <- min.eucl.dist(dat.e.asia,world.e.asia)
	min.dist.se.asia <- min.eucl.dist(dat.se.asia,world.se.asia)
	min.dist.n.afr <- min.eucl.dist(dat.n.afr,world.n.afr)
	min.dist.s.afr <- min.eucl.dist(dat.s.afr,world.s.afr)
	min.dist.aus <- min.eucl.dist(dat.aus,world.aus)
    min.dist.other <- min.eucl.dist(dat.other,world.other)
cat("min.dist done \n")

# ------------------------- #
# Plot best worldclim point #
# ------------------------- #
#	points(world$xcoord[as.numeric(min.dist.n.am)],world$ycoord[as.numeric(min.dist.n.am)],col=3,pch=16, cex=0.5)
#	points(world$xcoord[as.numeric(min.dist.s.am)],world$ycoord[as.numeric(min.dist.s.am)],col=3,pch=16, cex=0.5)
#	points(world$xcoord[as.numeric(min.dist.n.am.i)],world$ycoord[as.numeric(min.dist.n.am.i)],col=3,pch=16, cex=0.5)
#	points(world$xcoord[as.numeric(min.dist.s.am.i)],world$ycoord[as.numeric(min.dist.s.am.i)],col=3,pch=16, cex=0.5)
#	points(world$xcoord[as.numeric(min.dist.eur)],world$ycoord[as.numeric(min.dist.eur)],col=3,pch=16, cex=0.5)
#	points(world$xcoord[as.numeric(min.dist.w.asia)],world$ycoord[as.numeric(min.dist.w.asia)],col=3,pch=16, cex=0.5)
#	points(world$xcoord[as.numeric(min.dist.e.asia)],world$ycoord[as.numeric(min.dist.e.asia)],col=3,pch=16, cex=0.5)
#	points(world$xcoord[as.numeric(min.dist.se.asia)],world$ycoord[as.numeric(min.dist.se.asia)],col=3,pch=16, cex=0.5)
#	points(world$xcoord[as.numeric(min.dist.n.afr)],world$ycoord[as.numeric(min.dist.n.afr)],col=3,pch=16, cex=0.5)
#	points(world$xcoord[as.numeric(min.dist.s.afr)],world$ycoord[as.numeric(min.dist.s.afr)],col=3,pch=16, cex=0.5)
#	points(world$xcoord[as.numeric(min.dist.aus)],world$ycoord[as.numeric(min.dist.aus)],col=3,pch=16, cex=0.5)
#   points(world$xcoord[as.numeric(min.dist.other)],world$ycoord[as.numeric(min.dist.other)],col=3,pch=16, cex=0.5)
#	title("Closest Worldclim Locations")


# -------------------------------------------- #
# Match worldclim values to presence locations #
# -------------------------------------------- #

	species.index<- as.numeric(c(min.dist.n.am,min.dist.s.am,min.dist.n.am.i,min.dist.s.am.i,
					min.dist.eur,min.dist.w.asia,min.dist.e.asia,min.dist.se.asia,
					min.dist.n.afr,min.dist.s.afr,min.dist.aus,min.dist.other))
	species.index <- species.index[order(species.index)]

	
#	dat <- get(rdata)
#	species.dat <- dat[species.index,]
#	write.table(species.dat,paste(sp.name,"_",rdata,".txt",sep=""),row.names=F)

#	dat <- get(rdata)
	species.dat <- world[species.index,]

# ------------------------- #
# Plot best worldclim point #
# ------------------------- #
	win.graph(width=70,height=40)
	plot(border, cex=0.2, xlab="Longitude", ylab="Latitude")
    points(species.dat$xcoord,species.dat$ycoord,col=3,pch=16, cex=0.5)
	savePlot(paste(sp.name,"_closest_worldclim_locations.png",sep=""), type="png")
cat("plots done \n")


	write.table(species.dat,paste(sp.name,"_",rdata,".txt",sep=""),row.names=F)
  return(species.dat)
	
	cat("output done \n")
}
