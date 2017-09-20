# Calculate minimum distance to coast 

# Author  :   Takayoshi Ikeda
# Updated :   October 2009

# obtain closest location with minimum euclidean distance

min.eucl.dist <- function(rd.data,world.data) {
	min.dist <- NULL
	if(is.null(nrow(rd.data))) 
		min.dist <- rownames(world.data[which.min(sqrt((world.data$xcoord-rd.data[1])^2+(world.data$ycoord-rd.data[2])^2)),])
	if(!is.null(nrow(rd.data))) {
		n <- nrow(rd.data)
		for (i in 1:n)
			min.dist[i] <- rownames(world.data[which.min(sqrt((world.data$xcoord-rd.data[i,1])^2+(world.data$ycoord-rd.data[i,2])^2)),])
		}
	min.dist
}
