# lwd<-col<-lty<-1
# endpt.vec <- c(60,10,100,40)


seg.square <- function(endpt.vec,lwd=1,col=1,lty=1) {
	# to draw a box in a plot given bottom left and top right corner coordinates
	# col, lwd, lty can be specified
	segments(endpt.vec[1],endpt.vec[2],endpt.vec[3],endpt.vec[2],lwd=lwd,col=col,lty=lty);
	segments(endpt.vec[1],endpt.vec[2],endpt.vec[1],endpt.vec[4],lwd=lwd,col=col,lty=lty)
	segments(endpt.vec[3],endpt.vec[2],endpt.vec[3],endpt.vec[4],lwd=lwd,col=col,lty=lty)
	segments(endpt.vec[1],endpt.vec[4],endpt.vec[3],endpt.vec[4],lwd=lwd,col=col,lty=lty)
}

