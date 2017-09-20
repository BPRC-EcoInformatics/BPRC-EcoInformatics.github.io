# Perform climate change adjustment

# Author: Tak Ikeda
# Date: Sep 19 2010

clim.change.adj <- function(dat,wc.var,v1,v2) {

#	# Return original data if variable not present, with warning
#	if (sum(colnames(dat)==wc.var)==0) {
#		dat <- dat
#		warning("Variable not in data!")
#	}

	# Adjustment if variable is present
#	if (sum(colnames(dat)==wc.var)!=0) {
		cc <- which(colnames(dat)==wc.var)
		var.adj <- rnorm(nrow(dat),v1,v2)
		var.new <- dat[,cc] + var.adj
		dat[,cc] <- var.new
#	}
	
return(dat)

}

