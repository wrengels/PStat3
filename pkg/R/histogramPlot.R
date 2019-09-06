
histoplot <- function(histo, obstat=5, tailstat=15, alph=.8, testname="The Test Statistic") {
	requireNamespace("ggplot2")
	require("ggplot2")
	hd <- data.frame(x=histo$mids, y=histo$density)
	if(is.numeric(histo$statObs)) obstat <- histo$statObs
	if(is.numeric(histo$statTail)) tailstat <- histo$statTail
	if(!is.null(histo$testname)) testname <- histo$testname
	if("asymp" %in% names(histo)) hd$asymp <- histo$asymp
	hd$co <- alpha("black", alph)
	if(obstat<tailstat) {
		leftob <- obstat
		leftcolor <- "red"
		rightob <- tailstat
		rightcolor <- "blue"
	} else {
		leftob <- tailstat
		leftcolor <- "blue"
		rightob <- obstat
		rightcolor <- "red"
	}
	hd$co <- ifelse(hd$x < leftob, ggplot2::alpha(leftcolor, alph), hd$co)
	hd$co <- ifelse(hd$x > rightob, ggplot2::alpha(rightcolor, alph), hd$co)
	lwd <- 200/nrow(hd)
	g <- ggplot2::ggplot(data=hd, aes(x=x, y=y)) + geom_segment(aes(xend=x, yend=0, color=co), size=lwd)
	g <- g + theme_classic()
	g <- g + theme(axis.line.x=element_line(size=0.6), axis.line.y=element_line(size=0.6))
	g <- g + scale_color_identity()
	g <- g + xlab(testname) + ylab("Probability")
	if(length(hd$asymp)) {g <- g + geom_line(aes(y=hd$asymp), color="darkorange", size=1.2)}
	g
}

tstat <- function(x, y) {
	xy <- c(x,y)
	sumxy <- sum(xy)
	sumx <- sum(x)
	sumy <- sumxy - sumx
	nx <- length(x)
	ny <- length(y)
	sumxx <- sum(x*x)
	sumxxyy <- sum(xy*xy)
	sumyy <- sumxxyy - sumxx
	vx <- sumxx/nx - (sumx/nx)^2
	vy <- sumyy/ny - (sumy/ny)^2
	t1 <- sumx/nx - sumy/ny
	t <- t1/sqrt(vx/(nx-1) + vy/(ny-1))
	
	t
}

dift <- function(xy,xset) {
#	xy <- c(x,y)
	x <- xy[xset]
#	print(x)
	y <- xy[-xset]
	sumxy <- sum(xy)
	sumx <- sum(x)
	sumy <- sumxy - sumx
	nx <- length(x)
	ny <- length(y)
	sumxx <- sum(x*x)
	sumxxyy <- sum(xy*xy)
	sumyy <- sumxxyy - sumxx
	vx <- sumxx/nx - (sumx/nx)^2
	vy <- sumyy/ny - (sumy/ny)^2
	t1 <- sumx/nx - sumy/ny
	t <- t1/sqrt(vx/(nx-1) + vy/(ny-1))
	c(t1,t)
}

tbatch <- function(x,y, B=1000, histobins = 500) {
	xy <- c(x,y)
	nx <- length(x)
	ny <- length(y)
	nxm1 <- nx - 1
	nym1 <- ny - 1
	nxy <- nx + ny
	xxyy <- xy^2
	sumxy <- sum(xy)
	sumxxyy <- sum(xxyy)
	
	tsta <- function(s) {
		sumx <- sum(xy[s])
		sumy <- sumxy - sumx
		sumxx <- sum(xxyy[s])
		sumyy <- sumxxyy - sumxx
		vx <- sumxx/nx - (sumx/nx)^2
		vy <- sumyy/ny - (sumy/ny)^2
		return((sumx/nx - sumy/ny)/sqrt(vx/nxm1 + vy/nym1))
	}
	cat("\nt = ",obstat <- tsta(1:nx), "\n")
	tb <- replicate(B, tsta(sample(nxy, nx)))
	h <- hist(tb, breaks=histobins, plot=F)
	
	# Degrees of freedom if unequal variances assumed
	vx <- var(x); vy <- var(y)
	dfnumerator <- (vx/nx + vy/ny)^2
	dfdenom <- (vx/nx)^2/nxm1 + (vy/ny)^2/nym1
	df <- dfnumerator/dfdenom
	cat("    df =", df, "\n" )
	
	h$asymp <- dt(h$mids, df=df)
	g <- histoplot(h, obstat=obstat, tailstat=-obstat, testname="t Statistic")
	g

}

tbatch.ev <- function(x, y, B=1000, histobins=500) {
	# t with equal variance is ALMOST LINEAR with x sum:   t = ~(const)(xbar-ybar)
	
}