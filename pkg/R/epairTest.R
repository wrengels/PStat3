# Exact test for comparing magnitude of values in sets of pairs
# (c) William R. Engels, 2019

#' 
#' 


#' @title
#' epairtest
#' @description
#' Compare x versus y when values are paired. 
#' Uses efficient algorithms for either complete enumeration or monte carlo test. 
#' Resulting distribution is displayed graphically. 
#' Test statistic is sum of x values or the mean of the x values
#' Trials consist of pairings

#' 
#' @details
#' 


#' @param x array of x values to compare with y
#' @param y array with same length as x. y[1] is compared with x[1] etc.
#' @param stat currently ignored
#' @param histobins number of bins for distribution
#' @param histoLeft left end of plotted distribution. Leave as NA to have it set automatically.
#' @param histoRight right end of distribution. NA means set automatically
#' @param safeSecs Stop the calculations after this many seconds
#' @param safeSets Do not attempt full enumeration if there are more than this many subsets
#' @param plot Display distribution histogram if TRUE.
#' @param detail How much detail to print after completion of test (currently 0 or 1)
#' @param rseed supplements current time to seed random numbers. Can be ignored for full enumeration or if no more than 1 run per second.
#' @param B The number of random trials to use. Set to zero to use full enumeration.


#' @return Returns a list which includes the following
#' \item{x}{vector of differences (\code{x-y})}
#' \item{n}{the number of pairs = length of x or y}
#' \item{nTodo}{the number of monte carlo trials requested. Not used in full enumeration}
#' \item{statLeft}{mimimum possible value of the sum of differences statistic}
#' \item{statRight}{maximum possible sum of differences. \code{statRight = -statLeft}}
#' \item{histobins}{number of bins used in histogram}
#' \item{histoLeft, histoRight}{boundaries of plotted histogram}
#' \item{safeSecs}{process aborted after this many seconds}
#' \item{nlower, nupper}{the numbers of trials falling into the left and right tails}
#' \item{histoData}{found numbers of trials in each bin of histogram}
#' \item{subsetsExamined}{The actual number of combinations examined, as opposed to the requested number}
#' \item{timespent}{the time in seconds spent in the C routine}
#' \item{statObs}{observed sum of differences}
#' \item{statTail}{negative of \code{statObs}}
#' \item{P2t}{two-tail P value}
#' \item{P1t}{one-tail P value. The tail is always chosed as the one in the same direction as observed}
#' \item{testname}{"paired:    mean(x-y)"}



#' 
#' 
#' @examples
#' x <- rnorm(10)
#' y <- rnorm(10)
#' result <- e.pairtest(x, y) # performs full enumeration and plots results
#' result <- e.pairtest(x, y, B=1e8) # same test, but uses 100 million random trials instead of full enumeration
#' 


#' @useDynLib PStat3
#' @export
e.pairtest <- function(x, y, stat=c("diff"), histobins=500, histoLeft=NA, histoRight=NA, safeSecs=30, safeSets=1e10, plot=T, detail=2, rseed=0, B=0) {
	n <- length(y)
	stat <- match.arg(stat)
	if(length(x) != n) stop("Values must be pairs")
	totalSubsets <- 2^n
	if(totalSubsets > safeSets & B==0){
		cat("\n")
		warning("\nThere are ", totalSubsets, " subsets to examine. The current limit is ", safeSets, ". Change parameter \"safeSets\" for different limit.\nSwitching to random subset method.")
		B <- 0 # 1e5
	} 
	z <- (x - y)/n;
	statObs <- sum(z)
	statTail <- -statObs
	difs <- abs(z)
	if(is.na(histoLeft)){
		histoLeft <- -sum(difs)
	}
	if(is.na(histoRight)){
		histoRight <- sum(difs)
	}
	nlower <- 0
	nupper <- 0
	rhistoData <- 1:histobins
	subsetsExamined <- -42
	timespent <- 0.0
	rhistoData[5] <- rseed
	almost1 <- 0.999999
	
	method <- ifelse(B==0, "epairTest", "rpairTest")
	b <- .C(method,
		x=as.double(difs),
		n=as.integer(n),
		nTodo=as.double(B),
		statLeft=as.double(-abs(statObs) * almost1),
		statRight=as.double(abs(statObs) * almost1),
		histobins=as.integer(histobins),
		histoLeft=as.double(histoLeft),
		histoRight=as.double(histoRight),
		safeSecs=as.integer(safeSecs),
		nlower=as.double(nlower),
		nupper=as.double(nupper),
		histoData=as.double(rhistoData),
		subsetsExamined=as.double(subsetsExamined),
		timespent=as.double(timespent)
		)
	
	b$statObs <- statObs;
	b$statTail <- statTail;
	b$P2t <- (b$nupper + b$nlower)/b$subsetsExamined
	b$P1t <- ifelse(statObs < 0, b$nlower/b$subsetsExamined, b$nupper/b$subsetsExamined)
	b$testname <- "paired:   mean(x-y)"
	
	if(plot) {
	h <- list(density = b$histobins * b$histoData/(b$subsetsExamined * (b$histoRight - b$histoLeft)),
		mids = seq(from=b$histoLeft, to=b$histoRight, length.out=b$histobins),
		statObs = statObs,
		statTail = statTail,
		testname= b$testname					
		)
	rv <- function(a) mean(a^2) - mean(a)^2
	denom <- sqrt(rv(x-y)/n)
	tvals <- h$mids/denom
	h$asymp <- dt(tvals, df=n-1)/denom
#	fud <- (tvals[histobins] - tvals[1])/(histoRight - histoLeft)
#	fud <- 1/sqrt(rv(x-y)/n)
#	h$asymp <- h$asymp /denom
	xh <<- h

	g <- histoplot(h)
	print(g)
	
		if(detail >= 1) {
			cat("\n\nPaired test using ", ifelse(method=="epairTest", "full enumeration", "random trials"))
			cat("\n\nP value (two tail) = \t", (b$nupper + b$nlower), "/", b$subsetsExamined, " = ", b$P2t, sep="")
			if(statObs < 0) cat("\nP value (one tail) = \t", b$nlower, "/", b$subsetsExamined, " = ", b$P1t, " (Test for X < Y)", sep="")
			if(statObs >= 0) cat("\nP value (one tail) = \t", b$nupper, "/", b$subsetsExamined, " = ", b$P1t, " (Test for X > Y)", sep="")
		}
		if(detail >= 2) {
#			cat("\n\nmean(X) = ", mean(x),"   \tn = ", k)
#			cat("\nmean(Y) = ", mean(y),"   \tn = ", ny,"\n")
		}
	}
	return(b)   

}

