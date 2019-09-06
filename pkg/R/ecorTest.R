# Exact test for comparing magnitude of values in sets of pairs
# (c) William R. Engels, 2019

#' 
#' 


#' @title
#' ecortest
#' @description
#' Compare x versus y  
#' Uses efficient algorithms for either complete enumeration or monte carlo test. 
#' Resulting distribution is displayed graphically. 
#' Test statistic is Pearson correlation
#' Trials consist of permutations



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
#' @param detail How much detail to print after completion of test
#' @param nreps Trials in monte carlo version are performed in groups of size nreps. Optimum seems to be around 20
#' @param rseed seeds random numbers for monte carlo. Default is current time, but can be set programmatically if more than run per second is being done.
#' @param B The number of random trials to use. Set to zero to use full enumeration.


#' @return Returns a list which includes the following
#' \item{x}{vector containing one member of each pair}
#' \item{y}{vector same length as x but with other member of each pair} 
#' \item{n}{the number of data pairs}
#' \item{nTodo}{the number of monte carlo trials requested. Not used in full enumeration}
#' \item{statLeft}{correlations less than this are in left tail}
#' \item{statRight}{correlations greater than this are in right tail}
#' \item{histobins}{number of bins used in distribution histogram}
#' \item{histoLeft, histoRight}{boundaries of plotted histogram}
#' \item{safeSecs}{process is aborted if it exceeds this many seconds}
#' \item{safeSets}{process is aborted if it exceeds this many trials}
#' \item{nlower, nupper}{the actual numbers of subsets found in the left and right tail respectively}
#' \item{histoData}{observed number of subsets in each bin. length = histobins}
#' \item{subsetsExamined}{the actual, as opposed to requested, number of subsets tested}
#' \item{timespent}{how long was spent in C code doing the test}
#' \item{statObs}{the correlation coefficient corresponding to the observed permutation -- will be equal to either statLeft or statRight}
#' \item{statTail}{negative of statObs}
#' \item{P2t}{two-tail P value}
#' \item{P1t}{one-tail P value -- direction depends on the direction of observed data}
#' \item{testname}{"Pearson Correlation"}


#' 
#' 
#' @examples
#' x <- rnorm(10)
#' y <- x + rnorm(10, 5)
#' result <- e.cortest(x, y) # performs full enumeration and plots result
#' result <- e.cortest(x, y, B=1e6) # tests one million random permutations
#' 


#' @useDynLib PStat3
#' @export
e.cortest <- function(x, y, stat="diff", histobins=500, histoLeft=NA, histoRight=NA, safeSecs=30, safeSets=1e10, plot=T, detail=2, nreps=20, rseed=0, B=0) {
	n <- length(y)
	if(length(x) != n) stop("Values must be pairs")
	logTotalSubsets <- lgamma(n+1)
	totalSubsets <- ifelse(logTotalSubsets< 100, factorial(n), Inf)
	if(logTotalSubsets > log(safeSets) & B==0){
		cat("\n")
			B <- 1e5
		warning("\nThere are ", totalSubsets, " permutations to examine. The current limit is ", safeSets, ". \nChange parameter \"safeSets\" for different limit.\nSwitching to random subset method with default B = ", B, "\n")
	} 
	statObs <- cor(x,y)
	statTail <- -statObs
	maxCorrelation <- cor(sort(x), sort(y))
	minCorrelation <- cor(sort(x), sort(y, decreasing=T))
	rv <- function(a) mean(a^2) - mean(a)^2
	acon <- n * sqrt(rv(x) * rv(y))
	if(acon > 0) acon <- 1/acon
	ccon <- -mean(x) * mean(y) * acon

	if(is.na(histoLeft)){
		histoLeft <- minCorrelation
	}
	if(is.na(histoRight)){
		histoRight <- maxCorrelation
	}
	nlower <- 0
	nupper <- nreps
	if(histobins < 5) histobins <- 5   #be sure there is space to send C some calculations
	rhistoData <- 1:histobins
	subsetsExamined <- -42
	timespent <- 0.0
	almost1 <- 0.999999
	
	rhistoData[1] <- maxCorrelation
	rhistoData[2] <- minCorrelation
	rhistoData[3] <- acon
	rhistoData[4] <- ccon
	rhistoData[5] <- rseed
	
	
	method <- ifelse(B==0, "ecorTest", "rcorTest")
	b <- .C(method,
		x=as.double(x),
		y=as.double(y),
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
	b$testname <- "Pearson Correlation"
	
	if(method=="ecorTest" && b$subsetsExamined < totalSubsets) {
		cat("\n\nFull enumeration test timed out before finishing.\nThere are ", totalSubsets, " permutations, but only time for ", b$subsetsExamined)
		cat("\nFull enumeration test is only valid when all permutations are generated")
		cat("\nTo perform this test, try random permutations by setting parameter \"B\", or increasing time with parameter \"safeSecs\"\n")
		return(b)
	}
	
	if(plot) {
	h <- list(density = b$histobins * b$histoData/(b$subsetsExamined * (b$histoRight - b$histoLeft)),
		mids = seq(from=b$histoLeft, to=b$histoRight, length.out=b$histobins),
		statObs = statObs,
		statTail = statTail,
		testname= b$testname					
		)
	
	tvals <- h$mids * sqrt((n-2)/(1-h$mids^2))
	drdt <- (n-2)/(n-2 + tvals^2)^(3/2)
	h$asymp <- dt(tvals, df=n-2)
	h$asymp <- h$asymp/drdt # Use derivative dr/dt to adjust for transformed widths

	h$tvals <- tvals

	fud <- (tvals[histobins] - tvals[1])/(histoRight - histoLeft)
	wdth <- (histoRight - histoLeft)/histobins
	xh <<- h

	g <- histoplot(h)
	print(g)
		} # if plot
	
		if(detail >= 1) {
			cat("\n\nTest for correlation ", ifelse(method=="ecorTest", "(full enumeration of permutations)", "(random permutations)"))
			cat("\n\nP value (two tail) = \t", (b$nupper + b$nlower), "/", b$subsetsExamined, " = ", b$P2t, sep="")
			if(statObs < 0) cat("\nP value (one tail) = \t", b$nlower, "/", b$subsetsExamined, " = ", b$P1t, " (Test for negative correlation)\n", sep="")
			if(statObs >= 0) cat("\nP value (one tail) = \t", b$nupper, "/", b$subsetsExamined, " = ", b$P1t, " (Test for positive correlation)\n", sep="")
		}
		if(detail >= 2) {
			cat("\nThere were ", b$n, " pairs of data\nObserved correlation = ", b$statObs)
			cat("\n", b$subsetsExamined," permutations were examined in ", b$timespent, " seconds.\n")

	}
	return(b)   

}

