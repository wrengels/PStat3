
# Exact test for comparing two groups and testing for differences in location
# (c) William R. Engels, 2017
#' 


#' @title
#' esumtest
#' @description
#' Compare sets of values, x and y, for mean differences with exact test. 
#' Uses efficient algorithms for either complete enumeration or monte carlo test. 
#' Resulting distribution is displayed graphically. 


#' @details
#' This function calls one of two C routines optimized for speed. One is "esumTest" to perform the full enumeration of all k-subsets of the pooled x and y values. 
#' The other, "rsumTest", selects random subsets.
#' The distribution of the test statistic is displayed as a histogram.
#' Subsets different in the same direction but at least as extreme as observed are shown in red.
#' Those in the opposite tail are colored blue
#' An orange curve shows the asymptotic t distribution. Differences between orange curve and true distribution provides evidence for usefulness of exact test.


#' @param x array of x values to compare with y
#' @param y array of y values
#' @param stat Set to \dQuote{diff} to use the difference between the mean of x and mean of y. Set to \dQuote{mean} to use just the x mean as test statistic. P value is not affected.
#' @param histobins number of bins for distribution
#' @param histoLeft left end of plotted distribution. Leave as NA to have it set automatically.
#' @param histoRight right end of distribution. NA means set automatically
#' @param safeSecs Stop the calculations after this many seconds
#' @param safeSets Do not attempt full enumeration if there are more than this many subsets
#' @param plot Display distribution histogram if TRUE.
#' @param detail How much detail to print after completion of test
#' @param B The number of random trials to use. Set to zero to use full enumeration.


#' @return Returns a list which includes the following
#' \item{x}{all values, including original x and y}
#' \item{n}{the total number of values}
#' \item{k}{the number of values in first set}
#' \item{B}{the number of requested random trials, or zero if full enumeration}
#' \item{statLeft, statRight}{the cutoffs for computing P values, one- or two-tail}
#' \item{histobins}{number of bins in distribution histogram}
#' \item{histoLeft, histoRight}{boundaries of histogram}
#' \item{safeSecs, safeSets}{time and subset limits to prevent getting stuck in C routine longer than intended}
#' \item{nlower, nupper}{the actual numbers of subsets found in the upper and lower tails}
#' \item{histoData}{observed number of subsets in each bin}
#' \item{subsetsExamined}{the actual, as opposed to requested, number of subsets tested}
#' \item{timespent}{the time in seconds spent in the C routine}
#' \item{P2t}{two-tail P value}
#' \item{P1t}{one-tail P value -- direction depends on the direction of observed data}
#' \item{testname}{which test statistic was used for histogram plot}


#' 
#' 
#' @examples
#' x <- rnorm(10, 5)
#' y <- rnorm(20, 7)
#' result <- e.sumtest(x,y) #performs full enumeration
#' result <- e.sumtest(x,y, B=1e6) #tests one million random subsets
#' 

#' @useDynLib PStat3
#' @export
e.sumtest <- function(x, y, stat=c("diff", "mean"), histobins=500, histoLeft=NA, histoRight=NA, safeSecs=30, safeSets=1e10, plot=T, detail=2, nreps=20, rseed=0, B=0) {
	k <- length(x)
	ny <- length(y)
	stat <- match.arg(stat)
	n <- k + ny
	totalSubsets <- choose(n, k)
	if(totalSubsets > safeSets & B==0){
		cat("\n")
		warning("\nThere are ", totalSubsets, " subsets to examine. The current limit is ", safeSets, ". Change parameter \"safeSets\" for different limit.\nSwitching to random subset method.")
		B <- 1e5
	} 
	if(totalSubsets > 1e9 & B==0) {
		warning(call.=F, immediate.=T, "\nStarting recursive examination of ", totalSubsets, " subsets.\n");
		print("\n")
		}
	if(stat=="diff") {
		z <- c(x,y) *(1/k + 1/ny) - sum(x,y)/(k*ny)
	}
	if(stat=="mean") {
		z <- c(x,y)/k
	}
	# Set up reasonable histogram boundaries if not specified
	statse <- sqrt(k * var(z))
	nse <- 4 # plus or minus nse standard errors
	statmean <- sum(z)*k/n
	if(is.na(histoLeft)) {
		histoLeft <- statmean - nse * statse
	}
	if(is.na(histoRight)) {
		histoRight <- statmean + nse * statse
	}

	statObs <- sum(z[1:k])
	statTail <- 2*sum(z)*k/n - statObs
	obsLeft <- statObs < statTail
	
	# Define left and right cutoff points
	rstatLeft <- ifelse(obsLeft, statObs, statTail)
	rstatRight <- ifelse(obsLeft, statTail, statObs)
	
	# Expand cutoff regions slightly to avoid floating point test for equality
	oneMinus <- 0.9999999
	rstatLeft <- ifelse(rstatLeft<0, rstatLeft * oneMinus, rstatLeft/oneMinus)
	rstatRight <- ifelse(rstatRight<0, rstatRight / oneMinus, rstatRight*oneMinus)
	
	# We don't expect bounds to be exactly 0, but just in case ...

    # Set up remaining variables
    nlower <- -10
    nupper <- nreps # just a way to pass nreps without another variable
    rhistoData <- 1:histobins
    rhistoData[5] <- rseed  # just a way to pass rseed without another variable
    subsetsExamined <- -42
    timespent <- 0.0
    rsafeSecs <- safeSecs
	method <- ifelse(B==0, "esumTest", "rsumTest")
    # Call C routine to do complete enumeration test!
    b <- .C(method,
		x=as.double(z),
		n=as.integer(n),
		k=as.integer(k),
		B=as.double(B),
		statLeft=as.double(rstatLeft),
		statRight=as.double(rstatRight),
		histobins=as.integer(histobins),
		histoLeft=as.double(histoLeft),
		histoRight=as.double(histoRight),
		safeSecs=as.integer(rsafeSecs),
		
		nlower=as.double(nlower),
		nupper=as.double(nupper),
		histoData=as.double(rhistoData),
		subsetsExamined=as.double(subsetsExamined),
		timespent=as.double(timespent)
		)
		
		b$P2t <- (b$nupper + b$nlower)/b$subsetsExamined
		b$P1t <- ifelse(obsLeft, b$nlower/b$subsetsExamined, b$nupper/b$subsetsExamined)
		b$statObs <- statObs
		b$statTail <- statTail
		b$testname <- "mean(X) - mean(Y)"
		b$testname <- ifelse(stat=="diff", "mean(X) - mean(Y)", "mean(X)")
		if(round(b$subsetsExamined) != totalSubsets & B==0) warning("Subsets examined not as predicted")
		if(b$subsetsExamined < 0) {cat("\nTime limit expired after about ", abs(b$subsetsExamined), " subsets. To change the time limit, set parameter \"safeSecs\" to the maximum number of seconds to allow. It is currently set at ", safeSecs,"\n")} else {
			if(plot) {
				h <- list(density = b$histobins * b$histoData/(b$subsetsExamined * (b$histoRight - b$histoLeft)),
					mids = seq(from=b$histoLeft, to=b$histoRight, length.out=b$histobins),
					statObs = statObs,
					statTail = statTail,
					testname= b$testname					
					)
# Find corresponding t values for mids
				rv <- function(z) mean(z^2) - mean(z)^2
				sp <- sqrt(((k-1)*rv(x) +(ny-1)*rv(y))/(n-2))
				denom <- sp * sqrt(1/k + 1/ny)
				if(stat == "diff") numer <- h$mids
				if(stat == "mean") numer <- h$mids - (sum(x,y) - k*h$mids)/ny
				tvals <- numer/denom
				if(stat == "diff"){
					h$asymp <- dt(tvals, df=n-2)/denom  # Use "/denom" as derivative of inverse function
				}
				if(stat == "mean") {
					h$asymp <- dt(tvals, df=n-2)
					dmdd <- (ny/(k + ny))*denom
					h$asymp <- h$asymp / dmdd
				}
				
				xh <<- h
				xh$histoData <<- b$histoData
				xh$tvals <<- tvals
						
				g <- histoplot(h)
				print(g)
			}
			if(detail >= 1) {
				cat("\n\nTest for mean difference ", ifelse(method=="esumTest", "(full enumeration of subsets)", "(random subsets)"))
				cat("\n\nP value (two tail) = \t", (b$nupper + b$nlower), "/", b$subsetsExamined, " = ", b$P2t, sep="")
				if(obsLeft) cat("\nP value (one tail) = \t", b$nlower, "/", b$subsetsExamined, " = ", b$P1t, " (Test for X < Y)", sep="")
				if(!obsLeft) cat("\nP value (one tail) = \t", b$nupper, "/", b$subsetsExamined, " = ", b$P1t, " (Test for X > Y)", sep="")
			}
			if(detail >= 2) {
				cat("\n\nmean(X) = ", mean(x),"   \tn = ", k)
				cat("\nmean(Y) = ", mean(y),"   \tn = ", ny)
				cat("\nTest included ", b$subsetsExamined, " subsets and took ", b$timespent, " seconds\n")
				
			}
		}
		return(b)   
}