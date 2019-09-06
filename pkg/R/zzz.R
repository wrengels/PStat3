# (c) William R. Engels, 2019

#' @importFrom utils packageDescription
#' @importFrom stats dchisq pchisq qchisq runif

.onAttach <- function(libname, pkgname){
	version <- packageDescription("PStat3", fields = "Version")
	packageStartupMessage("------------------------------------------------")
	packageStartupMessage("                  PStat3")
	packageStartupMessage(paste("               version", version))
	packageStartupMessage("translated from two older PStats")
	packageStartupMessage("------------------------------------------------")
}