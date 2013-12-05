.onAttach <- function(lib, pkg) {
	packageStartupMessage(paste("Package sdcTable",packageVersion("sdcTable"),"has been loaded!\n"))
}
