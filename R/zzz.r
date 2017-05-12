.onAttach <- function(lib, pkg) {
  options(useFancyQuotes=FALSE)
  packageStartupMessage(paste("Package sdcTable",utils::packageVersion("sdcTable"),"has been loaded!\n"))
}
