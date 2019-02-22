#' query information for a specific cell in \code{\link{safeObj-class}} objects
#'
#' Function \code{\link{cellInfo}} is used to query information for a single table cell
#' for objects of class \code{\link{safeObj-class}}.
#'
#' @param object an object of class \code{\link{safeObj-class}}
#' @param characteristics a character vector specifying characteristics of the table cell that should be identified for each dimensional variable defining the table
#' @param varNames a character vector specifying variable names of dimensional variables defining the tables
#' @param verbose logical vector of length 1 defining verbosity, defaults to 'FALSE'
#'
#' @return a list containing the following calculated information
#' \itemize{
#' \item \code{cellID}: numeric vector of length 1 specifying the index of the cell within the final result dataset
#' \item \code{data}: a data.frame containing a single row with the index of the table cell of interest
#' \item \code{primSupp}: logical vector of length 1 that is 'TRUE' if the cell is a primary sensitive cell and 'FALSE' otherwise
#' \item \code{secondSupp}: logical vector of length 1 that is 'TRUE' if the cell is a secondary suppressed cell and 'FALSE' otherwise
#' }
#'
#' @examples
#' # load protected data (as created in the example
#' # of \code{\link{protectTable}})
#' sp <- searchpaths()
#' fn <- paste(sp[grep("sdcTable", sp)], "/data/protectedData.RData", sep="")
#' protectedData <- get(load(fn))
#' characteristics <- c('male', 'D')
#' varNames <- c('gender', 'region')
#' info <- cellInfo(protectedData, characteristics, varNames, verbose=FALSE)
#'
#' # show the info about this cell
#' str(info)
#'
#' @rdname cellInfo
#' @export cellInfo
#' @note Important: the \code{i}-th element of argument \code{characteristics} is uses as the desired characteristic for the dimensional variable specified at the \code{i}-th position of argument \code{varNames}!
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
cellInfo <- function(object, characteristics, varNames, verbose=FALSE) {
  paraList <- list()
  paraList[[1]] <- varNames
  paraList[[2]] <- characteristics
  paraList[[3]] <- verbose
  g_getCellInfo(object, input = paraList)
}
