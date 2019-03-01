#' change anonymization status of a specific cell
#'
#' Function \code{\link{changeCellStatus}} allows to change|modify the anonymization state
#' of single table cells for objects ofs class \code{\link{sdcProblem-class}}.
#'
#' @param object an object of class \code{\link{sdcProblem-class}}
#' @param characteristics a character vector specifying characteristics of the table cell that should be identified for each dimensional variable defining the table
#' @param varNames a character vector specifying variable names of dimensional variables defining the tables
#' @param rule character vector of length 1 specifying a valid anonymization code ('u', 'z', 'x', 's') to which the the cell under consideration should be set.
#' @param verbose logical vector of length 1 defining verbosity, defaults to 'FALSE'
#'
#' @return a \code{\link{sdcProblem-class}} object
#'
#' @examples
#' # load primary suppressed data (as created in the example
#' # of \code{\link{primarySuppression}})
#' sp <- searchpaths()
#' fn <- paste(sp[grep("sdcTable", sp)], "/data/problemWithSupps.RData", sep="")
#' problem <- get(load(fn))
#'
#' # we want to mark the cell region='D' and gender='male' primary sensitive
#' characteristics <- c('D', 'male')
#' varNames <- c('region', 'gender')
#' verbose <- TRUE
#' rule <- 'u'
#'
#' # looking at the distribution of anonymization states before...
#' print(table(getInfo(problem, 'sdcStatus')))
#'
#' # setting the specific cell as primary sensitive
#' problem <- changeCellStatus(problem, characteristics, varNames, rule, verbose)
#'
#' # having a second look at the anonymization states
#' print(table(getInfo(problem, 'sdcStatus')))
#'
#' @rdname changeCellStatus
#' @export changeCellStatus
#' @note Important: the \code{i}-th element of argument \code{characteristics} is uses as the desired characteristic for the dimensional variable specified at the \code{i}-th position of argument \code{varNames}!
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
changeCellStatus <- function(object, characteristics, varNames, rule, verbose=FALSE) {
  if (class(object) != "sdcProblem") {
    stop("changeCellStatus:: argument 'object' must be of class 'sdcProblem'!\n")
  }
  
  paraList <- list()
  paraList$names <- varNames
  paraList$codes <- characteristics
  paraList$verbose <- verbose
  
  cellID <- c_cellID(object, input = paraList)
  
  pI <- g_problemInstance(object)
  s_sdcStatus(pI) <- list(index = cellID, vals = rule)
  s_problemInstance(object) <- pI
  
  if (paraList$verbose) {
    freq <- g_freq(pI)[cellID]
    cat("--> The cell with ID=", cellID, "and Frequency", freq, "has been set to", rule, ".\n")
  }
  object
}
