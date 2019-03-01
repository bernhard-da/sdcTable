#' set information of \code{\link{sdcProblem-class}}- or \code{\link{problemInstance-class}} objects
#'
#' Function \code{\link{getInfo}} is used to query information from
#' \code{\link{sdcProblem-class}}- or \code{\link{problemInstance-class}} objects
#'
#' @param object an object of class \code{\link{sdcProblem-class}} or \code{\link{problemInstance-class}}
#' @param type a character vector of length 1 specifying the the information that should be changed or modified, valid choices are:
#' \itemize{
#' \item \code{lb}: slot 'lb' of input \code{object} if it is of class \code{\link{problemInstance-class}} or slot 'lb' within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}}
#' \item \code{ub}: slot 'ub' of input \code{object} if it is of class \code{\link{problemInstance-class}} or slot 'ub' within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}}
#' \item \code{LPL}: slot 'LPL' of input \code{object} if it is of class \code{\link{problemInstance-class}} or slot 'LPL' within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}}
#' \item \code{SPL}: slot 'SPL' of input \code{object} if it is of class \code{\link{problemInstance-class}} or slot 'SPL' within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}}
#' \item \code{UPL}: slot 'UPL' of input \code{object} if it is of class \code{\link{problemInstance-class}} or slot 'UPL' within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}}
#' \item \code{sdcStatus}:  slot 'sdcStatus' of input \code{object} if it is of class \code{\link{problemInstance-class}} or slot 'sdcStatus' within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}} }
#' @param index numeric vector defining cell-indices for which which values in a specified slot should be changed|modified
#' @param input numeric or character vector depending on argument \code{type} with its length matching the length of argument \code{index}
#' \itemize{
#' \item character vector if type matches 'sdcStatus'
#' \item a numeric vector if type matches 'lb', 'ub', 'LPL', 'SPL' or 'UPL'
#' }
#'
#' @return a \code{\link{sdcProblem-class}}- or \code{\link{problemInstance-class}} object
#'
#' @examples
#' # load primary suppressed data (created in the example of \code{\link{primarySuppression}})
#' sp <- searchpaths()
#' fn <- paste(sp[grep("sdcTable", sp)], "/data/problemWithSupps.RData", sep="")
#' problem <- get(load(fn))
#'
#' # which is the overall total?
#' index.tot <- which.max(getInfo(problem, 'freq'))
#' index.tot
#'
#' # we see that the cell with index.tot==1 is the overall total and its
#' # anonymization state of the total can be extracted as follows:
#' print(getInfo(problem, type='sdcStatus')[index.tot])
#'
#' # we want this cell to never be suppressed
#' problem <- setInfo(problem, type='sdcStatus', index=index.tot, input='z')
#'
#' # we can verify this:
#' print(getInfo(problem, type='sdcStatus')[index.tot])
#'
#' # changing slot 'UPL' for all cells
#' inp <- data.frame(strID=getInfo(problem,'strID'), UPL_old=getInfo(problem,'UPL'))
#' inp$UPL_new <- inp$UPL_old+1
#' problem <- setInfo(problem, type='UPL', index=1:nrow(inp), input=inp$UPL_new)
#'
#' @rdname setInfo
#' @export setInfo
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setInfo <- function(object, type, index, input) {
  if (!class(object) %in% c("sdcProblem", "problemInstance")) {
    stop("setInfo:: argument 'object' must be of class 'sdcProblem' or 'problemInstance'!\n")
  }
  
  if (!type %in% c("lb", "ub", "LPL", "SPL", "UPL", "sdcStatus")) {
    stop("setInfo:: check argument 'type'!\n")
  }
  
  if (class(object) == "sdcProblem") {
    pI <- g_problemInstance(object)
  } else {
    pI <- object
  }
  
  pI <- set.problemInstance(
    pI,
    type = type,
    input = list(index = index, values = input)
  )
  
  if (class(object) == "sdcProblem") {
    s_problemInstance(object) <- pI
  } else {
    object <- pI
  }
  object
}
