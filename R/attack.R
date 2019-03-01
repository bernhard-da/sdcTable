#' attacking primary suppressed cells and calculating current lower and upper bounds
#'
#' Function \code{\link{attack}} is used to calculate lower and upper bounds for a given
#' sdcProblem (stored as object of class \code{\link{sdcProblem-class}}).
#' For all calculations the current suppression pattern is used when calculating solutions of the
#' attacker's problem.
#'
#' @param object an object of class \code{\link{sdcProblem-class}}
#' @param verbose a logical vector specifying if output should be verbose (TRUE) or not (FALSE)
#' @return a data.frame with column 'index' holding indices of primary suppressed cells and columns
#' 'bounds_min' and 'bounds_max' featuring calculated lower and upper bounds for each cell.
#' Column 'protected' shows if a given cell is accordingly protected (TRUE) or not (FALSE).
#'
#' @examples
#' # load problem (as it was created after performing primary suppression
#' # in the example of \code{\link{primarySuppression}})
#' sp <- searchpaths()
#' fn <- paste(sp[grep("sdcTable", sp)], "/data/problemWithSupps.RData", sep="")
#' problem <- get(load(fn))
#'
#' # calculate current lower|upper bounds given current suppression pattern
#' # (in this case consisting of primary suppressions only)
#' attack(problem, verbose=FALSE)
#' @rdname attack
#' @export attack
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
attack <- function(object, verbose=FALSE) {
  csp_cpp(sdcProblem = object, attackonly = TRUE, verbose = verbose)
}
