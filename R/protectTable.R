#' protecting \code{\link{sdcProblem-class}} objects
#'
#' Function \code{\link{protectTable}} is used to protect primary sensitive table cells
#' (that usually have been identified and set using
#' \code{\link{primarySuppression}}). The function protects primary
#' sensitive table cells according to the method that has been chosen and the
#' parameters that have been set. Additional parameters that are used to control
#' the protection algorithm are set using parameter \code{...}.
#'
#' The implemented heuristic method may have bugs that yield in not-fully protected tables. In case you encounter any problems,
#' please report them by filing an issue or use Tau-Argus (\url{http://neon.vb.cbs.nl/casc/tau.htm}) using \code{\link{createArgusInput}},
#' and \code{\link{runArgusBatchFile}}.
#'
#' @param object a \code{\link{sdcProblem-class}} object that has created using \code{\link{makeProblem}} and has been modified by \code{\link{primarySuppression}}
#' @param verbose (logical) defining if verbose output should be produced.
#' @param detectSingletons logical, should a singleton-detection procedure be run before protecting the data, defaults to \code{FALSE}.
#' @param ... parameters used in the protection algorithm that has been selected. Parameters that can be changed are:
#' @return an \code{\link{safeObj-class}} object
#' @examples
#' # load problem (as it was created after performing primary suppression
#' # in the example of \code{\link{primarySuppression}})
#' sp <- searchpaths()
#' fn <- paste(sp[grep("sdcTable", sp)], "/data/problemWithSupps.RData", sep="")
#' problem <- get(load(fn))
#'
#' # protect the table using the heuristic algorith and verbose output
#' protectedData <- protectTable(problem, verbose=TRUE, detectSingletons=TRUE)
#'
#' # showing a summary
#' summary(protectedData)
#'
#' # looking at the final table with result suppression pattern
#' print(getInfo(protectedData, type='finalData'))
#' @rdname protectTable
#' @export protectTable
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
protectTable <- function(object, detectSingletons=TRUE, verbose=FALSE, ...) {
  paraList <- genParaObj(selection='control.secondary', detectSingletons=detectSingletons, verbose=verbose, ...)
  supps_u <- length(g_primSupps(object@problemInstance))
  supps_x <- length(g_secondSupps(object@problemInstance))
  if (supps_u + supps_x == 0) {
    return(c_finalize(object=object, input=paraList))
  }
  print(paraList)
  out <- c_quick_suppression(object, input=paraList)
  out <- out$object
  invisible(c_finalize(object=out, input=paraList))
}