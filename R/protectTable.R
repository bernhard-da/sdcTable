#' protecting \code{\link{sdcProblem-class}} objects
#'
#' Function \code{\link{protectTable}} is used to protect primary sensitive table cells
#' (that usually have been identified and set using
#' \code{\link{primarySuppression}}). The function protects primary
#' sensitive table cells according to the method that has been chosen and the
#' parameters that have been set. Additional parameters that are used to control
#' the protection algorithm are set using parameter \code{...}.
#'
#' The implemented methods may have bugs that yield in not-fully protected tables. Especially
#' the usage of \code{OPT}, \code{HITAS} and \code{HYPERCUBE} in production is not
#' suggested as these methods may eventually be removed completely. In case you encounter any problems,
#' please report it or use Tau-Argus (\url{http://research.cbs.nl/casc/tau.htm}).
#'
#' @param object a \code{\link{sdcProblem-class}} object that has created using \code{\link{makeProblem}} and has been modified by \code{\link{primarySuppression}}
#' @param method a character vector of length 1 specifying the algorithm that should be used to protect the primary sensitive table cells. Allowed values are:
#' \itemize{
#' \item \code{OPT}: protect the complete problem at once using a cut and branch algorithm. The optimal algorithm should be used for small problem-instances only.
#' \item \code{HITAS}: split the overall problem in smaller problems. These problems are protected using a top-down approach.
#' \item \code{HYPERCUBE}: protect the complete problem by protecting sub-tables with a fast heuristic that is based on finding and suppressing geometric structures (n-dimensional cubes) that are required to protect primary sensitive table cells.
#' \item \code{SIMPLEHEURISTIC}: heuristic, quick procedure which might be applied to very large problem instances
#' }
#' @param ... parameters used in the protection algorithm that has been selected. Parameters that can be changed are:
#' \itemize{
#' \item general parameters include:
#' \itemize{
#' \item \code{verbose}: logical vector of length 1 defining if verbose output should be produced. Parameter \code{verbose} defaults to 'FALSE'
#' \item \code{save}: logical vector of length 1 defining if temporary results should be saved in the current working directory (TRUE) or not (FALSE). Parameter \code{save} defaults to 'FALSE' }
#' \item parameters used for HITAS|OPT procedures:
#' \itemize{
#' \item \code{solver}: character vector of length 1 defining the solver to be used. Currently available choices are limited to 'glpk'.
#' \item \code{timeLimit}: numeric vector of length 1 (or NULL) defining a time limit in minutes after which the cut and branch algorithm should stop and return a possible non-optimal solution. Parameter \code{safe} has a default value of 'NULL'
#' \item \code{maxVars}: a numeric vector of length 1 (or NULL) defining the maximum problem size in terms of decision variables for which an optimization should be tried. If the number of decision variables in the current problem are larger than parameter \code{maxVars}, only a possible non-optimal, heuristic solution is calculated. Parameter \code{safe} has a default value of 'NULL'
#' \item \code{fastSolution}: logical vector of length 1 defining if or if not the cut and branch algorithm will be started or if the possibly non-optimal heuristic solution is returned independent of parameter \code{maxVars}. Parameter \code{fastSolution} has a default value of 'FALSE'
#' \item \code{fixVariables}: logical vector of length 1 defining whether or not it should be tried to fix some variables to zero or one based on reduced costs early in the cut and branch algorithm. Parameter \code{fixVariables} has a default value of 'TRUE'
#' \item \code{approxPerc}: numeric vector of length 1 that defines a percentage for which a integer solution of the cut and branch algorithm is accepted as optimal with respect to the upper bound given by the (relaxed) solution of the master problem. Its default value is set to '10'
#' \item \code{useC}: boolean vector of length 1 defining if c++ implementation of the secondary cell suppression problem should be used, defaults to FALSE}
#' \item parameters used for HYPERCUBE procedure:
#' \itemize{
#' \item \code{protectionLevel}: numeric vector of length 1 specifying the required protection level for the HYPERCUBE-procedure. Its default value is 80
#' \item \code{suppMethod}: character vector of length 1 defining the rule on how to select the 'optimal' cube to protect a single sensitive cells. Possible choices are:
#' \itemize{
#' \item \code{minSupps}: minimize the number of additional secondary suppressions (this is also the default setting).
#' \item \code{minSum}: minimize the sum of counts of additional suppressed cells
#' \item \code{minSumLogs}: minimize the log of the sum of additional suppressed cells}
#' \item suppAdditionalQuader: logical vector of length 1 specfifying if additional cubes should be suppressed if any secondary suppressions in the 'optimal' cube are 'singletons'. Parameter \code{suppAdditionalQuader} has a default value of 'FALSE'}
#' \item parameter used for protectLinkedTables():
#' \itemize{
#' \item \code{maxIter}: numeric vector of length 1 specifying the maximal number of interations that should be make while trying to protect common cells of two different tables. The default value of parameter \code{maxIter} is 10}
#' \item parameters used for SIMPLEHEURISTIC procedure:
#' \itemize{
#' \item \code{detectSingletons}: logical, should a singleton-detection procedure be run before protecting the data, defaults to \code{FALSE}.}
#' }
#'
#' @return an \code{\link{safeObj-class}} object
#' @examples
#' # load problem (as it was created after performing primary suppression
#' # in the example of \code{\link{primarySuppression}})
#' sp <- searchpaths()
#' fn <- paste(sp[grep("sdcTable", sp)], "/data/problemWithSupps.RData", sep="")
#' problem <- get(load(fn))
#'
#' # protect the table using the 'HITAS' algorithm with verbose output
#' protectedData <- protectTable(problem, method='HITAS', verbose=TRUE, useC=TRUE)
#'
#' # showing a summary
#' summary(protectedData)
#'
#' # looking at the final table with result suppression pattern
#' print(getInfo(protectedData, type='finalData'))
#' @rdname protectTable
#' @export protectTable
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
protectTable <- function(object, method, ...) {
  if (!method %in% c("HITAS", "OPT", "HYPERCUBE", "SIMPLEHEURISTIC")) {
    stop("valid methods are 'SIMPLEHEURISTIC', 'HITAS', 'HYPERCUBE' or 'OPT'!\n")
  }
  
  paraList <- genParaObj(
    selection = "control.secondary",
    method = method, ...
  )
  supps_u <- length(g_primSupps(object@problemInstance))
  supps_x <- length(g_secondSupps(object@problemInstance))
  if (supps_u + supps_x == 0) {
    return(c_finalize(object = object, input = paraList))
  }
  
  if (method == "SIMPLEHEURISTIC") {
    out <- c_quick_suppression(object, input = paraList)
    out <- out$object
  } else {
    if (paraList$useC) {
      if (method == "OPT") {
        out <- c_opt_cpp(object = object, input = paraList)
      }
      if (method == "HITAS") {
        out <- c_hitas_cpp(object = object, input = paraList)
      }
    } else {
      out <- c_anon_worker(object, input = paraList)
    }
  }
  invisible(c_finalize(object = out, input = paraList))
}
