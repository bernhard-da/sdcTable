#' perform primary suppression in \code{\link{sdcProblem-class}}-objects
#'
#' Function \code{\link{primarySuppression}} is used to identify and suppress primary
#' sensitive table cells in \code{\link{sdcProblem-class}} objects.
#' Argument \code{type} allows to select a rule that should be used to identify
#' primary sensitive cells. At the moment it is possible to identify and
#' suppress sensitive table cells using the frequency-rule, the nk-dominance
#' rule and the p-percent rule.
#'
#' @param object a \code{\link{sdcProblem-class}} object
#' @param type character vector of length 1 defining the primary suppression rule. Allowed types are:
#' \itemize{
#' \item \code{freq}: apply frequency rule with parameters \code{maxN} and \code{allowZeros}
#' \item \code{nk}: apply nk-dominance rule with parameters \code{n}, \code{k} and \code{numVarInd}
#' \item \code{p}: apply p-percent rule with parameters \code{p} and \code{numVarInd}
#' \item \code{pq}: apply pq-rule with parameters \code{p} and \code{q}
#' }
#' @param ... parameters used in the identification of primary sensitive cells. Parameters that can be modified|changed are:
#' \itemize{
#' \item \code{maxN}: numeric vector of length 1 used when applying the frequency rule. All cells having counts <= \code{maxN} are set as primary suppressed. The default value of \code{maxN} is 3.
#' \item \code{allowZeros}: logical vector of length 1 specifying if empty cells (count==0) should be considered sensitive when using the frequency rule. The default value of \code{allowZeros} is 'FALSE' so that empty cells are not considered primary sensitive by default.
#' \item \code{p}: numeric vector of length 1 specifying parameter \code{p} that is used when applying the p-percent rule with default value of 80.
#' \item \code{pq}: numeric vector of length 2 specifying parameters \code{p} and \code{q} that are used when applying the pq-rule with the default being c(25, 50).
#' \item \code{n}: numeric vector of length 1 specifying parameter \code{n} that is used when applying the nk-dominance rule. Parameter \code{n} is set to 2 by default.
#' \item \code{k}: numeric vector of length 1 specifying parameter \code{k} that is used when applying the nk-dominance rule. Parameter \code{n} is set to 85 by default.
#' \item \code{numVarInd}: numeric vector of length 1 specifying the index of the numerical variable that should be used to identify cells that are dominated by 2 (p-percent rule) or n (nk-dominance)-rule. If \code{type} is either 'nk', 'p' or 'pq', it is mandatory to specify \code{numVarInd}.
#' }
#' @return a \code{\link{sdcProblem-class}} object
#'
#' @examples
#' # load micro data
#' sp <- searchpaths()
#' fn <- paste(sp[grep("sdcTable", sp)], "/data/microData1.RData", sep="")
#' microData <- get(load(fn))
#'
#' # load problem (as it was created in the example in \code{\link{makeProblem}})
#' fn <- paste(sp[grep("sdcTable", sp)], "/data/problem.RData", sep="")
#' problem <- get(load(fn))
#'
#' # we have a look at the frequency table by gender and region
#' xtabs(rep(1, nrow(microData)) ~ gender + region, data=microData)
#'
#' # cell with region=='A' and gender=='female' has 2 units contributing to it
#' # this cell should be considered sensitive according the the freq-rule with 'maxN' equal to 2!
#' p1 <- primarySuppression(problem, type='freq', maxN=2)
#'
#' # we can also apply a p-percent rule with parameter 'p' being 30 as below.
#' # This is only possible if we are dealing with micro data and we also have to specify the index of
#' # a numeric variable.
#' p2 <- primarySuppression(problem, type='p', p=30, numVarInd=1)
#'
#' # looking at anonymization states we see, that one cell is primary suppressed (sdcStatus=='u')
#' # and the remaining cells are possible candidates for secondary suppression (sdcStatus=='s') given
#' # the frequency rule with parameter 'maxN=2'.
#' # Applying the p-percent rule with parameter 'p=30' resulted in two primary suppressions.
#' data.frame(p1.sdc=getInfo(p1, type='sdcStatus'), p2.sdc=getInfo(p2, type="sdcStatus"))
#'
#' @rdname primarySuppression
#' @export primarySuppression
#' @note the nk-dominance rule, the p-percent rule and the pq-rule can only be applied if micro data have been used as input data to function \code{\link{makeProblem}}.
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
primarySuppression <- function(object, type, ...) {
  start.time <- proc.time()
  if (!type %in% c("nk", "freq", "p", "pq")) {
    stop("valid types are 'nk', 'freq', 'p' or 'pq'!\n")
  }
  
  numVarsIndices <- g_numvar_ind(g_dataObj(object))
  paraList <- genParaObj(
    selection = "control.primary",
    numVarIndices = numVarsIndices, ...
  )
  
  if (type == "freq") {
    object <- c_rule_freq(object, input = paraList)
  }
  
  if (type == "nk") {
    if (is.na(paraList$numVarInd)) {
      stop("argument 'numVarInd' must be specified!\n")
    }
    object <- c_rule_nk(object, input = paraList)
  }
  
  if (type == "p") {
    if (is.na(paraList$numVarInd)) {
      stop("argument 'numVarInd' must be specified!\n")
    }
    object <- c_rule_p(object, input = paraList)
  }
  
  if (type == "pq") {
    if (is.na(paraList$numVarInd)) {
      stop("argument 'numVarInd' must be specified!\n")
    }
    object <- c_rule_pq(object, input = paraList)
  }
  
  elapsed.time <- g_elapsedTime(object) + (proc.time() - start.time)[3]
  s_elapsedTime(object) <- elapsed.time
  return(object)
}
