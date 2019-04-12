#' Apply primary suppression
#'
#' Function [primarySuppression()] is used to identify and suppress primary
#' sensitive table cells in [sdcProblem-class] objects.
#' Argument `type` allows to select a rule that should be used to identify
#' primary sensitive cells. At the moment it is possible to identify and
#' suppress sensitive table cells using the frequency-rule, the nk-dominance
#' rule and the p-percent rule.
#'
#' @param object a [sdcProblem-class] object
#' @param type character vector of length 1 defining the primary suppression
#' rule. Allowed types are:
#' - `freq`: apply frequency rule with parameters `maxN` and `allowZeros`
#' - `nk`: apply nk-dominance rule with parameters `n`, `k` and `numVarInd`
#' - `p`: apply p-percent rule with parameters `p` and `numVarInd`
#' - `pq`: apply pq-rule with parameters `p` and `q`
#' @param ... parameters used in the identification of primary sensitive cells.
#' Parameters that can be modified|changed are:
#' - `maxN`: numeric vector of length 1 used when applying the frequency rule.
#' All cells having counts <= `maxN` are set as primary suppressed. The default
#' value of `maxN` is `3`.
#' - `allowZeros`: logical vector of length 1 specifying if empty cells
#' (count==0) should be considered sensitive when using the frequency rule.
#' The default value of `allowZeros` is `FALSE` so that empty cells are not
#' considered primary sensitive by default.
#' - `p`: numeric vector of length 1 specifying parameter `p` that is used
#' when applying the p-percent rule with default value of `80`.
#' - `pq`: numeric vector of length 2 specifying parameters `p` and `q` that
#' are used when applying the pq-rule with the default being c(`25`, `50`).
#' - `n`: numeric vector of length 1 specifying parameter `n` that is used
#' when applying the nk-dominance rule. Parameter `n` is set to `2` by default.
#' - `k`: scalar numeric specifying parameter `k` that is used
#' when applying the nk-dominance rule. Parameter `n` is set to `85` by default.
#' - `numVarName`: character scalar specifying the name
#' of the numerical variable that should be used to identify cells that are
#' dominated by dominance rules (`p-rule`, `pq-rule` or `nk-rule`).
#' If `type` is either 'nk', 'p' or 'pq', it is mandatory to
#' specify either `numVarInd` or `numVarName`.
#' - `numVarInd`: same as `numVarName` but a scalar numeric
#' specifying the index of the variable is expected. If both `numVarName`
#' and `numVarInd` are specified, `numVarName` is used.
#' @return a [sdcProblem-class] object
#' @md
#' @export
#' @note the nk-dominance rule, the p-percent rule and the pq-rule can only
#' be applied if micro data have been used as input data to function [makeProblem()]
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
#' @md
#' @examples
#' # load micro data
#' data("microData1", package = "sdcTable")
#'
#' # load problem (as it was created in the example in ?makeProblem
#' data("problem", package = "sdcTable")
#'
#' # we have a look at the frequency table by gender and region
#' xtabs(rep(1, nrow(microData1)) ~ gender + region, data = microData1)
#'
#' # 2 units contribute to cell with region=='A' and gender=='female'
#' # --> this cell is considered sensitive according the the
#' # freq-rule with 'maxN' equal to 2!
#' p1 <- primarySuppression(
#'   object = problem,
#'   type = "freq",
#'   maxN = 2
#' )
#'
#' # we can also apply a p-percent rule with parameter "p" being 30 as below.
#' # This is only possible if we are dealing with micro data and we also
#' # have to specify the name of a numeric variable.
#' p2 <- primarySuppression(
#'   object = problem,
#'   type = "p",
#'   p = 30,
#'   numVarName = "val"
#' )
#'
#' # looking at anonymization states we see, that one cell is primary
#' # suppressed (sdcStatus == "u")
#' # the remaining cells are possible candidates for secondary cell
#' # suppression (sdcStatus == "s") given the frequency rule with
#' # parameter "maxN = 2".
#' #
#' # Applying the p-percent rule with parameter 'p = 30' resulted in
#' # two primary suppressions.
#' data.frame(
#'   p1_sdc = getInfo(p1, type = "sdcStatus"),
#'   p2_sdc = getInfo(p2, type = "sdcStatus")
#' )
primarySuppression <- function(object, type, ...) {
  start.time <- proc.time()
  if (!type %in% c("nk", "freq", "p", "pq")) {
    stop("valid types are 'nk', 'freq', 'p' or 'pq'!\n")
  }

  data_obj <- g_dataObj(object)
  dt <- g_raw_data(data_obj)

  numVarsIndices <- g_numvar_ind(g_dataObj(object))

  paraList <- genParaObj(
    selection = "control.primary",
    numVarIndices = numVarsIndices, ...
  )

  if (type == "freq") {
    object <- c_rule_freq(object, input = paraList)
  }

  if (type %in% c("nk", "p", "pq")) {
    pp <- list(...)
    if (is.na(paraList$numVarInd) & is.null(pp$numVarName)) {
      stop("Please specify either `numVarInd` or `numVarName`.", call. = FALSE)
    }

    cn <- colnames(dt)
    if (!is.null(pp$numVarName)) {
      if (!pp$numVarName %in% cn) {
        e <- "Variable specified in `numVarName` does not exist in the data."
        stop(e, call. = FALSE)
      }
      paraList$numVarInd <- match(pp$numVarName, cn)
    }
  }

  if (type == "nk") {
    object <- c_rule_nk(object, input = paraList)
  }

  if (type == "p") {
    object <- c_rule_p(object, input = paraList)
  }

  if (type == "pq") {
    object <- c_rule_pq(object, input = paraList)
  }

  elapsed.time <- g_elapsedTime(object) + (proc.time() - start.time)[3]
  s_elapsedTime(object) <- elapsed.time
  return(object)
}
