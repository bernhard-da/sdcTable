#' sdcProb2df
#'
#' returns a \code{data.frame} or \code{data.table} from an \code{\link{sdcProblem-class}}-object which contains the
#' current state of the problem instance.
#'
#' @param obj an object of class \code{\link{sdcProblem-class}} from \code{sdcTable}
#' @param addDups (logical), if \code{TRUE}, duplicated cells are included in the output
#' @param addNumVars (logical), if \code{TRUE}, numerical variables (if defined in \code{\link{makeProblem}}) will
#' be included in the output.
#' @param dimCodes (character) allows to specify in which coding the dimensional variables should be returned. Possible choices are:
#' \itemize{
#' \item "both": both original and internally used, standardized codes are included in the output
#' \item "original": only original codes of dimensional variables are included in the output
#' \item "default": only internally used, standardized codes are included in the output
#' }
#'
#' @return a \code{data.table} containing information about all cells of the given sdc problem
#' instance is returned.
#' @export
#' @examples
#' ## have a look at ?makeProblem
sdcProb2df <- function(obj, addDups=TRUE, addNumVars=FALSE, dimCodes="both") {
  if (length(dimCodes) != 1) {
    stop("Length of argument 'codes' must equal 1!\n")
  }
  if (!dimCodes %in% c("both", "original", "default")) {
    stop("Argument 'dimCodes'  must be either 'both', 'original' or 'default'!\n")
  }
  dt <- g_df(obj, addDups = addDups, addNumVars = addNumVars)
  
  dimV_d <- obj@dimInfo@vNames
  dimV_o <- paste0(obj@dimInfo@vNames, "_o")
  
  if (dimCodes == "original") {
    dt <- dt[, -c(match(dimV_d, names(dt))), with = F]
    cn <- names(dt)
    cn[match(dimV_o, cn)] <- dimV_d
    setnames(dt, cn)
  }
  if (dimCodes == "default") {
    dt <- dt[, -c(match(dimV_o, names(dt))), with = F]
  }
  dt
}
