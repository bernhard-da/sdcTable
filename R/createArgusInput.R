# Documentation: http://research.cbs.nl/casc/Software/TauManualV4.1.pdf

#' createArgusInput
#'
#' create required input-files and batch-file for tau-argus given an \code{sdcProblem-obj} object
#'
#' @param obj an object of class \code{\link{sdcProblem-class}} from \code{sdcTable}
#' @param typ (character) either \code{microdata} or \code{tabular}
#' @param verbose (logical) if TRUE, the contents of the batch-file are written to the prompt
#' @param path path, into which (temporary) files will be written to (amongst them being the batch-files).
#' Each file written to this folder belonging to the same problem contains a random id in its filename.
#' @param solver which solver should be used. allowed choices are
#' \itemize{
#' \item "FREE"
#' \item "CPLEX"
#' \item "XPRESS"
#' }
#' In case "CPLEX" is used, it is also mandatory to specify argument 'licensefile' which needs to be
#' the absolute path the the cplex license file
#' @param method secondary cell suppression algorithm, possible choices include:
#' \itemize{
#' \item MOD: modular approach. If specified, the following arguments in \code{...} can additionally be set:
#' \itemize{
#' \item MaxTimePerSubtable: number specifiying max. time (in minutes) spent for each subtable
#' \item SingleSingle: 0/1 (default=1)
#' \item SingleMultiple: 0/1 (default=1)
#' \item MinFreq: 0/1 (default=1)
#' }
#' \item GH: hypercube. If specified, the following arguments in \code{...} can additionally be set:
#' \itemize{
#' \item BoundPercentage: Default percentage to proctect primary suppressed cells, default 75
#' \item ModelSize: are we dealing with a small (0) or large (1) model? (default=1)
#' \item ApplySingleton: should singletons be additionally protected? 0/1 (default=1)
#' }
#' \item OPT: optimal cell suppression. If specified, the following arguments in \code{...} can additionally be set:
#' \itemize{
#' \item MaxComputingTime: number specifiying max. allowed computing time (in minutes)
#' }
#' }
#' @param primSuppRules rules for primary suppression, provided as a
#' \code{list}. For details, please have a look at the examples below.
#' @param responsevar which variable should be tabulated (defaults to frequencies). For details see tau-argus manual section 4.4.4.
#' @param shadowvar if specified, this variable is used to apply the safety rules, defaults to \code{responsevar}. For details see tau-argus manual section 4.4.4.
#' @param costvar if specified, this variable describes the costs of suppressing each individual cell. For details see tau-argus manual section 4.4.4.
#' @param requestvar if specified, this variable (0/1-coded) contains information about records that request protection. Records with 1 will be protected in case a corresponding request rule matches. It is ignored, if tabular input is used.
#' @param holdingvar if specified, this variable contains information about records that should be grouped together. It is ignored, if tabular input is used.
#' @param ... allows to specify additional parameters for selected suppression-method as described above
#' as well as \code{licensefile} in clase "CPLEX" was specified in argument \code{solver}.
#' @return the filepath to the batch-file
#' @export
#' @examples
#' # loading micro data from sdcTable
#' data("microData1", package="sdcTable")
#' microData1$num1 <- rnorm(mean=100, sd=25, nrow(microData1))
#' microData1$num2 <- round(rnorm(mean=500, sd=125, nrow(microData1)),2)
#' microData1$weight <- sample(10:100, nrow(microData1), replace=TRUE)
#'
#' dim.region <- data.frame(
#'   levels=c('@','@@@','@@@','@@@','@@@'),
#'   codes=c('Total', 'A','B','C','D'),
#'   stringsAsFactors=FALSE)
#'
#' dim.region_dupl <- data.frame(
#'   levels=c('@','@@@','@@@','@@@@@','@@@','@@@','@@@@@'),
#'   codes=c('Total', 'A','B','b1','C','D','d1'),
#'   stringsAsFactors=FALSE)
#'
#' dim.gender <- data.frame(
#'   levels=c('@','@@@','@@@'),
#'   codes=c('Total', 'male','female'),
#'   stringsAsFactors=FALSE)
#'
#' dimList <- list(region=dim.region, gender=dim.gender)
#' dimList_dupl  <- list(region=dim.region_dupl, gender=dim.gender)
#' dimVarInd <- c(1,2)
#' numVarInd <- 3:5
#' sampWeightInd <- 6
#'
#' freqVarInd <- weightInd <- NULL
#'
#' # creating an object of class \code{\link{sdcProblem-class}}
#' obj <- makeProblem(
#'   data=microData1,
#'   dimList=dimList,
#'   dimVarInd=dimVarInd,
#'   freqVarInd=freqVarInd,
#'   numVarInd=numVarInd,
#'   weightInd=weightInd,
#'   sampWeightInd=sampWeightInd)
#'
#' # creating an object of class \code{\link{sdcProblem-class}} containing "duplicated" codes
#' obj_dupl <- makeProblem(
#'   data=microData1,
#'   dimList=dimList_dupl,
#'   dimVarInd=dimVarInd,
#'   freqVarInd=freqVarInd,
#'   numVarInd=numVarInd,
#'   weightInd=weightInd,
#'   sampWeightInd=sampWeightInd)
#'
#' ## create primary suppression rules
#' primSuppRules <- list()
#' primSuppRules[[1]] <- list(type="freq", n=5, rg=20)
#' primSuppRules[[2]] <- list(type="p", n=5, p=20)
#' # other supported formats are:
#' # list(type="nk", n=5, k=20)
#' # list(type="zero", rg=5)
#' # list(type="mis", val=1)
#' # list(type="wgt", val=1)
#' # list(type="man", val=20)
#'
#' ## create batchInput object
#' bO_md1 <- createArgusInput(obj, typ="microdata", path=getwd(), solver="FREE", method="OPT",
#'   primSuppRules=primSuppRules, responsevar="num1")
#' bO_td1 <- createArgusInput(obj, typ="tabular", path=getwd(), solver="FREE", method="OPT")
#' bO_td2 <- createArgusInput(obj_dupl, typ="tabular", path=getwd(), solver="FREE", method="OPT")
#' \dontrun{
#' ## in case CPLEX should be used, it is required to specify argument licensefile
#' bO_md2 <- createArgusInput(obj, typ="microdata", path=getwd(), solver="CPLEX", method="OPT",
#'   primSuppRules=primSuppRules, responsevar="num1", licensefile="/path/to/my/cplexlicense")
#' }

createArgusInput <- function(
  obj, typ="microdata", verbose=FALSE, path=getwd(), solver="FREE", method,
  primSuppRules=NULL, responsevar=NULL, shadowvar=NULL, costvar=NULL,
  requestvar=NULL, holdingvar=NULL, ...) {

  if (class(obj) != "sdcProblem") {
    stop("argument 'obj' must be of class 'sdcProblem'.\n")
  }
  if (!typ %in% c("microdata", "tabular")) {
    stop("argument 'type' must be either 'microdata' or 'tabular'.\n")
  }

  args <- list(...)
  if (solver == "CPLEX" & is.null(args$licensefile)) {
    e <- c(
      dQuote("CPLEX"), "should be used, but",
      shQuote("licensefile"), "was not set."
    )
    stop(paste(e, collapse = " "), call. = FALSE)
  }

  if (typ == "microdata") {
    if (is.null(primSuppRules)) {
      e <- c(
        "primary suppression rules (argument 'primSuppRules') must be",
        "specified when using microdata as input!"
      )
      stop(paste(e, collapse = " "), call. = FALSE)
    }
    batchObj <- tauBatchInput_microdata(
      obj = obj,
      verbose = verbose,
      path = path,
      solver = solver,
      method = method,
      primSuppRules = primSuppRules,
      responsevar = responsevar,
      shadowvar = shadowvar,
      costvar = costvar,
      requestvar = requestvar,
      holdingvar = holdingvar, ...
    )
  }
  if (typ == "tabular") {
    if (!is.null(primSuppRules)) {
      message("ignoring argument 'primSuppRules'!\n")
      primSuppRules <- NULL
    }
    if (!is.null(requestvar)) {
      message("ignoring argument 'requestvar'!\n")
      requestvar <- NULL
    }
    if (!is.null(holdingvar)) {
      message("ignoring argument 'holdingvar'!\n")
      holdingvar <- NULL
    }

    batchObj <- tauBatchInput_table(
      obj = obj,
      verbose = verbose,
      path = path,
      solver = solver,
      method = method,
      responsevar = responsevar,
      shadowvar = shadowvar,
      costvar = costvar, ...
    )
  }
  ## write required files
  batchF <- writeBatchFile(batchObj)
  batchF <- normalizePath(batchF, winslash = "/", mustWork = TRUE)
  if (verbose) {
    ff <- paste("The batch-file", dQuote(batchF), "has the following content:")
    cat(paste0("\n", ff, "\n"))
    cat(paste(rep("-", nchar(ff)), collapse = ""), "\n")
    rr <- readLines(batchF, warn = FALSE)[-c(1:3)]
    cat(rr, sep = "\n")
  }
  return(invisible(batchF))
}
