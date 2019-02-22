#' protect two \code{\link{sdcProblem-class}} objects that have common cells
#'
#' \code{\link{protectLinkedTables}} can be used to protect tables, that have
#' common cells. It is of course required that after the anonymization process
#' has finished, all common cells have the same anonymization state in both
#' tables.
#'
#' @param objectA a \code{\link{sdcProblem-class}} object
#' @param objectB a \code{\link{sdcProblem-class}} object
#' @param commonCells a list object defining common cells in \code{objectA} and \code{objectB}. For each variable that has one or more common codes in both tables, a list element needs to be specified.
#' \itemize{
#' \item List-elements of length 3: Variable has exact same levels and structure in both tables
#' \itemize{
#' \item \code{first element}: character vector of length 1 specifying the variable name in argument \code{objectA}
#' \item \code{second element}: character vector of length 1 specifying the variable name in argument \code{objectB}
#' \item \code{third element}: character vector of length 1 being with keyword \code{ALL} }
#' \item List-elements of length 4: Variable has different codes and levels in tables \code{objectA} and \code{objectB}
#' \itemize{
#' \item \code{first element}: character vector of length 1 specifying the variable name in argument \code{objectA}
#' \item \code{second element}: character vector of length 1 specifying the variable name in argument \code{objectB}
#' \item \code{third element}: character vector defining codes within \code{objectA}
#' \item \code{fourth element}: character vector with length that equals the length of the third list-element. The vector defines codes of the variable in \code{objectB} that match the codes given in the third list-element for \code{objectA}.
#' }
#' }
#' @param method a character vector of length 1 specifying the algorithm that should be used to protect the primary sensitive table cells. Allowed values are:
#' \itemize{
#' \item \code{HITAS}:
#' \item \code{SIMPLEHEURISTIC}:
#' \item \code{OPT}: }
#' @param ... additional arguments to control the secondary cell suppression algorithm. For details, see \code{\link{protectTable}}.
#'
#' @return a list of length 2 with each list-element being an \code{\link{safeObj-class}} object
#'
#' @examples
#' \dontrun{
#' # load micro data for further processing
#' sp <- searchpaths()
#' fn <- paste(sp[grep("sdcTable", sp)], "/data/microData2.RData", sep="")
#' microData <- get(load(fn))
#'
#' # table1: defined by variables 'gender' and 'ecoOld'
#' microData1 <- microData[,c(2,3,5)]
#'
#' # table2: defined by variables 'region', 'gender' and 'ecoNew'
#' microData2 <- microData[,c(1,2,4,5)]
#'
#' # we need to create information on the hierarchies
#' # variable 'region': exists only in microDat2
#' d_region <- hier_create(root = "Tot", nodes = c("R1", "R2"))
#'
#' # variable 'gender': exists in both datasets
#' d_gender <- hier_create(root = "Tot", nodes = c("m", "f"))
#'
#' # variable 'eco1': exists only in microDat1
#' d_eco1 <- hier_create(root = "Tot", nodes = c("A", "B"))
#' d_eco1 <- hier_add(d_eco1, root = "A", nodes = c("Aa", "Ab"))
#' d_eco1 <- hier_add(d_eco1, root = "B", nodes = c("Ba", "Bb"))
#' 
#' # variable 'ecoNew': exists only in microDat2
#' d_eco2 <- hier_create(root = "Tot", nodes = c("C", "D"))
#' d_eco2 <- hier_add(d_eco2, root = "C", nodes = c("Ca", "Cb", "Cc"))
#' d_eco2 <- hier_add(d_eco2, root = "D", nodes = c("Da", "Db", "Dc"))
#' 
#' # creating objects holding information on dimensions
#' dl1 <- list(gender = d_gender, ecoOld = d_eco1)
#' dl2 <- list(region = d_region, gender = d_gender, ecoNew = d_eco2)
#'
#' # creating input objects for further processing. For details have a look at
#' # \code{\link{makeProblem}}.
#' p1 <- makeProblem(
#'   data = microData1, 
#'   dimList = dl1, 
#'   dimVarInd = 1:2,
#'   numVarInd = 3
#' )
#' p2 <- makeProblem(
#'   data = microData2, 
#'   dimList = dl2, 
#'   dimVarInd = 1:3,
#'   numVarInd = 4
#' )
#'
#' # the cell specified by gender == "Tot" and ecoOld == "A"
#' # is one of the common cells! -> we mark it as primary suppression
#' p1 <- changeCellStatus(
#'   object = p1, 
#'   characteristics = c("Tot", "A"),
#'   varNames = c("gender", "ecoOld"), 
#'   rule = "u", 
#'   verbose = FALSE
#' )
#'
#' # the cell specified by region == "Tot" and gender == "f" and ecoNew == "C"
#' # is one of the common cells! -> we mark it as primary suppression
#' p2 <- changeCellStatus(
#'   object = p2, 
#'   characteristics = c("Tot", "f", "C"),
#'   varNames = c("region", "gender", "ecoNew"), 
#'   rule = "u", 
#'   verbose = FALSE
#' )
#'
#' # specifying input to define common cells
#' common_cells <- list()
#'
#' # variable "gender"
#' common_cells$v.gender <- list()
#' common_cells$v.gender[[1]] <- "gender" # variable name in "p1"
#' common_cells$v.gender[[2]] <- "gender" # variable name in "p2"
#' # "gender" has equal characteristics on both datasets -> keyword "ALL"
#' common_cells$v.gender[[3]] <- "ALL"
#'
#' # variables: "ecoOld" and "ecoNew"
#' common_cells$v.eco <- list()
#' common_cells$v.eco[[1]] <- "ecoOld"   # variable name in "p1"
#' common_cells$v.eco[[2]] <- "ecoNew"   # variable name in "p2"
#'
#' # vector of common characteristics: "A" and "B" in variable "ecoOld" in 'p1'
#' common_cells$v.eco[[3]] <- c("A", "B")
#' # correspond to characteristics "C" and "D" in variable "ecoNew" in "p2"
#' common_cells$v.eco[[4]] <- c("C", "D")
#'
#' # protect the linked data
#' result <- protectLinkedTables(
#'   objectA = p1, 
#'   objectB = p2, 
#'   commonCells = common_cells, 
#'   method = "HITAS", 
#'   verbose = TRUE
#' )
#'
#' # having a look at the results
#' result.tab1 <- result[[1]]
#' result.tab2 <- result[[2]]
#' summary(result.tab1)
#' summary(result.tab2)
#' }
#'
#' @rdname protectLinkedTables
#' @export protectLinkedTables
#' @seealso \code{\link{protectTable}}
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
protectLinkedTables <- function(objectA, objectB, commonCells, method, ...) {
  f.calcCommonCellIndices <- function(input1, input2, commonCells) {
    id1 <- id2 <- NULL
    dat1 <- g_df(input1); dat1$id1 <- 1:nrow(dat1)
    dat2 <- g_df(input2); dat2$id2 <- 1:nrow(dat2)
    dat1$strID <- dat1$sdcStatus <- NULL
    dat2$strID <- dat2$sdcStatus <- NULL
    dI1 <- g_dimInfo(input1)
    dI2 <- g_dimInfo(input2)

    # restrict to totals in non-overlapping variables in dataset1
    totvars <- setdiff(dI1@vNames, sapply(commonCells, function(x) x[[1]]))
    if (length(totvars) > 0) {
      for (i in 1:length(totvars)) {
        cmd <- paste0("dat1 <- dat1[", totvars[i], "=='", dI1@dimInfo[[totvars[i]]]@codesDefault[1], "']")
        eval(parse(text = cmd))
      }
    }
    # restrict to totals in non-overlapping variables in dataset2
    totvars <- setdiff(dI2@vNames, sapply(commonCells, function(x) x[[2]]))
    if (length(totvars) > 0) {
      for (i in 1:length(totvars)) {
        cmd <- paste0("dat2 <- dat2[", totvars[i], "=='", dI2@dimInfo[[totvars[i]]]@codesDefault[1], "']")
        eval(parse(text = cmd))
      }
    }

    for (i in 1:length(commonCells)) {
      if (length(commonCells[[i]]) == 4) {
        cmd1 <- paste0("dat1 <- dat1[,tmpxxx_V", i, ":=", commonCells[[i]][[1]], "_o]")
        cmd2 <- paste0("dat2 <- dat2[,tmpxxx_V", i, ":=", commonCells[[i]][[2]], "_o]")
        eval(parse(text = cmd1))
        eval(parse(text = cmd2))

        # recode different codes to those of dataset1
        codes <- commonCells[[i]][[4]]
        codesX <- commonCells[[i]][[3]]
        for  (z in 1:length(codes)) {
          if (codes[z] != codesX[z]) {
            cmd <- paste0("dat2 <- dat2[tmpxxx_V", i, "=='", codes[z], "',tmpxxx_V", i, ":='", codesX[z], "']")
            eval(parse(text = cmd))
          }
        }
      } else {
        # nothing to do, codes are the same!
        cmd1 <- paste0("dat1[,tmpxxx_V", i, ":=", commonCells[[i]][[1]], "_o]")
        cmd2 <- paste0("dat2[,tmpxxx_V", i, ":=", commonCells[[i]][[2]], "_o]")
        eval(parse(text = cmd1))
        eval(parse(text = cmd2))
      }
    }

    kV1 <- names(dat1)[grep("tmpxxx", names(dat1))]; setkeyv(dat1, kV1)
    kV2 <- names(dat2)[grep("tmpxxx", names(dat2))]; setkeyv(dat2, kV2)
    mm <- merge(dat1, dat2); setkey(mm, id1)
    if (any(mm$freq.x != mm$freq.y)) {
      stop("Error: common cells must have same values!\n")
    }
    return(list(commonInd1 = mm$id1, commonInd2 = mm$id2))
  }

  f.checkCommonCells <- function(suppPattern1, suppPattern2, commonCellIndices) {
    indOK <- TRUE
    if (any(suppPattern1[commonCellIndices[[1]]] != suppPattern2[commonCellIndices[[2]]]))
      indOK <- FALSE
    return(indOK)
  }

  ### arguments ###
  if (!method %in% c("HITAS", "OPT", "SIMPLEHEURISTIC")) {
    stop("valid methods are 'HITAS', 'OPT' or 'SIMPLEHEURISTIC'!\n")
  }

  paraList <- genParaObj(selection = "control.secondary", method = method, ...)

  ### first run
  if (method == "SIMPLEHEURISTIC") {
    outA <- c_quick_suppression(objectA, input = paraList)$object
    outB <- c_quick_suppression(objectB, input = paraList)$object
  } else {
    if (paraList$useC) {
      if (method == "OPT") {
        outA <- c_opt_cpp(object = objectA, input = paraList)
        outB <- c_opt_cpp(object = objectB, input = paraList)
      }
      if (method == "HITAS") {
        outA <- c_hitas_cpp(object = objectA, input = paraList)
        outB <- c_hitas_cpp(object = objectB, input = paraList)
      }
    } else {
      outA <- c_anon_worker(object = objectA, input = paraList)
      outB <- c_anon_worker(object = objectB, input = paraList)
    }
  }

  pI.A <- g_problemInstance(outA)
  pI.B <- g_problemInstance(outB)
  # calc original primary suppressions

  origPrimSupp1Index <- g_primSupps(pI.A)
  origPrimSupp2Index <- g_primSupps(pI.B)

  # no primary suppressions
  if (length(origPrimSupp1Index) + length(origPrimSupp2Index) == 0) {
    if (paraList$verbose) {
      cat("\n===> no primary suppressions. All common cells have the same anonymity-status! [Finished]\n")
    }
    outA <- c_finalize(object = outA, input = paraList)
    outB <- c_finalize(object = outB, input = paraList)
    return(list(outObj1 = outA, outObj2 = outB))
  }

  # calculate commonCells:
  commonCellIndices <- f.calcCommonCellIndices(outA, outB, commonCells)

  # suppression patterns after the first run
  suppPatternA <- g_suppPattern(pI.A)
  suppPatternB <- g_suppPattern(pI.B)

  indOK <- f.checkCommonCells(suppPatternA, suppPatternB, commonCellIndices)
  counter <- 1
  if (!indOK) {
    if (paraList$verbose) {
      cat("we have to start the iterative procedure!\n")
    }
    runInd <- TRUE
    while (runInd) {
      x <- cbind(suppPatternA[commonCellIndices[[1]]], suppPatternB[commonCellIndices[[2]]])
      index <- list()
      i1 <- which(x[, 1] == 0 & x[, 2] == 1)
      i2 <- which(x[, 1] == 1 & x[, 2] == 0)
      index[[1]] <- commonCellIndices[[1]][i1]
      index[[2]] <- commonCellIndices[[2]][i2]

      for (j in 1:2) {
        if (length(index[[j]]) > 0) {
          if (j == 1) {
            pI.A <- g_problemInstance(outA)
            s_sdcStatus(pI.A) <- list(
              index = index[[j]],
              vals = rep("u", length(index[[j]]))
            )
            s_problemInstance(outA) <- pI.A
            s_indicesDealtWith(outA) <- NULL
            s_startJ(outA) <- 1
            s_startI(outA) <- 1
            outA <- c_anon_worker(outA, input = paraList)
          } else {
            pI.B <- g_problemInstance(outB)
            s_sdcStatus(pI.B) <- list(
              index = index[[j]],
              vals = rep("u", length(index[[j]]))
            )
            s_problemInstance(outB) <- pI.B
            s_indicesDealtWith(outB) <- NULL
            s_startJ(outB) <- 1
            s_startI(outB) <- 1
            outB <- c_anon_worker(outB, input = paraList)
          }
        }
      }

      suppPatternA <- g_suppPattern(g_problemInstance(outA))
      suppPatternB <- g_suppPattern(g_problemInstance(outB))

      cbind(suppPatternA[commonCellIndices[[1]]], suppPatternB[commonCellIndices[[2]]])
      indOK <- f.checkCommonCells(suppPatternA, suppPatternB, commonCellIndices)
      if (indOK)
        runInd <- FALSE
      if (counter > paraList$maxIter) {
        runInd <- FALSE
        warning("iterative procedure did not converge! --> returning NULL")
        return(NULL)
      }
      counter <- counter + 1
    }
  }
  if (paraList$verbose) {
    cat("\n===> all common cells have the same anonymity-state in both tables after", counter, "iterations! [Finished]\n")
  }
  outA <- c_finalize(object = outA, input = paraList)
  outB <- c_finalize(object = outB, input = paraList)
  return(list(outObj1 = outA, outObj2 = outB))
}