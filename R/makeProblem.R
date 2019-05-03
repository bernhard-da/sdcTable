#' Create a problem instance
#'
#' Function [makeProblem()] is used to create [sdcProblem-class] objects.
#'
#' @param data a data frame featuring at least one column for each desired
#' dimensional variable. Optionally the input data can feature variables
#' that contain information on cell counts, weights that should be used during
#' the cut and branch algorithm, additional numeric variables or variables that
#' hold information on sampling weights.
#' @param dimList a (named) list where the names refer to variable names in
#' input `data`. If the list is not named, it is required to specify argument
#' `dimVarInd`. Each list element can be one of:
#' - `tree`: generated with `hier_*()` functions from package `sdcHierarchies`
#' - `data.frame`: a two column `data.frame` containing the full hierarchy of
#' a dimensional variable using a top-to-bottom approach. The format of this
#' `data.frame` is as follows:
#'     * **first column:** a character vector specifying levels with each vector
#' element being a string only containing of `@@`s from length 1 to n.
#' If a vector element consists of `i`-chars, the corresponding code
#' is of level `i`. The code `@@` (one character) equals the grand
#' total (level=1), the code `@@@@` (two characters) is of level 2 (directly
#' below the overall total).
#'    * **second column:** a character vector specifying level codes
#' - `path`: absolute or relative path to a `.csv` file that
#' contains two columns seperated by semicolons (`;`) having the same structure
#' as the `"@@;levelname"`-structure described above
#' @param dimVarInd if `dimList` is a named list, this argument is
#' ignored (`NULL`). Else either a numeric or character vector
#' defining the column indices or names of dimensional variables
#' (specifying the table) within argument `data` are expected.
#' @param freqVarInd if not `NULL`, a scalar numeric or character vector
#' defining the column index or variable name of a variable holding counts
#' in `data`
#' @param numVarInd if not `NULL`, a numeric or character vector
#' defining the column indices or variable names of additional numeric
#' variables with respect to `data`
#' @param weightInd if not `NULL`, a scalar numeric or character vector
#' defining the column index or variable name holding costs within `data`
#' that should be used as objective coefficients when solving the secondary
#' cell suppression problem.
#' @param sampWeightInd if not `NULL`, a scalar numeric or character vector
#' defining the column index or variable name of a variable holding sampling
#' weights within `data`
#'
#' @return a [sdcProblem-class] object
#' @rdname makeProblem
#' @export
#' @author Bernhard Meindl
#' @md
#' @examples
#' # loading micro data
#' data("microData1", package="sdcTable")
#' # having a look at the data structure
#' str(microData1)
#'
#' # we can observe that we have a micro data set consisting
#' # of two spanning variables ('region' and 'gender') and one
#' # numeric variable ('val')
#'
#' # specify structure of hierarchical variable 'region'
#' # levels 'A' to 'D' sum up to a Total
#' dim.region <- data.frame(
#'  levels=c('@@','@@@@','@@@@','@@@@','@@@@'),
#'  codes=c('Total', 'A','B','C','D'),
#'  stringsAsFactors=FALSE)
#'
#' # specify structure of hierarchical variable 'gender'
#' # using create_node() and add_nodes() (see ?manage_hierarchies)
#' dim.gender <- hier_create(root = "Total", nodes = c("male", "female"))
#' hier_display(dim.gender)
#'
#' # create a named list with each element being a data-frame
#' # containing information on one dimensional variable and
#' # the names referring to variables in the input data
#' dimList <- list(region = dim.region, gender = dim.gender)
#'
#' # third column containts a numeric variable
#' numVarInd <- 3
#'
#' # no variables holding counts, numeric values, weights or sampling
#' # weights are available in the input data
#' # creating an problem instance using numeric indices
#' p1 <- makeProblem(
#'   data = microData1,
#'   dimList = dimList,
#'   numVarInd = 3 # third variable in `data`
#' )
#'
#' # using variable names is also possible
#' p2 <- makeProblem(
#'   data = microData1,
#'   dimList = dimList,
#'   numVarInd = "val"
#' )
#'
#' # what do we have?
#' print(class(p1))
#'
#' # have a look at the data
#' df1 <- sdcProb2df(p1, addDups = TRUE,
#'   addNumVars = TRUE, dimCodes = "original")
#' df2 <- sdcProb2df(p2, addDups=TRUE,
#'   addNumVars = TRUE, dimCodes = "original")
#' print(df1)
#'
#' identical(df1, df2)
makeProblem <- function(data,
                        dimList,
                        dimVarInd = NULL,
                        freqVarInd = NULL,
                        numVarInd = NULL,
                        weightInd = NULL,
                        sampWeightInd = NULL) {

  # returns an object of class 'sdcProblem'
  # 'doPrep()' is the old function 'newDimInfo()'
  # since it also recodes inputData eventually, it was renamed
  doPrep <- function(inputData, inputDims) {
    if (any(sapply(inputDims, class) != "dimVar")) {
      stop("Error: all elements of 'inputDims' must be of class 'dimVar'!\n")
    }
    if (class(inputData) != "dataObj") {
      stop("Error: 'inputData' be of class 'dataObj'!\n")
    }

    varNames <- g_var_name(inputData)
    varNamesInDims <- sapply(1:length(dimList), function(x) {
      g_varname(dimList[[x]])
    })

    if (!all(varNamesInDims %in% varNames)) {
      stop("makeProblem::doPrep() mismatch in variable names in 'inputData' and 'inputDims'!\n")
    }

    rawData <- g_raw_data(inputData)

    # variable names in dataObj
    vNamesInData <- g_var_name(inputData)

    # vNames in inputDims
    vNamesInDimList <- sapply(1:length(inputDims), function(x) {
      g_varname(inputDims[[x]])
    })

    # variables not used
    vNotUsed <- setdiff(vNamesInDimList, varNames)
    if (length(vNotUsed) > 0) {
      removeIndex <- match(vNotUsed, vNamesInDimList)
      inputDims <- inputDims[-c(removeIndex)]
      vNamesInDimList <- sapply(1:length(inputDims), function(x) {
        g_varname(inputDims[[x]])
      })

      if (any(vNamesInDimList != varNames)) {
        stop("Error: Matching failed!\n")
      }
    }

    posIndex <- match(vNamesInData, vNamesInDimList)
    dimVarInd <- g_dimvar_ind(inputData)
    if (length(posIndex) < 1) {
      stop("Error: matching of variable names failed. Please check 'inputData' and/or 'inputDims'!\n")
    } else {
      if (any(is.na(posIndex))) {
        dimVarInd <- setdiff(dimVarInd, which(is.na(posIndex)))
        vNamesInData <- vNamesInData[dimVarInd]
        inputDims <- inputDims[na.omit(posIndex)]
      } else {
        # correct order
        inputDims <- inputDims[posIndex]
      }
    }

    ss <- list()
    for (i in seq_along(dimVarInd)) {
      remove.vals <- FALSE
      remove_ind <- NULL
      if (!c_has_default_codes(inputDims[[i]], input = rawData[[dimVarInd[i]]])) {
        dups <- g_dups(inputDims[[i]])
        if (length(dups) > 0) {
          dupsUp <- g_dups_up(inputDims[[i]])
          for (k in length(dups):1) {
            ind <- which(rawData[[dimVarInd[i]]] == dups[k])
            if (length(ind) > 0) {
              if (length(which(rawData[[dimVarInd[i]]] == dupsUp[k])) > 0) {
                remove.vals <- TRUE
                remove_ind <- c(remove_ind, ind)
              } else {
                rawData[[dimVarInd[i]]][ind] <- dupsUp[k]
              }
            }
          }
          if (remove.vals) {
            rawData <- rawData[-unique(remove_ind)]
          }
          s_raw_data(inputData) <- list(rawData)
        }
        ss[[i]] <- c_standardize(inputDims[[i]], input = rawData[[dimVarInd[i]]])
      } else {
        ss[[i]] <- rawData[[dimVarInd[i]]]
      }
      # remove entries in ss[[1]...ss[[i-1]]
      if (remove.vals) {
        if (i > 1) {
          for (z in 1:(i - 1)) {
            ss[[z]] <- ss[[z]][-remove_ind]
          }
        }
      }
    }
    strID <- pasteStrVec(as.vector(unlist(ss)), length(posIndex))

    info <- lapply(inputDims, function(x) {
      sum(g_structure(x))
    })
    strInfo <- list()
    for (i in 1:length(inputDims)) {
      sumCur <- info[[i]]
      if (i == 1) {
        strInfo[[i]] <- c(1, sumCur)
      } else {
        strInfo[[i]] <- c(1 + max(strInfo[[c(i - 1)]]), max(strInfo[[c(i - 1)]]) + sumCur)
      }
    }

    dimInfoObj <- new(
      Class = "dimInfo",
      dimInfo = inputDims,
      strID = strID,
      strInfo = strInfo,
      vNames = vNamesInData,
      # because of ordering
      posIndex = dimVarInd # because dimVars are re-ordered according to input data!
    )
    return(list(inputData = inputData, dimInfoObj = dimInfoObj))
  }

  cn <- names(data)
  # cn: column names of data
  # make sure, indices are returned even when variable names
  # were provided
  .convert_to_ind <- function(cn, v) {
    if (!is.character(v)) {
      return(v)
    }
    ind <- match(v, cn)
    if (any(is.na(ind))) {
      stop("A non-existing variable name was provided.", call. = FALSE)
    }
    ind
  }

  reserved <- c("id", "freq", "Freq", "sdcStatus")
  if (any(reserved %in% names(dimList))) {
    stop("please do not use either 'id','freq','Freq' or 'sdcStatus' as names for dimensional variables!\n")
  }

  # check/calculate dimVarInd
  if (!all(names(dimList) %in% names(data))) {
    stop("For at least one dimensional variable specified in 'dimList', we do not have a corresponding variable in 'data'!\n")
  }

  # convert from tree- to standard format
  for (i in 1:length(dimList)) {
    if (inherits(dimList[[i]], "sdc_hierarchy")) {
      dimList[[i]] <- hier_convert(dimList[[i]], as = "df")
    }
  }

  if (is.null(dimVarInd)) {
    # we need to calculate dimVarInd from names in dimList
    dimVarInd <- .convert_to_ind(cn, v = names(dimList))
  } else {
    # we just need to check the names match
    if (!all(names(dimList) == names(data)[dimVarInd])) {
      stop("Names of dimensional variables specified in 'dimList' do not match with variables names in 'data' specified in 'dimVarInd'!\n")
    }
  }

  for (i in seq_along(dimList)) {
    dimList[[i]] <- init.dimVar(
      input = list(
        input = dimList[[i]],
        vName = names(dimList)[i]
      )
    )
  }

  ## generate inputData from data
  inputData <- init.dataObj(
    input = list(
      inputData = data,
      dimVarInd = dimVarInd,
      freqVarInd = .convert_to_ind(cn, v = freqVarInd),
      numVarInd = .convert_to_ind(cn, v = numVarInd),
      weightInd = .convert_to_ind(cn, v = weightInd),
      sampWeightInd = .convert_to_ind(cn, v = sampWeightInd)
    )
  )

  ## check if all variable names listed in inputDims exist in the
  ## specified dimensions of the input data
  varNames <- g_var_name(inputData)
  varNamesInDims <- sapply(1:length(dimList), function(x) {
    g_varname(dimList[[x]])
  })

  if (!all(varNamesInDims %in% varNames)) {
    stop("makeProblem:: mismatch in variable names in 'inputData' and 'inputDims'!\n")
  }

  ## calculate the dimInfoObj and eventually recode inputData
  ## (eventually recode rawData slot of inputData if "rawData" contains "wrong" dups)
  out <- doPrep(inputData, dimList)

  ## compute the full sdcProblem object
  c_calc_full_prob(
    input = list(
      objectA = out$inputData,
      objectB = out$dimInfoObj
    )
  )
}
