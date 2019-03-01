#' create \code{\link{sdcProblem-class}}-objects
#'
#' Function \code{\link{makeProblem}} is used to create
#' \code{\link{sdcProblem-class}}-objects.
#'
#' @param data a data frame featuring at least one column for each desired
#' dimensional variable. Optionally the input data can feature variables
#' that contain information on cell counts, weights that should be used during
#' the cut and branch algorithm, additional numeric variables or variables that
#' hold information on sampling weights.
#' @param dimList a named list where the names refer to variable names in 
#' input \code{data}. Each list element can be one of: 
#' \itemize{
#' \item{\code{tree}: }{generated with \code{hier_*()} functions from the 
#' sdcHierarchies package}
#' \item{\code{data.frame}: }{a two column \code{data.frame} containing 
#' the complete level hierarchy of a dimensional variable using a 
#' top-to-bottom approach. The format of the \code{data.frame} is as follows:
#' \itemize{
#' \item first column: a character vector specifying levels with each vector 
#' element being a string only containing of '@@'s from length 1 to n. 
#' If a vector element consists of \code{i}-chars, the corresponding code 
#' is of level \code{i}. The code '@@' (one character) equals the grand 
#' total (level=1), the code `@@@@` (two characters) is of level 2 (directly
#' below the overall total).
#' \item second column: a character vector specifying level codes
#' }}
#' \item{\code{path}: }{absolute or relative path to a \code{.csv} file that
#' contains two columns seperated by semicolons (;) having the same structure 
#' as the \code{"@@;levelname"}-structure described above}
#' }
#' @param dimVarInd numeric vector (or NULL) defining the column-indices of 
#' dimensional variables (defining the table) within argument \code{data}. 
#' If \code{NULL}, the names of argument \code{dimList} are used to calculate 
#' the indices of the dimensional variables within \code{data} internally.
#' @param freqVarInd numeric vector (or NULL) defining the column-indices 
#' of a variable holding counts within argument \code{data}
#' @param numVarInd numeric vector (or NULL) defining the column-indices 
#' of additional numeric variables available in argument \code{data}
#' @param weightInd numeric vector of length 1 (or NULL) defining the 
#' column-index of a variable holding weights that should be used during 
#' as objective coefficients during the cut and branch algorithm to 
#' protect primary sensitive cells within argument \code{data}
#' @param sampWeightInd numeric vector of length 1 (or NULL) defining 
#' the column-index of a variable holding sampling weights within 
#' argument \code{data}
#'
#' @return a \code{\link{sdcProblem-class}}-object
#' @examples
#' # loading micro data
#' data("microData1", package="sdcTable")
#' microData <- microData1; rm(microData1)
#' # having a look at the data structure
#' str(microData)
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
#' dimList <- list(region=dim.region, gender=dim.gender)
#'
#' # third column containts a numeric variable
#' numVarInd <- 3
#'
#' # no variables holding counts, numeric values, weights or sampling
#' # weights are available in the input data
#' freqVarInd <- weightInd <- sampWeightInd <- NULL
#'
#' # creating an object of class \code{\link{sdcProblem-class}}
#' problem <- makeProblem(
#'  data=microData,
#'  dimList=dimList,
#'  freqVarInd=freqVarInd,
#'  numVarInd=numVarInd,
#'  weightInd=weightInd,
#'  sampWeightInd=sampWeightInd)
#'
#' # what do we have?
#' print(class(problem))
#'
#' # have a look at the data
#' sdcProb2df(problem, addDups=TRUE,
#'   addNumVars=TRUE, dimCodes="original")
#' @rdname makeProblem
#' @export makeProblem
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
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
    dimVarInd <- match(names(dimList), names(data))
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
      freqVarInd = freqVarInd,
      numVarInd = numVarInd,
      weightInd = weightInd,
      sampWeightInd = sampWeightInd
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
