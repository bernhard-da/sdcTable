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
#' @param dimList a named list with each list element being either a data-frame or a link to a .csv-file containing the complete level-hierarchy of a dimensional variable using a top-to-bottom approach. The list names correspond to variable names that must exist in argument \code{data}. The level-hierarchy must be specified as follows:
#' \itemize{
#' \item list-element is a data-frame that must contain exactly 2 columns with the first column specifying levels and the second column holding variable-codes.
#' \itemize{
#' \item first column: a character vector specifying levels with each vector element being a string only containing of '@@'s from length 1 to n. If a vector element consists of \code{i}-chars, the corresponding code is of level \code{i}. The code '@@' (one character) equals the grand total (level=1).
#' \item second column: a character vector specifying level codes
#' }
#' \item list-element is full path to a .csv-file with two columns seperated by semicolons (;) having the same structure as the data.frame described above
#' }
#' @param dimVarInd numeric vector (or NULL) defining the column-indices of dimensional variables (defining the table) within argument \code{data} 
#' @param freqVarInd numeric vector (or NULL) defining the column-indices of a variable holding counts within argument \code{data} 
#' @param numVarInd numeric vector (or NULL) defining the column-indices of additional numeric variables available in argument \code{data} 
#' @param weightInd numeric vector of length 1 (or NULL) defining the column-index of a variable holding weights that should be used during as objective coefficients during the cut and branch algorithm to protect primary sensitive cells within argument \code{data} 
#' @param sampWeightInd numeric vector of length 1 (or NULL) defining the column-index of a variable holding sampling weights within argument \code{data} 
#' 
#' @return a \code{\link{sdcProblem-class}}-object
#' 
#' @examples
#' \dontrun{ 
#' # loading micro data
#' sp <- searchpaths()
#' fn <- paste(sp[grep("sdcTable", sp)], "/data/microData1.RData", sep="")
#' microData <- get(load(fn))
#' 
#' having a look at the data structure
#' str(microData)
#' 
#' # specify structure of hierarchical variable 'region'
#' # levels 'A' to 'D' sum up to a Total
#' dim.region <- data.frame(
#' 	levels=c('@@','@@@@','@@@@','@@@@','@@@@'),
#' 	codes=c('Total', 'A','B','C','D'),
#' 	stringsAsFactors=FALSE)
#' 	
#' # specify structure of hierarchical variable 'gender'
#' # levels 'male' and 'female' sum up to a Total
#' dim.gender <- data.frame(
#' 	levels=c('@@','@@@@','@@@@'),
#' 	codes=c('Total', 'male','female'),
#' 	stringsAsFactors=FALSE)
#' 
#' # create a list with each element being a data-frame containing information
#' # on a dimensional variables
#' dimList <- list(dim.region, dim.gender)
#' 
#' # name the list: 
#' # - first list-element: corresponds to variable 'region'
#' # - second list-element: corresponds to variable 'gender'
#' names(dimList) <- c('region', 'gender')
#' 
#' # specify the indices where dimensional variables are located 
#' # within the input data
#' 
#' # - variable 'region': first column
#' # - variable 'gender': second column
#' dimVarInd <- c(1,2) 
#' 
#' # no variables holding counts, numeric values, weights or sampling 
#' # weights are available in the input data
#' freqVarInd <- numVarInd <- weightInd <- sampWeightInd <- NULL
#' 
#' # we are dealing with micro data!
#' 
#' # creating an object of class \code{\link{sdcProblem-class}}
#' problem <- makeProblem(
#' 	data=microData, 
#' 	dimList=dimList, 
#' 	dimVarInd=dimVarInd, 
#' 	freqVarInd=freqVarInd, 
#' 	numVarInd=numVarInd, 
#' 	weightInd=weightInd,
#' 	sampWeightInd=sampWeightInd) 
#' 
#' what do we have?
#' print(class(problem))
#' }
#' @rdname makeProblem
#' @export makeProblem
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
makeProblem <- function(data, dimList, dimVarInd, freqVarInd=NULL, numVarInd=NULL, weightInd=NULL,sampWeightInd=NULL) {
	# returns an object of class 'sdcProblem'
	# 'doPrep()' is the old function 'newDimInfo()'
	# since it also recodes inputData eventually, it was renamed
	doPrep <- function(inputData, inputDims) {
		if ( any(sapply(inputDims, class) != "dimVar") ) {
			stop("Error: all elements of 'inputDims' must be of class 'dimVar'!\n")
		} 	
		if ( class(inputData) != "dataObj") {
			stop("Error: 'inputData' be of class 'dataObj'!\n")
		}
		
		varNames <- get.dataObj(inputData, type='varName')
		varNamesInDims <- sapply(1:length(dimList), function(x) { get.dimVar(dimList[[x]], type='varName') })
		
		if ( !all(varNamesInDims %in% varNames) ) {
			stop("makeProblem::doPrep() mismatch in variable names in 'inputData' and 'inputDims'!\n")
		}	
		
		rawData <- get.dataObj(inputData, type='rawData')
		
		# variable names in dataObj
		vNamesInData <- get.dataObj(inputData, type='varName')
		
		# vNames in inputDims
		vNamesInDimList <- sapply(1:length(inputDims), function(x) { get.dimVar(inputDims[[x]], type='varName') })
		
		# variables not used
		vNotUsed <- setdiff(vNamesInDimList, varNames)
		if ( length(vNotUsed) > 0 ) {
			removeIndex <- match(vNotUsed, vNamesInDimList)
			inputDims <- inputDims[-c(removeIndex)]
			vNamesInDimList <- sapply(1:length(inputDims), function(x) { get.dimVar(inputDims[[x]], type='varName') })
			
			if ( any(vNamesInDimList != varNames) ) {
				stop("Error: Matching failed!\n")
			} 
		}
		
		posIndex <- match(vNamesInData, vNamesInDimList)
		dimVarInd <- get.dataObj(inputData, type='dimVarInd')
		if ( length(posIndex) < 1 ) {
			stop("Error: matching of variable names failed. Please check 'inputData' and/or 'inputDims'!\n")
		} else {
			if ( any(is.na(posIndex)) ) {
				dimVarInd <- setdiff(dimVarInd, which(is.na(posIndex)))
				vNamesInData <- vNamesInData[dimVarInd]
				inputDims <- inputDims[na.omit(posIndex)]
			} else {
				# correct order
				inputDims <- inputDims[posIndex]			
			}
		}
		
		ss <- list()
		for ( i in seq_along(dimVarInd) ) {
			if ( !calc.dimVar(inputDims[[i]], type='hasDefaultCodes', input=rawData[[dimVarInd[i]]]) ) {
				dups <- get.dimVar(inputDims[[i]], type='dups')
				if ( length(dups) > 0 ) {
					dupsUp <- get.dimVar(inputDims[[i]], type='dupsUp')
					for ( k in length(dups):1 ) {
						ind <- which(rawData[[dimVarInd[i]]]==dups[k])
						if ( length(ind) > 0 ) {
							rawData[[dimVarInd[i]]][ind] <- dupsUp[k]
						}
					}
					inputData <- set.dataObj(inputData, type='rawData', input=list(rawData))
				}
				ss[[i]] <- calc.dimVar(inputDims[[i]], type='standardize', input=rawData[[dimVarInd[i]]])
			} else {
				ss[[i]] <- rawData[[dimVarInd[i]]]
			}
		}
		strID <- pasteStrVec(as.vector(unlist(ss)), length(posIndex))		
		
		info <- lapply(inputDims, function(x) {sum(get.dimVar(x, type='structure'))} )
		strInfo <- list()
		for ( i in 1:length(inputDims) ) {
			sumCur <- info[[i]]	
			if ( i == 1 )
				strInfo[[i]] <- c(1, sumCur)
			else 
				strInfo[[i]] <- c(1+max(strInfo[[c(i-1)]]), max(strInfo[[c(i-1)]])+sumCur)
		}	
		
		dimInfoObj <- new("dimInfo", 
				dimInfo=inputDims,
				strID=strID,
				strInfo=strInfo,
				vNames=vNamesInData,# because of ordering
				posIndex=dimVarInd # because dimVars are re-ordered according to input data!
		)	
		return(list(inputData=inputData, dimInfoObj=dimInfoObj))
	}	
	
	for ( i in seq_along(dimList) ) {
		dimList[[i]] <- init.dimVar(input=list(input=dimList[[i]], vName=names(dimList)[i])) 
	}
	
	## generate inputData from data
	inputData <- init.dataObj(input=list(inputData=data, dimVarInd=dimVarInd, freqVarInd=freqVarInd, numVarInd=numVarInd, weightInd=weightInd,sampWeightInd=sampWeightInd))
	
	## check if all variable names listed in inputDims exist in the 
	## specified dimensions of the input data
	varNames <- get.dataObj(inputData, type='varName')
	varNamesInDims <- sapply(1:length(dimList), function(x) { get.dimVar(dimList[[x]], type='varName') })
	
	if ( !all(varNamesInDims %in% varNames) ) {
		stop("makeProblem:: mismatch in variable names in 'inputData' and 'inputDims'!\n")
	}
	
	## calculate the dimInfoObj and eventually recode inputData
	## (eventually recode rawData slot of inputData if "rawData" contains "wrong" dups)
	out <- doPrep(inputData, dimList)
	
	## use output of doPrep() to calculate an object of class "sdcProblem"
	prob <- calc.multiple(type='calcFullProblem', input=list(objectA=out$inputData, objectB=out$dimInfoObj))
	prob
}

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
#' }
#' @param ... parameters used in the identification of primary sensitive cells. Parameters that can be modified|changed are:
#' \itemize{
#' \item \code{maxN}: numeric vector of length 1 used when applying the frequency rule. All cells having counts <= \code{maxN} are set as primary suppressed. The default value of \code{maxN} is 3.
#' \item \code{allowZeros}: logical vector of length 1 specifying if empty cells (count==0) should be considered sensitive when using the frequency rule. The default value of \code{allowZeros} is 'FALSE' so that empty cells are not considered primary sensitive by default.
#' \item \code{p}: numeric vector of length 1 specifying parameter \code{p} that is used when applying the p-percent rule with default value of 80.
#' \item \code{n}: numeric vector of length 1 specifying parameter \code{n} that is used when applying the nk-dominance rule. Parameter \code{n} is set to 2 by default.
#' \item \code{k}: numeric vector of length 1 specifying parameter \code{k} that is used when applying the nk-dominance rule. Parameter \code{n} is set to 85 by default.
#' \item \code{numVarInd}: numeric vector of length 1 specifying the index of the numerical variable that should be used to identify cells that are dominated by 2 (p-percent rule) or n (nk-dominance)-rule. If \code{type} is either 'nk' or 'p' it is mandatory to specify p.
#' }
#' @return a \code{\link{sdcProblem-class}} object
#'  
#' @examples
#' \dontrun{ 
#' # load micro data
#' sp <- searchpaths()
#' fn <- paste(sp[grep("sdcTable", sp)], "/data/microData1.RData", sep="")
#' microData <- get(load(fn))
#' 
#' # load problem (as it was created in the example in \code{\link{makeProblem}})
#' fn <- paste(sp[grep("sdcTable", sp)], "/data/problem.RData", sep="")
#' problem <- get(load(fn))
#' 
#' # we have a look at the table
#' print(table(microData))
#' 
#' # cell with region=='A' and gender=='female' has 2 units contributing to it
#' # this cell should be considered senstive!
#' problem <- primarySuppression(problem, type='freq', maxN=3)
#' 
#' # looking at anonymization states
#' print(table(getInfo(problem, type='sdcStatus')))
#' 
#' # we see that exactly one cell is primary suppressed (sdcStatus=='u') and
#' # the remaining cells are possible candidates for secondary suppression ('s')
#' }
#' @rdname primarySuppression
#' @export primarySuppression
#' @note the nk-dominance rule and the p-percent rule can only be applied if micro data have been used as input data to function \code{\link{makeProblem}}.
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
primarySuppression <- function(object, type, ...) {
	start.time <- proc.time()
	if ( !type %in% c('nk', 'freq', 'p') ) {
		stop("valid types are 'nk', 'freq' or 'p'!\n")
	}
	
	numVarsIndices <- get.dataObj(get.sdcProblem(object, type='dataObj'), type='numVarInd')
	paraList <- genParaObj(selection='control.primary', numVarIndices=numVarsIndices, ...)
	
	if ( type == "freq") {
		object <- calc.sdcProblem(object, type="rule.freq", input=paraList)
	}	
	
	if ( type == "nk" ) {
		if ( is.na(paraList$numVarInd) ) {
			stop("argument 'numVarInd' must be specified!\n")
		}
		object <- calc.sdcProblem(object, type="rule.nk", input=paraList)
	}
	
	if ( type == "p") {
		if ( is.na(paraList$numVarInd) ) {
			stop("argument 'numVarInd' must be specified!\n")
		}
		object <- calc.sdcProblem(object, type="rule.p", input=paraList)
	}	
	
	elapsed.time <- get.sdcProblem(object, type='elapsedTime') + (proc.time() - start.time)[3]
	object <- set.sdcProblem(object, type='elapsedTime', input=list(elapsed.time)) 
	return(object)	
}

#' protecting \code{\link{sdcProblem-class}} objects
#'
#' Function \code{\link{protectTable}} is used to protect primary sensitive table cells 
#' (that usually have been identified and set using 
#' \code{\link{primarySuppression}}). The function protects primary 
#' sensitive table cells according to the method that has been chosen and the 
#' parameters that have been set. Additional parameters that are used to control
#' the protection algorithm are set using parameter \code{...}.
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
#' \item \code{save}: logical vector of length 1 defining if temporary results should be saved in the current working directory (TRUE) or not (FALSE). Parameter \code{safe} defaults to 'FALSE' }
#' \item parameters used for HITAS|OPT procedures:
#' \itemize{
#' \item \code{solver}: character vector of length 1 defining the solver to be used. Currently available choices are limited to 'glpk'.
#' \item \code{timeLimit}: numeric vector of length 1 (or NULL) defining a time limit in minutes after which the cut and branch algorithm should stop and return a possible non-optimal solution. Parameter \code{safe} has a default value of 'NULL'
#' \item \code{maxVars}: a numeric vector of length 1 (or NULL) defining the maximum problem size in terms of decision variables for which an optimization should be tried. If the number of decision variables in the current problem are larger than parameter \code{maxVars}, only a possible non-optimal, heuristic solution is calculated. Parameter \code{safe} has a default value of 'NULL'
#' \item \code{fastSolution}: logical vector of length 1 defining if or if not the cut and branch algorithm will be started or if the possibly non-optimal heuristic solution is returned independent of parameter \code{maxVars}. Parameter \code{fastSolution} has a default value of 'FALSE'
#' \item \code{fixVariables}: logical vector of length 1 defining whether or not it should be tried to fix some variables to zero or one based on reduced costs early in the cut and branch algorithm. Parameter \code{fixVariables} has a default value of 'TRUE'
#' \item \code{approxPerc}: numeric vector of length 1 that defines a percentage for which a integer solution of the cut and branch algorithm is accepted as optimal with respect to the upper bound given by the (relaxed) solution of the master problem. Its default value is set to '10'}
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
#' }
#' 
#' @return an \code{\link{safeObj-class}} object
#' @examples
#' \dontrun{ 
#' # load problem (as it was created after performing primary suppression
#' # in the example of \code{\link{primarySuppression}})
#' sp <- searchpaths()
#' fn <- paste(sp[grep("sdcTable", sp)], "/data/problemWithSupps.RData", sep="")
#' problem <- get(load(fn))
#' 
#' # protect the table using the 'HITAS' algorithm with verbose output
#' protectedData <- protectTable(problem, method='HITAS', verbose=TRUE)
#' 
#' # showing a summary
#' summary(protectedData)
#' 
#' # looking at the final table with result suppression pattern
#' print(getInfo(protectedData, type='finalData'))
#' }
#' @rdname protectTable
#' @export protectTable
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
protectTable <- function(object, method, ...) {
	if ( !method %in% c('HITAS', 'OPT', 'HYPERCUBE', 'SIMPLEHEURISTIC') ) {
		stop("valid methods are 'SIMPLEHEURISTIC', 'HITAS', 'HYPERCUBE' or 'OPT'!\n")
	}
	
	paraList <- genParaObj(selection='control.secondary', method=method, ...)
	if ( method == 'SIMPLEHEURISTIC' ) {
		out <- performQuickSuppression(object, input=paraList)	
	} else {
		out <- calc.sdcProblem(object, type='anonWorker', input=paraList)
	}
	
	safeObj <- calc.sdcProblem(out, type='finalize', input=paraList)
	safeObj
}

#' query information from objects
#'
#' Function \code{\link{getInfo}} is used to query information from objects of class 
#' \code{\link{sdcProblem-class}}, \code{\link{problemInstance-class}} or \code{\link{safeObj-class}}
#'
#' @param object a \code{\link{sdcProblem-class}} object, \code{\link{problemInstance-class}} object or \code{\link{safeObj-class}} object. 
#' @param type a character vector of length 1 specifying the information which should be returned.
#' \itemize{
#' \item if argument \code{object} is of class \code{sdcProblem-class} or \code{\link{problemInstance-class}}, valid choices are:
#' \itemize{
#' \item \code{lb}: slot 'lb' of input \code{object} if it is of class \code{\link{problemInstance-class}} or this slot within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}}
#' \item \code{ub}: slot 'ub' of input \code{object} if it is of class \code{\link{problemInstance-class}} or this slot within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}}
#' \item \code{LPL}: slot 'LPL' of input \code{object} if it is of class \code{\link{problemInstance-class}} or this slot within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}}
#' \item \code{SPL}: slot 'SPL' of input \code{object} if it is of class \code{\link{problemInstance-class}} or this slot within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}}
#' \item \code{UPL}: slot 'UPL' of input \code{object} if it is of class \code{\link{problemInstance-class}} or this slot within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}}
#' \item \code{sdcStatus}:  slot 'sdcStatus' of input \code{object} if it is of class \code{\link{problemInstance-class}} or this slot within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}}
#' \item \code{freq}: slot 'freq' of input \code{object} if it is of class \code{\link{problemInstance-class}} or this slot within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}}
#' \item \code{strID}: slot 'strID' of input \code{object} if it is of class \code{\link{problemInstance-class}} or this slot within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}}
#' \item \code{numVars}: slot 'numVars' of input \code{object} if it is of class \code{\link{problemInstance-class}} or this slot within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}}
#' \item \code{w}: slot 'w' of input \code{object} if it is of class \code{\link{problemInstance-class}} or this slot within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}} }
#' \item if argument \code{object} is of class \code{\link{safeObj-class}}, valid choices are:
#' \itemize{
#' \item \code{finalData}: slot 'finalData' of input \code{object} of class \code{\link{safeObj-class}}
#' \item \code{nrNonDuplicatedCells}: slot 'nrNonDuplicatedCells' of input \code{object} of class \code{\link{safeObj-class}}
#' \item \code{nrPrimSupps}: slot 'nrPrimSupps' of input \code{object} of class \code{\link{safeObj-class}}
#' \item \code{nrSecondSupps}: slot 'nrSecondSupps' of input \code{object} of class \code{\link{safeObj-class}}
#' \item \code{nrPublishableCells}: slot 'nrPublishableCells' of input \code{object} of class \code{\link{safeObj-class}}
#' \item \code{suppMethod}: slot 'suppMethod' of input \code{object} of class \code{\link{safeObj-class}}}
#' }
#' 
#' @return manipulated data dependend on arguments \code{object} and \code{type}
#' 
#' @examples
#' \dontrun{ 
#' # load problem (as it was created in the example 
#' of \code{\link{makeProblem}})
#' sp <- searchpaths()
#' fn <- paste(sp[grep("sdcTable", sp)], "/data/problem.RData", sep="")
#' problem <- get(load(fn))
#' 
#' # problem is an object of class \code{\link{sdcProblem-class}}
#' print(class(problem))
#' 
#' for ( slot in c('lb','ub','LPL','SPL','UPL','sdcStatus', 
#' 		'freq', 'strID', 'numVars', 'w') ) {
#' cat('slot', slot,':\n')
#' 	print(getInfo(problem, type=slot))
#' }
#'
#' # extracting information for objects of class \code{\link{safeObj-class}} 
#' fn <- paste(sp[grep("sdcTable", sp)], "/data/protectedData.RData", sep="")
#' protectedData <- get(load(fn)) 
#' for ( slot in c('finalData', 'nrNonDuplicatedCells', 'nrPrimSupps', 
#' 		'nrSecondSupps', 'nrPublishableCells', 'suppMethod') ) {
#' 	cat('slot', slot,':\n')
#' 	print(getInfo(protectedData, type=slot))
#' }
#' }
#' @rdname getInfo
#' @export getInfo
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
getInfo <- function(object, type) {
	if ( !class(object) %in% c('sdcProblem', 'problemInstance', 'safeObj') ) {
		stop("getInfo:: argument 'object' must be of class 'sdcProblem', 'problemInstance' or 'safeObj'!\n")
	}
	
	if ( class(object) == 'safeObj' ) {
		if ( !type %in% c('finalData', 'nrNonDuplicatedCells', 'nrPrimSupps', 'nrSecondSupps', 'nrPublishableCells', 'suppMethod') ) {
			stop("getInfo:: type must be one of 'finalData', 'nrNonDuplicatedCells', 'nrPrimSupps', 'nrSecondSupps', 'nrPublishableCells' or 'suppMethod'!\n")
		}
		return(get.safeObj(object, type=type, input=list()))
	}
	else {
		if ( !type %in% c('lb','ub','LPL','SPL','UPL','sdcStatus', 'freq', 'strID', 'numVars', 'w') ) {
			stop("getInfo:: check argument 'type'!\n")
		}		
		if ( class(object) == 'sdcProblem' ) {
			pI <- get.sdcProblem(object, type='problemInstance')
		} else {
			pI <- object
		}	
		return(get.problemInstance(pI, type=type))		
	}
}

#' set information of \code{\link{sdcProblem-class}}- or \code{\link{problemInstance-class}} objects
#'
#' Function \code{\link{getInfo}} is used to query information from 
#' \code{\link{sdcProblem-class}}- or \code{\link{problemInstance-class}} objects 
#'
#' @param object an object of class \code{\link{sdcProblem-class}} or \code{\link{problemInstance-class}}
#' @param type a character vector of length 1 specifying the the information that should be changed or modified, valid choices are:
#' \itemize{
#' \item \code{lb}: slot 'lb' of input \code{object} if it is of class \code{\link{problemInstance-class}} or slot 'lb' within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}}
#' \item \code{ub}: slot 'ub' of input \code{object} if it is of class \code{\link{problemInstance-class}} or slot 'ub' within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}}
#' \item \code{LPL}: slot 'LPL' of input \code{object} if it is of class \code{\link{problemInstance-class}} or slot 'LPL' within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}}
#' \item \code{SPL}: slot 'SPL' of input \code{object} if it is of class \code{\link{problemInstance-class}} or slot 'SPL' within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}}
#' \item \code{UPL}: slot 'UPL' of input \code{object} if it is of class \code{\link{problemInstance-class}} or slot 'UPL' within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}}
#' \item \code{sdcStatus}:  slot 'sdcStatus' of input \code{object} if it is of class \code{\link{problemInstance-class}} or slot 'sdcStatus' within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}} }
#' @param index numeric vector defining cell-indices for which which values in a specified slot should be changed|modified
#' @param input numeric or character vector depending on argument \code{type} with its length matching the length of argument \code{index}
#' \itemize{
#' \item character vector if type matches 'sdcStatus'
#' \item a numeric vector if type matches 'lb', 'ub', 'LPL', 'SPL' or 'UPL'
#' }
#' 
#' @return a \code{\link{sdcProblem-class}}- or \code{\link{problemInstance-class}} object
#' 
#' @examples
#' \dontrun{ 
#' # load primary suppressed data (as created in the example 
#' of \code{\link{primarySuppression}})
#' sp <- searchpaths()
#' fn <- paste(sp[grep("sdcTable", sp)], "/data/problemWithSupps.RData", sep="")
#' problem <- get(load(fn))
#' 
#' # which is the overall total?
#' index.tot <- which.max(getInfo(problem, 'freq')
#' index.tot
#' 
#' # the anonymization state of the total is
#' print(getInfo(problem, type='sdcStatus')[index.tot])
#' 
#' # we want this cell to never be suppressed
#' problem <- setInfo(problem, type='sdcStatus', index=index.tot, input='z')
#' 
#' print(getInfo(problem, type='sdcStatus')[index.tot])
#' }
#' @rdname setInfo
#' @export setInfo
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setInfo <- function(object, type, index, input) {
	if ( !class(object) %in% c('sdcProblem', 'problemInstance') ) {
		stop("setInfo:: argument 'object' must be of class 'sdcProblem' or 'problemInstance'!\n")
	}
	
	if ( !type %in% c('lb','ub','LPL','SPL','UPL','sdcStatus') ) {
		stop("setInfo:: check argument 'type'!\n")
	}	
	
	if ( class(object) == "sdcProblem" ) {
		pI <- get.sdcProblem(object, type='problemInstance')
	} else {
		pI <- object
	}
	
	pI <- set.problemInstance(pI, type=type, input=list(index=index, values=input))
	
	if ( class(object) == "sdcProblem" ) {
		object <- set.sdcProblem(object, type='problemInstance', input=list(pI))
	} else {
		object <- pI
	}	
	object			
}

#' change anonymization status of a specific cell
#'
#' Function \code{\link{changeCellStatus}} allows to change|modify the anonymization state
#' of single table cells for objects ofs class \code{\link{sdcProblem-class}}.
#'
#' @param object an object of class \code{\link{sdcProblem-class}}
#' @param characteristics a character vector specifying characteristics of the table cell that should be identified for each dimensional variable defining the table 
#' @param varNames a character vector specifying variable names of dimensional variables defining the tables
#' @param rule character vector of length 1 specifying a valid anonymization code ('u', 'z', 'x', 's') to which the the cell under consideration should be set.
#' @param verbose logical vector of length 1 defining verbosity, defaults to 'FALSE'
#' 
#' @return a \code{\link{sdcProblem-class}} object
#' 
#' @examples
#' \dontrun{ 
#' # load primary suppressed data (as created in the example 
#' of \code{\link{primarySuppression}})
#' sp <- searchpaths()
#' fn <- paste(sp[grep("sdcTable", sp)], "/data/problemWithSupps.RData", sep="")
#' problem <- get(load(fn))
#' 
#' # we want to mark the cell region='D' and gender='male' primary sensitive
#' characteristics <- c('D', 'male')
#' varNames <- c('region', 'gender')
#' verbose <- TRUE
#' rule <- 'u' 
#' 
#' # looking at the distribution of anonymization states before...
#' print(table(getInfo(problem, 'sdcStatus')))
#' 
#' # setting the specific cell as primary sensitive
#' problem <- changeCellStatus(problem, characteristics, varNames, rule, verbose)
#' 
#' # having a second look at the anonymization states 
#' print(table(getInfo(problem, 'sdcStatus')))
#' }
#' @rdname changeCellStatus
#' @export changeCellStatus
#' @note Important: the \code{i}-th element of argument \code{characteristics} is uses as the desired characteristic for the dimensional variable specified at the \code{i}-th position of argument \code{varNames}!
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
changeCellStatus <- function(object, characteristics, varNames, rule, verbose=FALSE) {
	if ( class(object) != 'sdcProblem' ) {
		stop("changeCellStatus:: argument 'object' must be of class 'sdcProblem'!\n")
	}
	
	paraList <- list()
	paraList$names <- varNames
	paraList$codes <- characteristics
	paraList$verbose <- verbose
	
	cellID <- calc.sdcProblem(object, type='cellID', input=paraList)
	
	pI <- get.sdcProblem(object, type='problemInstance')
	pI <- set.problemInstance(pI, type='sdcStatus', input=list(index=cellID, values=rule))
	object <- set.sdcProblem(object, type='problemInstance', input=list(pI))
	
	if ( paraList$verbose ) {
		cat('--> The cell with ID=', cellID,'and Frequency',get.problemInstance(pI, type='freq')[cellID], 'has been set to', rule,'.\n')
	}
	object
}

#' query information for a specific cell in \code{\link{safeObj-class}} objects
#'
#' Function \code{\link{cellInfo}} is used to query information for a single table cell
#' for objects of class \code{\link{safeObj-class}}.
#'
#' @param object an object of class \code{\link{safeObj-class}}
#' @param characteristics a character vector specifying characteristics of the table cell that should be identified for each dimensional variable defining the table 
#' @param varNames a character vector specifying variable names of dimensional variables defining the tables
#' @param verbose logical vector of length 1 defining verbosity, defaults to 'FALSE'
#' 
#' @return a list containing the following calculated information
#' \itemize{
#' \item \code{cellID}: numeric vector of length 1 specifying the index of the cell within the final result dataset
#' \item \code{data}: a data.frame containing a single row with the index of the table cell of interest
#' \item \code{primSupp}: logical vector of length 1 that is 'TRUE' if the cell is a primary sensitive cell and 'FALSE' otherwise
#' \item \code{secondSupp}: logical vector of length 1 that is 'TRUE' if the cell is a secondary suppressed cell and 'FALSE' otherwise
#' }
#' 
#' @examples
#' \dontrun{ 
#' # load protected data (as created in the example 
#' of \code{\link{protectTable}})
#' sp <- searchpaths()
#' fn <- paste(sp[grep("sdcTable", sp)], "/data/protectedData.RData", sep="")
#' protectedData <- get(load(fn))
#' characteristics <- c('male', 'D')
#' varNames <- c('gender', 'region')
#' verbose <- TRUE
#' info <- cellInfo(protectedData, characteristics, varNames, verbose)
#' 
#' # show the info about this cell
#' str(info)
#' }
#' @rdname cellInfo
#' @export cellInfo
#' @note Important: the \code{i}-th element of argument \code{characteristics} is uses as the desired characteristic for the dimensional variable specified at the \code{i}-th position of argument \code{varNames}!
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
cellInfo <- function(object, characteristics, varNames, verbose=FALSE) {
	paraList <- list()
	paraList[[1]] <- varNames
	paraList[[2]] <- characteristics
	paraList[[3]] <- verbose
	get.safeObj(object, type='cellInfo', input=paraList)	
}

#' protect two \code{\link{sdcProblem-class}} objects that have common cells
#'
#' \code{\link{protectLinkedTables}} can be used to protect tables, that have
#' common cells. It is of course required that after the anonymization process
#' has finished, all common cells have the same anonymization state in both 
#' tables.
#'
#' @param objectA a \code{\link{sdcProblem-class}} object
#' @param objectB a \code{\link{sdcProblem-class}} object
#' @param commonCells a list object defining common cells in code{objectA} and \code{objectB}. For each variable that has one or more common codes in both tables, a list element needs to be specified. 
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
#' \item \code{HYPERCUBE}:
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
#' dim.region <- data.frame(h=c('@@','@@@@','@@@@'), l=c('Tot', 'R1','R2'))
#' 
#' # variable 'gender': exists in both datasets
#' dim.gender <- data.frame(h=c('@@','@@@@','@@@@'), l=c('Tot', 'm','f'))
#' 
#' # variable 'ecoOld': exists only in microDat1
#' dim.ecoOld <- data.frame(
#' 	h=c('@@','@@@@','@@@@@@','@@@@@@','@@@@','@@@@@@','@@@@@@'), 
#' 	l=c('Tot','A','Aa','Ab','B','Ba','Bb'))
#' 
#' # variable 'ecoNew': exists only in microDat2
#' dim.ecoNew <- data.frame(
#' 	h=c('@@','@@@@','@@@@@@','@@@@@@','@@@@@@','@@@@','@@@@@@','@@@@@@','@@@@@@'), 
#' 	l=c('Tot','C','Ca','Cb','Cc','D','Da','Db','Dc'))
#' 
#' # creating objects holding information on dimensions
#' dimList1 <- list(gender=dim.gender, ecoOld=dim.ecoOld)			
#' dimList2 <- list(region=dim.region, gender=dim.gender, ecoNew=dim.ecoNew)
#' 
#' # creating input objects for further processing. For details have a look at
#' # \code{\link{makeProblem}}.
#' problem1 <- makeProblem(data=microData1, dimList=dimList1, dimVarInd=c(1,2),  
#' 			numVarInd=3, isMicroData=TRUE)
#' problem2 <- makeProblem(data=microData2, dimList=dimList2, dimVarInd=c(1,2,3), 
#' 			numVarInd=4, isMicroData=TRUE)
#' 
#' # the cell specified by gender=='Tot' and ecoOld=='A'
#' # is one of the common cells! -> we mark it as primary suppression
#' problem1 <- changeCellStatus(problem1, characteristics=c('Tot', 'A'), 
#' 		varNames=c('gender','ecoOld'), rule='u', verbose=FALSE) 
#' 
#' # the cell specified by region=='Tot' and gender=='f' and ecoNew=='C'
#' # is one of the common cells! -> we mark it as primary suppression
#' problem2 <- changeCellStatus(problem2, characteristics=c('Tot', 'f', 'C'), 
#' 	varNames=c('region','gender', 'ecoNew'), rule='u', verbose=FALSE) 
#' 
#' # specifying input to define common cells
#' commonCells <- list()
#' 
#' # variable "gender"
#' commonCells$v.gender <- list()
#' commonCells$v.gender[[1]] <- 'gender' # variable name in 'problem1'
#' commonCells$v.gender[[2]] <- 'gender' # variable name in 'problem2'
#' # 'gender' has equal characteristics on both datasets -> keyword 'ALL'
#' commonCells$v.gender[[3]] <- 'ALL' 
#' 
#' # variable: ecoOld and ecoNew
#' commonCells$v.eco <- list()
#' commonCells$v.eco[[1]] <- 'ecoOld'	# variable name in 'problem1'
#' commonCells$v.eco[[2]] <- 'ecoNew'	# variable name in 'problem2'
#' 
#' # vector of common characteristics: A and B in variable 'ecoOld' in 'problem1'
#' commonCells$v.eco[[3]] <- c("A","B")	
#' # correspond to characteristics 'C' and 'D' in variable 'ecoNew' in 'problem2'
#' commonCells$v.eco[[4]] <- c("C","D")
#' 
#' # protect the linked data
#' result <- protectLinkedTables(problem1, problem2, commonCells, method='HITAS', verbose=TRUE)
#' 
#' # having a look at the results
#' result.tab1 <- result[[1]]
#' result.tab2 <- result[[2]]
#' summary(result.tab1)
#' summary(result.tab2)
#' }
#' @rdname protectLinkedTables
#' @export protectLinkedTables
#' @seealso \code{\link{protectTable}}
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
protectLinkedTables <- function(objectA, objectB, commonCells, method, ...) {
	f.calcCommonCellIndices <- function(input1, input2, commonCells) {
		pI1 <- get.sdcProblem(input1, type='problemInstance')
		pI2 <- get.sdcProblem(input2, type='problemInstance')
		
		commonInd1 <- 1:get.problemInstance(pI1, type='nrVars')
		commonInd2 <- 1:get.problemInstance(pI2, type='nrVars')
			
		dI1 <- get.sdcProblem(input1, type='dimInfo')
		dI2 <- get.sdcProblem(input2, type='dimInfo')
		
		strInfo1 <- get.dimInfo(dI1, type='strInfo')
		strInfo2 <- get.dimInfo(dI2, type='strInfo')
		
		vNames1 <- get.dimInfo(dI1, type='varName')
		vNames2 <- get.dimInfo(dI2, type='varName')
		
		varsUsed1 <- as.character(sapply(commonCells, function(x) x[[1]]))
		varsUsed2 <- as.character(sapply(commonCells, function(x) x[[2]]))

		varsUsed1 <- which(vNames1 %in% varsUsed1)
		varsUsed2 <- which(vNames2 %in% varsUsed2)
		
		varsNotUsed1 <- setdiff(1:length(vNames1), varsUsed1)
		varsNotUsed2 <- setdiff(1:length(vNames2), varsUsed2)

		### for each common variable -> get original labels from dI
		codesDefault1 <- lapply(1:length(strInfo1), function(x) { mySplit(get.problemInstance(pI1, type='strID'), strInfo1[[x]][1]:strInfo1[[x]][2]) } )
		codesDefault2 <- lapply(1:length(strInfo2), function(x) { mySplit(get.problemInstance(pI2, type='strID'), strInfo2[[x]][1]:strInfo2[[x]][2]) } )
		
		codesOrig1 <-  list()
		for ( i in 1:length(codesDefault1) ) {
			codesOrig1[[i]] <- calc.dimVar(object=get.dimInfo(dI1, type='dimInfo')[[i]], type='matchCodeOrig', input=codesDefault1[[i]] )
		}
		codesOrig2 <-  list()
		for ( i in 1:length(codesDefault2) ) {
			codesOrig2[[i]] <- calc.dimVar(object=get.dimInfo(dI2, type='dimInfo')[[i]], type='matchCodeOrig', input=codesDefault2[[i]] )
		}
		
		### find matching indices
		for ( i in 1:length(commonCells) ) {
			# it is not the same variable --> different characterisics
			if ( length(commonCells[[i]]) != 3 ) {
				commonInd1 <- setdiff(commonInd1, which(!codesOrig1[[varsUsed1[i]]] %in% commonCells[[i]][[3]]))			
				commonInd2 <- setdiff(commonInd2, which(!codesOrig2[[varsUsed2[i]]] %in% commonCells[[i]][[4]]))			
			}		
		}	

		# eliminate 'subtotals' from variables that are not used!
		if ( length(varsNotUsed1) > 0 ) {
			for ( i in varsNotUsed1 )  {
				subTotals <- get.dimVar(get.dimInfo(dI1, type='dimInfo')[[varsNotUsed1[i]]], type='codesOriginal')[get.dimVar(get.dimInfo(dI1, type='dimInfo')[[varsNotUsed1[i]]], type='codesMinimal')==FALSE]
				commonInd1 <- setdiff(commonInd1, which(codesOrig1[[varsNotUsed1[i]]] %in% subTotals))				
			}	
		}
		
		if ( length(varsNotUsed2) > 0 ) {
			for ( i in varsNotUsed2 )  {
				subTotals <- get.dimVar(get.dimInfo(dI2, type='dimInfo')[[varsNotUsed2[i]]], type='codesOriginal')[get.dimVar(get.dimInfo(dI2, type='dimInfo')[[varsNotUsed2[i]]], type='codesMinimal')==FALSE]
				commonInd2 <- setdiff(commonInd2, which(!codesOrig2[[varsNotUsed2[i]]] %in% subTotals))				
			}	
		}
		
		if ( length(commonInd1) != length(commonInd2) ) {
			stop("Error: length of common cell indices must be equal!\n")
		}
		
		if ( any(get.problemInstance(pI1, type='freq')[commonInd1] != get.problemInstance(pI2, type='freq')[commonInd2]) ) {
			stop("Error: common cells must have same values!\n")
		} 		
		return(list(commonInd1=commonInd1, commonInd2=commonInd2))
	}	

	f.checkCommonCells <- function(suppPattern1, suppPattern2, commonCellIndices) {
		indOK <- TRUE
		if ( any(suppPattern1[commonCellIndices[[1]]] != suppPattern2[commonCellIndices[[2]]]) ) 
			indOK <- FALSE
		return(indOK)
	}
	
	### arguments ###
	if ( !method %in% c('HITAS', 'OPT', 'HYPERCUBE') ) {
		stop("valid methods are 'HITAS', 'HYPERCUBE' or 'OPT'!\n")
	}

	paraList <- genParaObj(selection='control.secondary', method=method, ...)
	
	### first run
	outA <- calc.sdcProblem(objectA, type='anonWorker', input=paraList)
	outB <- calc.sdcProblem(objectB, type='anonWorker', input=paraList)
		
	pI.A <- get.sdcProblem(outA, type='problemInstance')	
	pI.B <- get.sdcProblem(outB, type='problemInstance')		
	# calc original primary suppressions
	
	origPrimSupp1Index <- get.problemInstance(pI.A, type='primSupps')
	origPrimSupp2Index <- get.problemInstance(pI.B, type='primSupps')
	
	# calculate commonCells:
	commonCellIndices <- f.calcCommonCellIndices(outA, outB, commonCells)	
	
	# suppression patterns after the first run
	suppPatternA <- get.problemInstance(pI.A, type='suppPattern')
	suppPatternB <- get.problemInstance(pI.B, type='suppPattern')
	
	indOK <- f.checkCommonCells(suppPatternA, suppPatternB, commonCellIndices)
	if ( indOK == FALSE ) {
		if ( paraList$verbose)
			cat("we have to start the iterative procedure!\n")
		runInd <- TRUE
		counter <- 1
		while ( runInd == TRUE ) {
			x <- cbind(suppPatternA[commonCellIndices[[1]]], suppPatternB[commonCellIndices[[2]]])
			index <- list()
			i1 <- which(x[,1] == 0 & x[,2]==1)
			i2 <- which(x[,1] == 1 & x[,2]==0)
			index[[1]] <- commonCellIndices[[1]][i1]
			index[[2]] <- commonCellIndices[[2]][i2]
			
			for ( j in 1:2 ) {
				if ( length(index[[j]]) > 0 ) {
					if ( j == 1 ) {
						pI.A <- get.sdcProblem(outA, type='problemInstance')
						pI.A <- set.problemInstance(pI.A, type='sdcStatus', input=list(index=index[[j]], values=rep("u", length(index[[j]]))))
						outA <- set.sdcProblem(outA, type='problemInstance', input=list(pI.A))
						outA <- set.sdcProblem(outA, type='indicesDealtWith', input=list(NULL))
						outA <- set.sdcProblem(outA, type='startJ', input=list(1))
						outA <- set.sdcProblem(outA, type='startI', input=list(1))
						outA <- calc.sdcProblem(outA, type='anonWorker', input=paraList)
					} else {
						pI.B <- get.sdcProblem(outB, type='problemInstance')
						pI.B <- set.problemInstance(pI.B, type='sdcStatus', input=list(index=index[[j]], values=rep("u", length(index[[j]]))))
						outB <- set.sdcProblem(outB, type='problemInstance', input=list(pI.B))
						outB <- set.sdcProblem(outB, type='indicesDealtWith', input=list(NULL))
						outB <- set.sdcProblem(outB, type='startJ', input=list(1))
						outB <- set.sdcProblem(outB, type='startI', input=list(1))
						outB <- calc.sdcProblem(outB, type='anonWorker', input=paraList)
					}					
				}				
			}
			
			suppPatternA <- get.problemInstance(get.sdcProblem(outA, type='problemInstance'), type='suppPattern')
			suppPatternB <- get.problemInstance(get.sdcProblem(outB, type='problemInstance'), type='suppPattern')	
						
			cbind(suppPatternA[commonCellIndices[[1]]], suppPatternB[commonCellIndices[[2]]])
			indOK <- f.checkCommonCells(suppPatternA, suppPatternB, commonCellIndices)
			if ( indOK == TRUE )
				runInd <- FALSE	
			if( counter > paraList$maxIter ) {
				runInd <- FALSE
				warning("iterative procedure did not converge! --> returning NULL")
				return(NULL)
			}				
			counter <- counter + 1			
		}		
	} else {
		if ( paraList$verbose)
			cat("\n===> all common cells have the same anonymity-state in both tables! [Finished]\n")
	}
	if ( paraList$verbose)
		cat("\n===> all common cells have the same anonymity-state in both tables after",counter,"iterations! [Finished]\n")
	
	outA <- calc.sdcProblem(outA, type='finalize', input=paraList)
	outB <- calc.sdcProblem(outB, type='finalize', input=paraList)
	return(list(outObj1=outA, outObj2=outB))
}
