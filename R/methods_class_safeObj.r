#' @aliases get.safeObj,safeObj,character,list-method
#' @rdname get.safeObj-method
setMethod(f='get.safeObj', signature=c('safeObj', 'character', 'list'),
	definition=function(object, type, input) {
		if ( !type %in% c('dimInfo', 'elapsedTime', 'finalData', 'nrNonDuplicatedCells', 
				'nrPrimSupps', 'nrSecondSupps', 'nrPublishableCells', 'suppMethod', 
				'cellInfo', 'cellID') ) {
			stop("get.safeObj:: argument 'type' is not valid!\n")
		}
		
		if ( type == 'dimInfo' ) {
			return(object@dimInfo)
		}
		if ( type == 'elapsedTime' ) {
			return(object@elapsedTime)
		}		
		if ( type == 'finalData' ) {
			return(object@finalData)
		}	
		if ( type == 'nrNonDuplicatedCells' ) {
			return(object@nrNonDuplicatedCells)
		}				
		if ( type == 'nrPrimSupps' ) {
			return(object@nrPrimSupps)
		}	
		if ( type == 'nrSecondSupps' ) {
			return(object@nrSecondSupps)
		}	
		if ( type == 'nrPublishableCells' ) {
			return(object@nrPublishableCells)
		}		
		if ( type == 'suppMethod' ) {
			return(object@suppMethod)
		}			
		
		if ( type == 'cellInfo' ) {
			cellID <- get.safeObj(object, type='cellID', input=input)
			
			primSupp <- secondSupp <- NULL
			finalData <- get.safeObj(object, type='finalData', input=list())
			cellStatus <- finalData[cellID, 'sdcStatus']
			
			if ( cellStatus == "u" ) {
				primSupp <- TRUE
				secondSupp <- FALSE
				if ( input[[3]] )
					cat ("The cell is a sensitive cell!\n")		
			}				
			if ( cellStatus == "s" ) {
				primSupp <- FALSE
				secondSupp <- FALSE				
				if ( input[[3]] )
					cat ("The cell can be published!\n")
			}
			if ( cellStatus == "z" ) {
				primSupp <- FALSE
				secondSupp <- FALSE		
				if ( input[[3]] )
					cat ("The cell will be enforced for publication!\n")			
			}	
			if ( cellStatus == "x" ) {
				primSupp <- FALSE
				secondSupp <- TRUE		
				if ( input[[3]] )
					cat ("The cell has been secondary suppressed!\n")			
			}				
			return(list(cellID=cellID, data=finalData[cellID,], primSupp=primSupp, secondSupp=secondSupp))
		}
		
		if ( type == 'cellID' ) {
			para.names <- input[[1]]
			para.codes <- input[[2]]
			para.verbose <- input[[3]]
			
			finalData <- get.safeObj(object, type='finalData', input=list())
			dimInfo <- get.safeObj(object, type='dimInfo', input=list())
			
			vNames <- get.dimInfo(dimInfo, type='varName')
			vIndex <- get.dimInfo(dimInfo, type='posIndex')		
			
			indexVar <- match(para.names, vNames)	
			
			if ( length(input) != 3 ) {
				stop("get.safeObj:: length of argument 'input' must equal 3!\n")
			}
			if ( length(para.names) != length(para.codes) ) {
				stop("get.safeObj:: check argument 'input'!\n")
			}
			if ( !all(para.names %in% vNames) ) {
				stop("get.safeObj:: check variable names in 'input[[1]]'!\n")
			}		
			if ( !is.logical(para.verbose) ) {
				stop("get.safeObj:: argument in 'input[[3]]' must be logical!\n")
			}		
			
			cellID <- 1:nrow(finalData)
			for ( i in seq_along(para.names) ) {
				cellID <- intersect(cellID, which(!is.na(match(as.character(finalData[,indexVar[i]]), para.codes[i]))))
			}
			if ( length(cellID) != 1) {
				stop("get.safeObj:: check argument 'input' -> 0 or > 1 cells identified!\n")
			}
			return(cellID)			
		}
	}
)

#' summarize \code{\link{safeObj-class}} objects
#' 
#' extract and show information stored in \code{\link{safeObj-class}} objects
#' 
#' @aliases summary,safeObj-method
#' @rdname summary-method
#' @export
#' @docType methods
setMethod(f='summary', signature='safeObj',
	definition=function (object, ...) {
		cat("\n######################################################\n")
		cat("### Summary of the result object of class 'safeObj' ###\n")
		cat("######################################################\n")
		cat(paste("--> The input data have been protected using algorithm ",get.safeObj(object, type='suppMethod', input=list()),".\n", sep=""))
		cat(paste("--> The algorithm ran for ",formatTime(get.safeObj(object, type='elapsedTime', input=list()))$time.str,".\n", sep=""))
		cat(paste("--> To protect ",get.safeObj(object, type='nrPrimSupps', input=list())," primary sensitive cells, ", get.safeObj(object, type='nrSecondSupps', input=list())," cells need to be additionally suppressed.\n", sep=""))
		cat(paste("--> A total of ",get.safeObj(object, type='nrPublishableCells', input=list())," cells may be published.\n", sep=""))
		
		finalData <- get.safeObj(object, type='finalData', input=list())
		nrCells <- nrow(finalData)
		nrNonDups <- get.safeObj(object, type='nrNonDuplicatedCells', input=list())
		if ( nrCells > nrNonDups ) {
			cat(paste("--> Duplicated cells: Only ",nrNonDups," table cells are unique, the remaining ",nrCells-nrNonDups," cells are duplicates.\n", sep=""))
		}
		cat("\n###################################\n")
		cat("### Structure of protected Data ###\n")
		cat("###################################\n")	
		print(str(finalData))
	}
)

#' show \code{\link{safeObj-class}} objects
#' 
#' extract and show information stored in \code{\link{safeObj-class}} objects
#' 
#' @aliases show,safeObj-method
#' @rdname show-method
#' @export
#' @docType methods
setMethod(f='show', signature='safeObj',
	definition=function (object) { print(str(object)) }
)
