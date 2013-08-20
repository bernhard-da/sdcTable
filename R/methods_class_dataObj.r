########################################
### methods only for class 'dataObj' ###
########################################
#' @aliases get.dataObj,dataObj,character-method
#' @rdname get.dataObj-method
setMethod(f='get.dataObj', signature=c('dataObj', 'character'),
	definition=function(object, type) { 
		if ( !type %in% c('rawData', 'dimVarInd', 'freqVarInd', 
				'numVarInd', 'weightVarInd', 'sampWeightInd', 
				'isMicroData', 'numVarNames', 'freqVarName', 'varName') ) {
			stop("get.dataObj:: argument 'type' is not valid!\n")
		}
		if ( type == 'rawData' ) {
			return(object@rawData) 	
		}
		if ( type == 'dimVarInd' ) {
			return(object@dimVarInd) 	
		}		
		if ( type == 'freqVarInd' ) {
			return(object@freqVarInd) 	
		}		
		if ( type == 'numVarInd' ) {
			return(object@numVarInd) 	
		}	
		if ( type == 'weightVarInd' ) {
			return(object@weightVarInd) 	
		}	
		if ( type == 'sampWeightInd' ) {
			return(object@sampWeightInd) 	
		}	
		if ( type == 'isMicroData' ) {
			return(object@isMicroData) 	
		}	
		if ( type == 'numVarNames' ) {
			return(names(get.dataObj(object, type='rawData'))[get.dataObj(object, type='numVarInd')])
		}		
		if ( type == 'freqVarName' ) {
			return(names(get.dataObj(object, type='rawData'))[get.dataObj(object, type='freqVarInd')])
		}	
		if ( type == 'varName' ) {
			return(names(get.dataObj(object, type='rawData'))[get.dataObj(object, type='dimVarInd')])
		}			
	}
)

#' @aliases set.dataObj,dataObj,character,listOrNULL-method
#' @rdname set.dataObj-method
setMethod(f='set.dataObj', signature=c('dataObj', 'character', 'listOrNULL'),
	definition=function(object, type, input) { 
		if ( !type %in% c('rawData') ) {
			stop("set.dataObj:: check argument 'type'!\n")
		}
		
		if ( type == 'rawData' ) {
			object@rawData <- input	
		}
		
		validObject(object)
		return(object)
	}
)

#' @aliases init.dataObj,list-method
#' @rdname init.dataObj-method
setMethod(f='init.dataObj', signature=c('list'),
	definition=function(input) { 
		inputData <- input$inputData
		dimVarInd <- input$dimVarInd
		freqVarInd <- input$freqVarInd
		numVarInd <- input$numVarInd
		weightInd <- input$weightInd
		sampWeightInd <- input$sampWeightInd
		isMicroData <- input$isMicroData
		
		rawData <- as.list(inputData)[c(dimVarInd, freqVarInd,numVarInd,weightInd,sampWeightInd)]
		if ( length(dimVarInd) >= 1 ) {
			dimVarInd <- 1:length(dimVarInd)	
			rawData[c(dimVarInd)] <- lapply(rawData[c(dimVarInd)], as.character)
		}
		last <- length(dimVarInd)	 
		if ( is.null(freqVarInd) ) {
			rawData$Freq <- rep(1, length(rawData[[1]]))
			freqVarInd <- length(rawData)
			if ( !is.null(sampWeightInd) ) {
				rawData[[freqVarInd]] <- rawData$Freq*rawData[[sampWeightInd]]
			}				
		} else {
			freqVarInd <- (last+1):((last)+length(freqVarInd))
			last <- last+1
			if ( !is.null(sampWeightInd) ) {
				rawData[[freqVarInd]] <- rawData[[freqVarInd]]*rawData[[sampWeightInd]]
			}				
		}
		if ( !is.null(numVarInd) ) {
			numVarInd <- (last+1):((last)+length(numVarInd))
		}	
		last <- last+length(numVarInd)	
		if ( !is.null(weightInd) ) {
			weightInd <- (last+1):((last)+length(weightInd))
		}
		last <- last+length(weightInd)	
		if ( !is.null(sampWeightInd) ) {
			sampWeightInd <- (last+1):((last)+length(sampWeightInd))
		}			
		out <- new("dataObj",
			rawData=rawData,
			dimVarInd=dimVarInd,
			freqVarInd=freqVarInd,
			numVarInd=numVarInd,
			weightVarInd=weightInd,
			sampWeightInd=sampWeightInd,
			isMicroData=isMicroData			
		)
		return(out)		
	}
)
