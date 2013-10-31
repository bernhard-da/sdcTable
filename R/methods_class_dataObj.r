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
			object@rawData <- input[[1]]
		}
		
		validObject(object)
		return(object)
	}
)

#' @aliases init.dataObj,list-method
#' @rdname init.dataObj-method
setMethod(f='init.dataObj', signature=c('list'),
	definition=function(input) { 
		freq <- N <- .N <- NULL
		inputData <- input$inputData
		dimVarInd <- input$dimVarInd
		freqVarInd <- input$freqVarInd
		numVarInd <- input$numVarInd
		weightInd <- input$weightInd
		sampWeightInd <- input$sampWeightInd
		
		isMicroData <- FALSE
			
		## aggregate data, use data.table
		datO <- data.table(inputData, key=colnames(inputData)[dimVarInd])
		rawData <- datO[,.N, by=key(datO)]
		
		if ( any(rawData$N != 1)  ) {
			isMicroData <- TRUE
			if ( is.null(freqVarInd) ) {
				rawData[, freq:=as.numeric(rawData$N)]
				freqVarInd <- which(colnames(rawData)=="freq")
			} else {
				f <- datO[,list(freq=sum(get(colnames(datO)[freqVarInd]))), by=key(datO)]$freq
				rawData[,freq:=as.numeric(f)]
			}			
		} else {
			# data already aggregated
			if ( is.null(freqVarInd) ) {
				rawData[, freq:=as.numeric(rawData$N)]
				freqVarInd <- which(colnames(rawData)=="freq")
			} else {
				rawData[, freq:=as.numeric(datO[,get(colnames(datO)[freqVarInd])])]
			}					
		}		
		
		rawData[,N:=NULL]
		
		if ( !is.null(sampWeightInd) ) {
			sw <- datO[,list(sw=sum(get(colnames(datO)[sampWeightInd]))), by=key(datO)]$sw
			rawData[, colnames(datO)[sampWeightInd]:=sw]
			sampWeightInd <- ncol(rawData)
			rawData[, freq:=as.numeric(sw*rawData$freq)]
		}	
		
		## numvars
		if ( !is.null(numVarInd) ) {
			c.start <- ncol(rawData)
			cols <- colnames(datO)[numVarInd]
			for ( j in 1:length(cols)) {
				v <- datO[,j=sum(get(cols[j])), by=key(datO)]$V1
				rawData[,j=cols[j]:=as.numeric(v)]
			}		
			numVarInd <- (c.start+1):ncol(rawData)			
		}	
		
		## weight var
		if ( !is.null(weightInd) ) {
			w <- datO[,list(w=sum(get(colnames(datO)[weightInd]))), by=key(datO)]$w
			rawData[, colnames(datO)[weightInd]:=as.numeric(w)]
			weightVarInd <- ncol(rawData)
		}
		
		## do not use factors
		cols <- colnames(rawData)[dimVarInd]
		for ( j in 1:length(cols)) {
			set(rawData,j=j,value=as.character(rawData[[cols[j]]]))
		}
		rm(datO)

		setkeyv(rawData, colnames(rawData)[dimVarInd])
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
