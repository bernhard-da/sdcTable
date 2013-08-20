#' @aliases get.cutList,cutList,character-method
#' @rdname get.cutList-method
setMethod(f='get.cutList', signature=c('cutList', 'character'),
	definition=function(object, type) { 
		if ( !type %in% c('constraints', 'direction', 'rhs', 'nrConstraints') ) {
			stop("get.cutList:: argument 'type' is not valid!\n")
		}
		if ( type == 'constraints' ) {
			return(object@con)
		}
		if ( type == 'direction' ) {
			return(object@direction)
		}
		if ( type == 'rhs' ) {
			return(object@rhs)
		}				
		if ( type == 'nrConstraints' ) {
			return(length(get.cutList(object, type='rhs')))
		}	
	}
)

#' @aliases set.cutList,cutList,character,list-method
#' @rdname set.cutList-method
setMethod(f='set.cutList', signature=c('cutList', 'character', 'list'),
	definition=function(object, type, input) { 
		if ( !type %in% c('addCompleteConstraint', 'removeCompleteConstraint') ) {
			stop("set.cutList:: argument 'type' is not valid!\n")
		}
		
		input <- input[[1]]
		
		if ( type == 'addCompleteConstraint' ) {
			if ( get.simpleTriplet(get.cutList(object, type='constraints'), type='nrCols', input=list()) != get.simpleTriplet(get.cutList(input, type='constraints'), type='nrCols', input=list()) ) {
				stop("set.cutList:: nrCols of 'object' and 'input' differ!\n")
			}
			if ( get.cutList(input, type='nrConstraints') > 0 ) {
				con <- get.cutList(input, type='constraints')
				for ( k in 1:get.simpleTriplet(con, type='nrRows', input=list()) ) {
					x <- get.simpleTriplet(con, type='getRow', input=list(k))
					object@con <- calc.simpleTriplet(get.cutList(object, type='constraints'), type='addRow', input=list(index=get.simpleTriplet(x, type='colInd', input=list()), values=get.simpleTriplet(x, type='values', input=list())))			
				}
				object@direction <- c(get.cutList(object, type='direction'), get.cutList(input, type='direction'))
				object@rhs <- c(get.cutList(object, type='rhs'), get.cutList(input, type='rhs'))			
			}
		}
		
		if ( type == 'removeCompleteConstraint' ) {
			if ( !all(input %in% 1:length(get.cutList(object, type='rhs'))) ) {
				stop("elements of argument 'input' must be >=1 and <=",length(get.cutList(object, type='rhs')),"!\n")
			}
			object@con <- calc.simpleTriplet(get.cutList(object, type='constraints'), type='removeRow', input=list(input))			
			object@direction <- get.cutList(object, type='direction')[-input]
			object@rhs <- get.cutList(object, type='rhs')[-input]	
		}		
		
		validObject(object)
		return(object) 
	}
)

#' @aliases calc.cutList,cutList,character,list-method
#' @rdname calc.cutList-method
setMethod(f='calc.cutList', signature=c('cutList', 'character', 'list'),
	definition=function(object, type, input) {
		if ( !type %in% c('strengthen', 'checkViolation', 'bindTogether') ) {
			stop("calc.cutList:: argument 'type' is not valid!\n")
		}	
		
		if ( type == 'strengthen' ) {
			con <- get.cutList(object, type='constraints')
			nrRows <- get.simpleTriplet(con, type='nrRows', input=list())
			nrCols <- get.simpleTriplet(con, type='nrCols', input=list())
			rhs <- get.cutList(object, type='rhs')
			vals <- lapply(1:nrRows, function(x) { get.simpleTriplet(get.simpleTriplet(con, type='getRow', input=list(x)), type='values', input=list()) } )
			vals <- lapply(1:nrRows, function(x) { sapply(vals[[x]], function(y) { min(y, rhs[x]) } )   })
			colInds <- lapply(1:nrRows, function(x) { get.simpleTriplet(get.simpleTriplet(con, type='getRow', input=list(x)), type='colInd', input=list()) } )
			
			st <- init.simpleTriplet(type='simpleTriplet', input=list(mat=matrix(0, 0, nrCols)))
			lapply(1:nrRows, function(x) { st <<- calc.simpleTriplet(st, type='addRow', input=list(index=colInds[[x]], values=vals[[x]]))} )
			
			object@con <- st
			validObject(object)
			return(object)					
		}

		if ( type == 'checkViolation' ) {
			sol <- input[[1]]
			w <- input[[2]]
			con <- get.cutList(object, type='constraints')
			if ( get.simpleTriplet(con, type='nrCols', input=list()) != length(sol) ) {
				stop("calc.cutList (type==checkViolation):: length of argument 'sol' and number of columns of 'someCuts' differ!\n")
			}
			if ( length(sol) != length(w) ) {
				stop("calc.cutList (type==checkViolation):: length of argument 'sol' 'w' differ!\n")
			}		
			res <- rep(NA, get.simpleTriplet(con, type='nrRows', input=list()))
			rhs <- get.cutList(object, type='rhs')
			dir <- get.cutList(object, type='direction')
			for ( z in 1:get.simpleTriplet(con, type='nrRows', input=list()) ) {
				colInd <- get.simpleTriplet(get.simpleTriplet(con, type='getRow', input=list(z)), type='colInd', input=list())
				expr <- paste(sum(sol[colInd]*w[colInd]), dir[z], rhs[z])
				res[z] <- eval(parse(text=expr))		
			}
			return(any(res == FALSE))				
		}
		
		if ( type == 'bindTogether' ) {
			object1 <- object
			object2 <- input[[1]]
			x <- new("cutList")
			x@con <- calc.simpleTriplet(object=get.cutList(object1, type='constraints'), type='bind', input=list(get.cutList(object2, type='constraints'), bindRow=TRUE)) 
			x@direction <- c(get.cutList(object1, type='direction'), get.cutList(object2, type='direction'))
			x@rhs <- c(get.cutList(object1, type='rhs'), get.cutList(object2, type='rhs'))
			validObject(x)
			return(x)			
		}
	}
)

#' @aliases init.cutList,character,list-method
#' @rdname init.cutList-method
setMethod(f='init.cutList', signature=c('character', 'list'),
	definition=function(type, input) {
		if ( !type %in% c('empty', 'singleCut', 'multipleCuts') ) {
			stop("initialize.cutList:: argument 'type' is not valid!\n")
		}	
		
		x <- new("cutList")
		if ( type == 'empty' ) {
			x@con <- init.simpleTriplet(type='simpleTriplet', input=list(mat=matrix(0, nrow=0, ncol=input$nrCols)))
		}
		
		if ( type == 'singleCut' ) {
			x@con <- init.simpleTriplet(type='simpleTriplet', input=list(mat=matrix(input$vals, nrow=1, ncol=length(input$vals))))
			x@direction <- input$dir
			x@rhs <- input$rhs		
		}
		
		if ( type == 'multipleCuts' ) {
			x@con <- input$mat
			x@direction <- input$dir
			x@rhs <- input$rhs		
		}
		validObject(x)
		x				
	}
)
