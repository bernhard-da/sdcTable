#' @aliases get.problemInstance,problemInstance,character-method
#' @rdname get.problemInstance-method
setMethod(f='get.problemInstance', signature=c('problemInstance', 'character'),
	definition=function(object, type) {
		if ( !type %in% c('strID', 'nrVars', 'freq', 'w', 'numVars', 'sdcStatus', 
				'lb', 'ub', 'LPL', 'SPL', 'UPL', 'primSupps', 'secondSupps',
				'forcedCells', 'hasPrimSupps', 'hasSecondSupps', 'hasForcedCells',
				'weight', 'suppPattern') ) {
			stop("get.problemInstance:: argument 'type' is not valid!\n")
		}
		
		if ( type == 'strID' ) {
			return(object@strID) 	
		}		
		if ( type == 'nrVars' ) {
			return(length(get.problemInstance(object, type='strID'))) 	
		}
		if ( type == 'freq' ) {
			return(object@Freq)	
		}		
		if ( type == 'w' ) {
			return(object@w)	
		}			
		if ( type == 'numVars' ) {
			return(object@numVars)	
		}			
		if ( type == 'sdcStatus' ) {
			return(object@sdcStatus)	
		}		
		if ( type == 'lb' ) {
			return(object@lb) 	
		}
		if ( type == 'ub' ) {
			return(object@ub) 	
		}
		if ( type == 'LPL' ) {
			return(object@LPL) 	
		}
		if ( type == 'UPL' ) {
			return(object@UPL) 	
		}
		if ( type == 'SPL' ) {
			return(object@SPL) 	
		}	
		if ( type == 'primSupps' ) {
			return(which(get.problemInstance(object, type='sdcStatus')=="u"))	
		}	
		if ( type == 'secondSupps' ) {
			return(which(get.problemInstance(object, type='sdcStatus')=="x"))	
		}	
		if ( type == 'forcedCells' ) {
			return(which(get.problemInstance(object, type='sdcStatus')=="z"))	
		}	
		if ( type == 'hasPrimSupps' ) {
			return(length(get.problemInstance(object, type='primSupps')) > 0)	
		}		
		if ( type == 'hasSecondSupps' ) {
			return(length(get.problemInstance(object, type='secondSupps')) > 0)	
		}			
		if ( type == 'hasForcedCells' ) {
			return(length(get.problemInstance(object, type='forcedCells')) > 0)	
		}		
		if ( type == 'weight' ) {
			w <- get.problemInstance(object, type='w')
			if ( !is.null(w) ) {
				return(w)
			} else {
				return(get.problemInstance(object, type='freq'))
			}	
		}
		if ( type == 'suppPattern' ) {
			suppPattern <- rep(0, get.problemInstance(object, type='nrVars'))
			if ( get.problemInstance(object, type='hasPrimSupps') ) {
				suppPattern[get.problemInstance(object, type='primSupps')] <- 1
			}
			
			secondSupps <- get.problemInstance(object, type='secondSupps')
			if ( length(secondSupps) > 0 ) {
				suppPattern[secondSupps] <- 1
			}
			return(suppPattern)			
		}
	}
)

#' @aliases set.problemInstance,problemInstance,character,list-method
#' @rdname set.problemInstance-method
setMethod(f='set.problemInstance', signature=c('problemInstance', 'character', 'list'),
	definition=function(object, type, input) {
		index <- input[[1]]
		values <- input[[2]]
		if ( !type %in% c('lb', 'ub', 'LPL', 'UPL', 'SPL', 'sdcStatus') ) {
			stop("set.problemInstance:: check argument 'type'!\n" )
		}
		if ( !is.null(index) & length(values) != length(index) ) {
			stop("set.problemInstance:: arguments 'values' and 'index' differ in length!\n")
		} 
		if ( !all(index %in% 1:get.problemInstance(object, type='nrVars')) ) {
			stop("set.problemInstance:: argument 'index' does not fit to given problem!\n")
		} 	
		if ( type == 'lb' ) {
			object@lb[index] <- values
		}
		if ( type == 'ub' ) {
			object@ub[index] <- values
		}
		if ( type == 'LPL' ) {
			object@LPL[index] <- values
		}	
		if ( type == 'UPL' ) {
			object@UPL[index] <- values
		}
		if ( type == 'SPL' ) {
			object@SPL[index] <- values
		}	
		if ( type == 'sdcStatus' ) {
			sdcStatus <- get.problemInstance(object, type='sdcStatus')
			if ( length(input) > length(sdcStatus) ) {
				stop("set.problemInstance (type==sdcStatus):: length of 'sdcVec' must be <=",length(sdcStatus),"\n")
			}
			if ( is.null(index) ) {
				indexVec <- 1:length(sdcStatus)
				if ( length(values) != length(sdcStatus) ) {
					stop("set.problemInstance (type==sdcStatus):: length of 'values' must be ==",length(sdcStatus), "if 'index' is NULL!\n")
				}	
				object@sdcStatus[indexVec] <- values
			}
			if ( !is.null(index) ) {
				if ( !all(index %in% 1:length(sdcStatus)) ) {
					stop("set.problemInstance (type==sdcStatus):: elements of 'index' must be in 1:",length(sdcStatus),"!\n")
				}
				object@sdcStatus[index] <- values
			}		
		}			
		validObject(object)
		object		
	}
)

#' @aliases calc.problemInstance,problemInstance,character,list-method
#' @rdname calc.problemInstance-method
setMethod(f='calc.problemInstance', signature=c('problemInstance', 'character','list'),
	definition=function(object, type, input) {
		if ( !type %in% c('makeMasterProblem', 'isProtectedSolution') ) {
			stop("calc.problemInstance:: check argument 'type'!\n" )
		}
		if ( type == 'makeMasterProblem' ) {
			mProb <- NULL
			if ( get.problemInstance(object, type='hasPrimSupps') ) {
				objective <- get.problemInstance(object, type='weight')
				primSupps <- get.problemInstance(object, type='primSupps')
				nrVars <- get.problemInstance(object, type='nrVars')
				
				M <- init.simpleTriplet(type='simpleTriplet', input=list(mat=matrix(0, nrow=0, ncol=nrVars)))
				direction <- rep("==", get.simpleTriplet(M, type='nrRows', input=list()))
				rhs <- rep(1, get.simpleTriplet(M, type='nrRows', input=list()))
				
				# cells with sdcStatus=="z" must be published
				if ( get.problemInstance(object, type='hasForcedCells') ) {
					forcedCells <- get.problemInstance(object, type='forcedCells')
					if ( length(forcedCells) > 0 ) {
						for ( i in seq_along(forcedCells) ) {
							M <- calc.simpleTriplet(M, type='addRow', input=list(index=forcedCells[i], values=1))
						}
						direction <- c(direction, rep("==", length(forcedCells)))
						rhs <- c(rhs, rep(0, length(forcedCells)))
					}
				}
				types <- rep("C", nrVars)
				boundsLower <- list(ind=1:nrVars, val=rep(0, nrVars))
				boundsUpper <- list(ind=1:nrVars, val=rep(1, nrVars))
				
				if ( length(primSupps) > 0 ) {
					boundsLower$val[primSupps] <- 1
				}
				mProb <- new("linProb", 
					objective=objective,
					constraints=M,
					direction=direction,
					rhs=rhs,
					boundsLower=boundsLower,
					boundsUpper=boundsUpper,
					types=types)		
			}
			return(mProb) 			
		}

		if ( type == 'isProtectedSolution' ) {
			input1 <- input$input1 # ~ limitsDown
			input2 <- input$input2 # ~ limitsUp
			primSupps <- get.problemInstance(object, type='primSupps')
			if ( length(input1) != length(input2) ) {
				stop("parameters 'input1 (~limitsDown)' and 'input2 (~limitsUp)' differ in length!\n")
			}
			if ( length(input1) != length(primSupps) ) {
				stop("parameter 'limits.x' and length of primary suppressed cells differ!\n")
			}			
			
			protected <- TRUE	
			weights <- get.problemInstance(object, type='weight')[primSupps]
			
			limits <- list()
			limits$LPL <- get.problemInstance(object, type='LPL')[primSupps]
			limits$UPL <- get.problemInstance(object, type='UPL')[primSupps]
			limits$SPL <- get.problemInstance(object, type='SPL')[primSupps]
			
			if ( any(weights - input1 < limits[[1]]) == TRUE ) {
				protected <- FALSE
			}
			
			if ( any(input2 - weights  < limits[[2]]) == TRUE ) {
				protected <- FALSE
			}		
			
			if ( any(input2 - input1  < limits[[3]]) == TRUE ) {
				protected <- FALSE
			}			
			return(protected)			
		}		
	}
)
