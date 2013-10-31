######################
### global methods ###
######################
setMethod(f='as.matrix', signature='simpleTriplet',
	definition=function(x, ...) { 		
		M <- matrix(0, nrow=get.simpleTriplet(x, type='nrRows', input=list()), ncol=get.simpleTriplet(x, type='nrCols', input=list()))
		i.x <- get.simpleTriplet(x, type='rowInd', input=list())
		j.x <- get.simpleTriplet(x, type='colInd', input=list())
		v.x <- get.simpleTriplet(x, type='values', input=list())
		for ( i in 1:get.simpleTriplet(x, type='nrCells', input=list()) ) {
			M[i.x[i], j.x[i]] <- v.x[i]
		}
		return(M) 
	}
)

##############################################
### methods only for class 'simpleTriplet' ###
##############################################
#' @aliases get.simpleTriplet,simpleTriplet,character,list-method
#' @rdname get.simpleTriplet-method
setMethod(f='get.simpleTriplet', signature=c('simpleTriplet','character', 'list'),
	definition=function(object, type, input) { 
		if ( !type %in% c('rowInd', 'colInd', 'values', 'nrRows', 
				'nrCols', 'nrCells', 'duplicatedRows', 'transpose',
				'getRow', 'getCol') ) {
			stop("get.simpleTriplet:: argument 'type' is not valid!\n")
		}
	
		if ( type == 'rowInd' ) {
			return(object@i)
		}
		
		if ( type == 'colInd' ) {
			return(object@j)
		}
		
		if ( type == 'values' ) {
			return(object@v)
		}
		
		if ( type == 'nrRows' ) {
			return(object@nrRows)
		}	
		
		if ( type == 'nrCols' ) {
			return(object@nrCols)
		}	
		
		if ( type == 'nrCells' ) {
			return(length(object@v))
		}	
		
		if ( type == 'duplicatedRows' ) {
			i <- get.simpleTriplet(object, type='rowInd', input=list())
			j <- get.simpleTriplet(object, type='colInd', input=list())
			v <- get.simpleTriplet(object, type='values', input=list())
			len <- get.simpleTriplet(object, type='nrRows', input=list())	
			o <- order(i, j)
			y <- split(paste(j[o], v[o], sep = "\r"), i[o])
			tmp <- character(len)
			names(tmp) <- seq_along(tmp)
			tmp[names(y)] <- sapply(y, paste, collapse = "\r")
			dupRows <- which(duplicated(tmp))
			if ( length(dupRows) == 0 ) {
				dupRows <- NULL
			}
			return(dupRows)		
		}
		
		if ( type == 'transpose' ) {
			out <- init.simpleTriplet(type='simpleTriplet', input=list(mat=matrix(0, nrow=get.simpleTriplet(object, type='nrCols', input=list()), ncol=0)))
			for ( i in 1:get.simpleTriplet(object, type='nrRows', input=list()) ) {
				r <- get.simpleTriplet(object, type='getRow', input=list(i))
				out <- calc.simpleTriplet(out, type='addCol', input=list(index=get.simpleTriplet(r, type='colInd', input=list()), values=get.simpleTriplet(r, type='values', input=list())))
			}
			return(out)		
		}
		
		if ( type == 'getRow' ) {
			## if somebody specifies a vector of length > 1
			## the row with the first index is returned
			index <- input[[1]][1]
			if ( !index %in% 1:get.simpleTriplet(object, type='nrRows', input=list()) ) {
				stop("get.simpleTriplet:: parameter 'index' must be >=1 and <=",get.simpleTriplet(object, type='nrRows', input=list()),"!\n")			
			}	

			out <- NULL
			indI <- which(get.simpleTriplet(object, type='rowInd', input=list()) == index)
			if ( length(indI) > 0 ) {
				out <- new("simpleTriplet",
					i=rep(1, length(indI)),
					j=get.simpleTriplet(object, type='colInd', input=list())[indI],
					v=get.simpleTriplet(object, type='values', input=list())[indI],
					nrRows=1,
					nrCols=get.simpleTriplet(object, type='nrCols', input=list())					
				)			
			}
			return(out)
		}
		
		if ( type == 'getCol' ) {
			## if somebody specifies a vector of length > 1
			## the column with the first index is returned
			index <- input[[1]][1]
			if ( !index %in% 1:get.simpleTriplet(object, type='nrCols', input=list()) ) {
				stop("get.simpleTriplet:: parameter 'index' must be >=1 and <=",get.simpleTriplet(object, type='nrCols', input=list()),"!\n")			
			}	
			out <- NULL
			indJ <- which(get.simpleTriplet(object, type='colInd', input=list()) == index)
			if ( length(indJ) > 0 ) {
				out <- new("simpleTriplet",
					i=get.simpleTriplet(object, type='rowInd', input=list())[indJ],
					j=rep(1, length(indJ)),
					v=get.simpleTriplet(object, type='values', input=list())[indJ],
					nrRows=get.simpleTriplet(object, type='nrRows', input=list()),
					nrCols=1					
				)			
			}		
			return(out)
		}			
	}
)

#' @aliases calc.simpleTriplet,simpleTriplet,character,list-method
#' @rdname calc.simpleTriplet-method
setMethod(f='calc.simpleTriplet', signature=c('simpleTriplet', 'character', 'list'),
	definition=function(object, type, input) { 
		if ( !type %in% c('removeRow', 'removeCol', 'addRow', 'addCol', 
			'modifyRow', 'modifyCol', 'modifyCell', 'bind') ) {
			stop("calc.simpleTriplet:: check argument 'type'!\n")
		}	
		
		if ( type == 'removeRow' ) {
			index <- input[[1]]
			if ( !all(index %in% 1:get.simpleTriplet(object, type='nrRows', input=list())) ) {
				stop("calc.simpleTriplet:: check dimensions of parameter 'index'!\n")			
			}
			ind <- which(get.simpleTriplet(object, type='rowInd', input=list()) %in% index)
			if ( length(ind) > 0 ) {
				object@i <- get.simpleTriplet(object, type='rowInd', input=list())[-ind]
				object@j <- get.simpleTriplet(object, type='colInd', input=list())[-ind]
				object@v <- get.simpleTriplet(object, type='values', input=list())[-ind]
			}
			object@nrRows <- get.simpleTriplet(object, type='nrRows', input=list())-length(index)
			object@i <- rep(1:length(unique(get.simpleTriplet(object, type='rowInd', input=list()))), table(get.simpleTriplet(object, type='rowInd', input=list())))		
		}

		if ( type == 'removeCol' ) {
			index <- input[[1]]
			if ( !all(index %in% 1:get.simpleTriplet(object, type='nrCols', input=list())) ) {		
				stop("calc.simpleTriplet:: check dimensions of parameter 'index'!\n")			
			}		
			ind <- which(get.simpleTriplet(object, type='colInd', input=list()) %in% index)
			if ( length(ind) > 0 ) {
				object@i <- get.simpleTriplet(object, type='rowInd', input=list())[-ind]
				object@j <- get.simpleTriplet(object, type='colInd', input=list())[-ind]
				object@v <- get.simpleTriplet(object, type='values', input=list())[-ind]
			}
			object@nrCols <- get.simpleTriplet(object, type='nrCols', input=list())-length(index)
			object@j <- rep(1:length(unique(get.simpleTriplet(object, type='colInd', input=list()))), table(get.simpleTriplet(object, type='colInd', input=list())))
		}			

		if ( type == 'addRow' ) {
			index <- input[[1]]
			values <- input[[2]]
			if ( !all(index %in% 1:get.simpleTriplet(object, type='nrCols', input=list())) ) {
				stop("calc.simpleTriplet (type==addRow):: check dimensions of parameter 'index'!\n")
			}		
			if ( length(index) != length(values) ) {
				stop("calc.simpleTriplet (type==addRow):: dimensions of 'index' and 'values' do not match!\n")
			}
			rowInd <- get.simpleTriplet(object, type='nrRows', input=list())+1
			object@nrRows <- rowInd
			ind <- which(values != 0)
			nrAddedCells <- length(ind)
			if ( nrAddedCells > 0 ) {
				object@i <- c(get.simpleTriplet(object, type='rowInd', input=list()), rep(rowInd, nrAddedCells))
				object@j <- c(get.simpleTriplet(object, type='colInd', input=list()), index[ind])
				object@v <- c(get.simpleTriplet(object, type='values', input=list()), values[ind])
			}				
		}

		if ( type == 'addCol' ) {
			index <- input[[1]]
			values <- input[[2]]
			if ( !all(index %in% 1:get.simpleTriplet(object, type='nrRows', input=list())) ) {
				stop("calc.simpleTriplet (type==addCol):: check dimensions of parameter 'index'!\n")
			}		
			if ( length(index) != length(values) ) {
				stop("calc.simpleTriplet (type==addCol):: dimensions of 'index' and 'values' do not match!\n")
			}
			colInd <- get.simpleTriplet(object, type='nrCols', input=list())+1
			object@nrCols <- colInd
			ind <- which(values != 0)
			nrAddedCells <- length(ind)
			if ( nrAddedCells > 0 ) {
				object@i <- c(get.simpleTriplet(object, type='rowInd', input=list()), index[ind])
				object@j <- c(get.simpleTriplet(object, type='colInd', input=list()), rep(colInd, nrAddedCells))
				object@v <- c(get.simpleTriplet(object, type='values', input=list()), values[ind])
			}			
		}
		
		if ( type == 'modifyRow' ) {
			rowInd <- input[[1]]
			colInd <- input[[2]]
			values <- input[[3]]
			if ( length(rowInd) != 1 ) {
				stop("calc.simpleTriplet (type==modifyRow):: length of parameter 'rowInd' must equal 1!\n")
			}
			if ( !rowInd %in% 1:get.simpleTriplet(object, type='nrRows', input=list()) ) {
				stop("calc.simpleTriplet (type==modifyRow):: check dimensions of parameter 'rowInd'!\n")
			}		
			if ( !all(colInd %in% 1:get.simpleTriplet(object, type='nrCols', input=list())) ) {
				stop("calc.simpleTriplet (type==modifyRow):: check dimensions of parameter 'colInd'!\n")
			}
			if ( length(colInd) != length(values) ) {
				stop("calc.simpleTriplet (type==modifyRow):: dimensions of 'colInd' and 'values' do not match!\n")
			}
			ind <- which(get.simpleTriplet(object, type='rowInd', input=list()) %in% rowInd)
			if ( length(ind) == 0 ) {
				stop("calc.simpleTriplet (type==modifyRow):: no row to modify!\n")
			}		
			
			for ( j in seq_along(colInd) ) {
				object <- calc.simpleTriplet(object, type='modifyCell', input=list(rowInd, colInd[j], values[j]))	
			}	
		}

		if ( type == 'modifyCol' ) {
			rowInd <- input[[1]]
			colInd <- input[[2]]
			values <- input[[3]]
			
			if ( length(colInd) != 1 ) {
				stop("calc.simpleTriplet (type==modifyCol):: length of parameter 'colInd' must equal 1!\n")
			}
			if ( !all(rowInd %in% 1:get.simpleTriplet(object, type='nrRows', input=list())) ) {
				stop("calc.simpleTriplet (type==modifyCol):: check dimensions of parameter 'rowInd'!\n")
			}		
			if ( !colInd %in% 1:get.simpleTriplet(object, type='nrCols', input=list()) ) {
				stop("calc.simpleTriplet (type==modifyCol):: check dimensions of parameter 'colInd'!\n")
			}
			if ( length(rowInd) != length(values) ) {
				stop("calc.simpleTriplet (type==modifyCol):: dimensions of 'rowInd' and 'values' do not match!\n")
			}
			ind <- which(get.simpleTriplet(object, type='colInd', input=list()) %in% colInd)
			if ( length(ind) == 0 ) {
				stop("calc.simpleTriplet (type==modifyCol):: no column to modify!\n")
			}		
			for ( i in seq_along(rowInd) ) {
				object <- calc.simpleTriplet(object, type='modifyCell', input=list(rowInd[i], colInd, values[i]))	
			}			
		}											

		if ( type == 'modifyCell' ) {
			rowInd <- input[[1]]
			colInd <- input[[2]]
			values <- input[[3]]
		
			if ( any(length(values), length(rowInd), length(colInd) != 1) ) {
				stop("calc.simpleTriplet (type==modifyCell):: length of all arguments 'rowInd', 'colInd', 'values' must equal 1!\n")
			}
			if ( !all(colInd %in% 1:get.simpleTriplet(object, type='nrCols', input=list())) ) {		
				stop("calc.simpleTriplet (type==modifyCell):: check dimensions of parameter 'colInd'!\n")			
			}			
			if ( !all(rowInd %in% 1:get.simpleTriplet(object, type='nrRows', input=list())) ) {		
				stop("calc.simpleTriplet (type==modifyCell):: check dimensions of parameter 'rowInd'!\n")			
			}				
			ind <- which(get.simpleTriplet(object, type='rowInd', input=list()) == rowInd & get.simpleTriplet(object, type='colInd', input=list()) == colInd)
			if ( length(ind) == 1 & values != 0 ) {
				object@v[ind] <- values
			} else if ( length(ind)==1 & values == 0 ) {
				object@i <- get.simpleTriplet(object, type='rowInd', input=list())[-ind]
				object@j <- get.simpleTriplet(object, type='colInd', input=list())[-ind]
				object@v <- get.simpleTriplet(object, type='values', input=list())[-ind]				
			} else if ( length(ind)==0 & values!=0 ) {
				object@i <- c(get.simpleTriplet(object, type='rowInd', input=list()), rowInd)
				object@j <- c(get.simpleTriplet(object, type='colInd', input=list()), colInd)
				object@v <- c(get.simpleTriplet(object, type='values', input=list()), values)	
			}		
		}					

		if ( type == 'bind' ) {
			object1 <- object
			object2 <- input[[1]]
			bindRow <- input[[2]]
			
			if ( bindRow == TRUE ) {
				### "rbind"
				if ( get.simpleTriplet(object1, type='nrCols', input=list()) != get.simpleTriplet(object2, type='nrCols', input=list()) ) {
					stop("calc.simpleTriplet (type==bind):: nr of columns of 'object1' and 'object2' differ!\n")
				}
				out <- init.simpleTriplet(type='simpleTriplet', input=list(mat=matrix(0, nrow=get.simpleTriplet(object1, type='nrRows', input=list())+ get.simpleTriplet(object2, type='nrRows', input=list()), ncol=get.simpleTriplet(object1, type='nrCols', input=list()))))
				object2@i <- object2@i + get.simpleTriplet(object1, type='nrRows', input=list())
			} else {
				### "cbind"
				if ( get.simpleTriplet(object1, type='nrRows', input=list()) != get.simpleTriplet(object2, type='nrRows', input=list()) ) {
					stop("calc.simpleTriplet (type==bind):: nr of rows of 'object1' and 'object2' differ!\n")
				}	
				out <- init.simpleTriplet(type='simpleTriplet', input=list(mat=matrix(0, nrow=get.simpleTriplet(object1, type='nrRows', input=list()), ncol=get.simpleTriplet(object1, type='nrCols', input=list())+get.simpleTriplet(object2, type='nrCols', input=list()))))
				object2@j <- object2@j + get.simpleTriplet(object1, type='nrCols', input=list())
			}
			out@i <- c(object1@i, object2@i)
			out@j <- c(object1@j, object2@j)
			out@v <- c(object1@v, object2@v)
			object <- out
		}
		validObject(object)
		return(object)	
	}
)

#' @aliases init.simpleTriplet,character,list-method
#' @rdname init.simpleTriplet-method
setMethod(f='init.simpleTriplet', signature=c('character', 'list'),
	definition=function(type, input) { 
		if ( !type %in% c('simpleTriplet', 'simpleTripletDiag') ) {
			stop("init.simpleTriplet:: check argument 'type'!\n")
		}	

		if ( type == 'simpleTriplet' ) {
			matA <- input$mat
			dims <- dim(matA)
			v <- as.vector(t(matA))
			ind <- v!=0
			i <- rep(1:dims[1], each=dims[2])[ind==TRUE]
			j <- rep(1:dims[2], length=dims[1]*dims[2])[ind]
			v <- v[ind]			
			out <- new("simpleTriplet",
				i=i,
				j=j,
				v=v,
				nrRows=dims[1],	
				nrCols=dims[2]
			)	
		}

		if ( type == 'simpleTripletDiag' ) {
			nrRows <- input$nrRows
			negative <- input$negative
			i <- j <- 1:nrRows
			if ( negative ) {
				v <- rep(-1, nrRows)
			} else {
				v <- rep(1, nrRows)
			}
			out <- new("simpleTriplet",
				i=i,
				j=j,
				v=v,
				nrRows=nrRows,	
				nrCols=nrRows
			)
		}
		validObject(out)
		return(out)
	}
)	
