#' @aliases calc.multiple,character,list-method
#' @rdname calc.multiple-method
setMethod(f='calc.multiple', signature=c('character', 'list'),
	definition=function(type, input) {
		.SD <- ID <- NULL
		if (!type %in% c('makePartitions', 'genMatMFull', 
				'makeAttackerProblem', 'calcFullProblem') ) {
			stop("calc.multiple:: argument 'type' is not valid!\n")
		}	

		if ( type == 'makePartitions' ) {
			pI <- input$objectA
			dimInfoObj <- input$objectB
			dimInfo <- get.dimInfo(dimInfoObj, type='dimInfo')
			strIDs <- get.problemInstance(pI, type='strID')
			
			## create classes and groups
			tmpDat <- expand.grid(lapply(1:length(dimInfo), function(x) { 1:get.dimVar(dimInfo[[x]], type='nrLevels') } ))
			groups <- apply(tmpDat, 1, function(x) { paste(x, collapse="-")})
			classes <- apply(tmpDat, 1, sum)
			sortOrder <- order(classes)
			classes <- classes[sortOrder]
			classesUnique <- unique(classes)
			groups <- groups[sortOrder]
			splitGroups <- split(groups, classes)
			
			## create tables for all classes and groups 
			final <- list()
			final$groups <- as.list(groups)
			final$indices <- list()
			for ( i in 1:length(groups) ) {
				final$indices[[i]] <- list()
				levs <- as.integer(unlist(sapply(groups[[i]], strsplit, "-")))	
				
				res <- list()
				for ( z in 1:length(dimInfo) ) {
					res[[z]] <- list()
					index <- which(get.dimVar(dimInfo[[z]], type='levels') %in% c(levs[z], levs[z]-1))
					codesDefault <- get.dimVar(dimInfo[[z]], type='codesDefault')[index]
					if ( levs[z] == 1 ) {
						res[[z]] <- codesDefault
					} else {
						levOrig <- get.dimVar(dimInfo[[z]], type='levels')[index]
						diffs <- c(0,diff(levOrig))
						checkInd <- which(diffs == 1)-1
						out <- data.frame(index=index, levOrig=levOrig, codesDefault=codesDefault, ind=NA)
						out$ind[checkInd] <- 1
						
						checkInd <- c(checkInd, length(index))
						splitVec <- rep(0, length(index))
						for ( j in 2:length(checkInd) ) {
							if ( j < length(checkInd) ) {
								splitVec[checkInd[j-1]:(checkInd[j]-1)] <- j-1
							} else {
								splitVec[checkInd[j-1]:(checkInd[j])] <- j-1
							}
						}
						spl <- split(index, splitVec)
						counter <- 1
						for ( k in 1:length(spl) ) {
							rowInd <- match(spl[[k]], out$index)
							tmp <- out[rowInd,]
							if ( any(tmp[,"levOrig"]==levs[z]) ) {					
								tmp <- tmp[1:(max(which(tmp$levOrig==levs[z]))),]
								res[[z]][[length(res[[z]])+1]] <- sort(unique(as.character(tmp$codesDefault)))
							}
						}	
					}
				}
				final$indices[[i]] <- list()
				combs <- expand.grid(lapply(1:length(res), function(x) { 1:length(res[[x]]) } ))
				for ( m in 1:nrow(combs) ) {
					final.strIDs <- pasteStrVec(expand(lapply(1:ncol(combs), function(x) { res[[x]][[combs[m,x]]] })), ncol(combs))
					final$indices[[i]][[m]] <- which(strIDs %in% final.strIDs)
				}			
			}
			final$nrGroups <- length(groups)
			final$nrTables <- sum(sapply(1:final$nrGroups, function(x) { length(final$indices[[x]]) } ))				
			return(final)		
		}		
		
		if ( type == 'genMatMFull' ) {
			x <- input$objectA
			y <- input$objectB
			
			levelObj <- get.dimInfo(y, type='dimInfo')
			strID <- get.problemInstance(x, type='strID')
			nrVars <- length(levelObj)
			nrCells <- get.problemInstance(x, type='nrVars')
			freqs <- get.problemInstance(x, type='freq')
			
			constraintM <- init.simpleTriplet(type='simpleTriplet', input=list(mat=matrix(0, nrow=0, ncol=nrCells)))
			for ( i in 1:nrVars ) {
				lO <- levelObj[[i]]
				keepList <- lapply(get.dimInfo(y, type='strInfo')[-i], function(k) { seq(k[1], k[2]) } )
				keepList2 <- lapply(get.dimInfo(y, type='strInfo')[i], function(k) { seq(k[1], k[2]) } )
				f1 <- f2 <- mySplitIndicesList(strID, keepList2)
				
				if ( nrVars > 1 ) { 
					f1 <- mySplitIndicesList(strID, keepList)
				} 
				
				dimlO <- get.dimVar(lO, type='dims')
				if ( length(unique(f2)) != 1 ) {
					dimInd <- sapply(1:length(dimlO), function(x) { identical( sort(unique(f2)), dimlO[[x]]) } )
					if ( sum(dimInd) == 0 ) {
						for ( j in 1:length(get.dimVar(lO, type='dims')) ) {
							splitInd <- which(f2 %in% get.dimVar(lO, type='dims')[[j]])
							spl <- split(splitInd, f1[splitInd])
							for ( z in 1:length(spl) ) {
								ind <- rep(1,length(spl[[z]])) 
								ind[which.max(freqs[spl[[z]]])] <- -1
								if ( !is.zero(sum(freqs[spl[[z]]]*ind)) ) {
									stop("something went wrong!\n")
								}
								constraintM <- calc.simpleTriplet(constraintM, type='addRow', input=list(index=spl[[z]], values=ind))
							}
						}				
					} else {
						splitInd <- which(f2 %in% get.dimVar(lO, type='dims')[[which(dimInd==TRUE)]])
						## only 1 dimension
						if ( nrVars > 1 ) {
							spl <- split(splitInd, f1[splitInd])	
						} else {
							spl <- split(splitInd, rep(1, length(splitInd)))
						}
						
						for ( z in 1:length(spl) ) {
							ind <- rep(1,length(spl[[z]])) 
							ind[which.max(freqs[spl[[z]]])] <- -1
							if ( !is.zero(sum(freqs[spl[[z]]]*ind)) ) {
								stop("something went wrong! (z=",z," und names(spl)[z]='",names(spl)[z],")\n")
							}
							constraintM <- calc.simpleTriplet(constraintM, type='addRow', input=list(index=spl[[z]], values=ind))
						}				
					}				
				}
			}
			return(constraintM)		
		}	
	
		if ( type == 'makeAttackerProblem' ) {
			x <- input$objectA
			y <- input$objectB
			nrVars <- get.problemInstance(x, type='nrVars')	
			A <- calc.multiple(type='genMatMFull', input=list(objectA=x, objectB=y))
			
			## calculating (logical) constraints for the master problem ##
			# idea: for each constraint at least 2 suppressions must 
			# exist if one xi != 0! (http://www.eia.doe.gov/ices2/missing_papers.pdf)
			newCutsMaster <- init.cutList(type='empty', input=list(nrCols=nrVars))
			#xx <- lapply(1:get.simpleTriplet(A, type='nrRows', input=list()), function(x) {
			#	cols <- get.simpleTriplet(get.simpleTriplet(A, type='getRow', input=list(x)), type='colInd')
			#	v <- rep(0, nrVars)
			#	v[cols] <- c(1, rep(-1, length(cols)))
			#	newCutsMaster <<- set.cutList(newCutsMaster, type='addCompleteConstraint', input=list(init.cutList(type='singleCut', input=list(vals=v, dir="<=", rhs=0))))
			#})
			################################################################
			
			nrConstraints <- get.simpleTriplet(A, type='nrRows', input=list())
			objective <- rep(0, length=2*nrVars+nrConstraints)	
			z1 <- init.simpleTriplet(type='simpleTripletDiag', input=list(nrRows=nrVars, negative=FALSE))
			z2 <- init.simpleTriplet(type='simpleTripletDiag', input=list(nrRows=nrVars, negative=TRUE))
			z <- calc.simpleTriplet(object=z1, type='bind', input=list(z2, bindRow=FALSE))
			A <- calc.simpleTriplet(object=z, type='bind', input=list(get.simpleTriplet(A, type='transpose', input=list()), bindRow=FALSE))
			direction <- rep("==", get.simpleTriplet(A, type='nrRows', input=list()))
			rhs <- rep(0, get.simpleTriplet(A, type='nrRows', input=list()))
			
			types <- rep("C", get.simpleTriplet(A, type='nrCols', input=list()))
			boundsLower <- list(ind=1:get.simpleTriplet(A, type='nrCols', input=list()), val=c(rep(0, 2*nrVars), rep(-Inf, nrConstraints)))
			boundsUpper <- list(ind=1:get.simpleTriplet(A, type='nrCols', input=list()), val=c(rep(Inf, 2*nrVars), rep(Inf,  nrConstraints)))
			
			aProb <- new("linProb", 
				objective=objective,
				constraints=A,
				direction=direction,
				rhs=rhs,
				boundsLower=boundsLower,
				boundsUpper=boundsUpper,
				types=types)
			return(list(aProb=aProb, newCutsMaster=newCutsMaster))		
		}		
		
		if ( type == 'calcFullProblem' ) {
			x <- input$objectA
			y <- input$objectB
			time.start <- proc.time()
			rawData <- get.dataObj(x, type='rawData')
			dimObj <- get.dimInfo(y, type='dimInfo')
			ind.dimvars <- get.dataObj(x, type='dimVarInd')
			ind.freq <- get.dataObj(x, type='freqVarInd')
			
			## no need to aggregate data
			## aggregation already done in 'init.dataObj()'
			
			codes <- list(); length(codes) <- length(ind.dimvars)
			for ( i in 1:length(codes) ) {
				codes[[i]] <- rawData[[ind.dimvars[i]]]
				cDefault <- get.dimVar(dimObj[[i]], type='codesDefault')
				cOriginal <- get.dimVar(dimObj[[i]], type='codesOriginal')	
				cOriginalDups <- get.dimVar(dimObj[[i]], type='dups')
				cOriginalDupsUp <- get.dimVar(dimObj[[i]], type='dupsUp')
				if ( all(codes[[i]] %in% c(cOriginal, cOriginalDups)) ) {
					mInd1 <- match(codes[[i]], cOriginalDups)
					mInd2 <- which(!is.na(mInd1))
					if ( length(mInd2) > 0 ) {
						codes[[i]][mInd2] <- cOriginalDupsUp[mInd1[mInd2]]
					}
					codes[[i]] <- calc.dimVar(object=dimObj[[i]], type='matchCodeDefault', input=rawData[[ind.dimvars[i]]])
				} else if ( all(codes[[i]] %in% cDefault) ) {
					# cat("no recoding necessary!\n")
				} else {
					stop("calc.multiple (type==calcFullProblem):: recoding not possible!\n")
				}
			}
			
			## calculate all possible combinations within the lowest levels of dim-vars
			## if any combinations are missing (missing.codes), we have to set them to 0 later
			strID <- as.character(pasteStrVec(unlist(codes), length(codes)))
			exDims <- pasteStrVec(unlist(codes), length(codes))
			possDims <- sort(pasteStrVec(as.character(expand(lapply(dimObj, function(x) { get.dimVar(x, type='minimalCodesDefault') }), vector=TRUE)), length(dimObj)))
			missing.codes <- setdiff(possDims, exDims)
			
			## fill the table
			nrIndexvars <- length(ind.dimvars)
			fullDims <- lapply(dimObj, get.dimVar, type='dims')
			
			allCodes <- expand(lapply(dimObj, get.dimVar, type='codesDefault'), vector=FALSE)
			fullTabObj <- data.table(ID=1:length(allCodes[[1]]))
			for ( i in 1:length(allCodes)) {
				fullTabObj[,colnames(rawData)[ind.dimvars][i]:=allCodes[[i]]]
			}
			setkeyv(fullTabObj, colnames(rawData)[ind.dimvars])	
			fullTabObj[,ID:=NULL]
			
			## revert rawData codes to default codes
			for ( j in seq_along(ind.dimvars) ) {
				v <- calc.dimVar(object=dimObj[[j]], type="matchCodeDefault", input=rawData[,get(names(dimObj)[j])])
				rawData[,names(dimObj)[j]:=v]
			}
			setkeyv(rawData, colnames(rawData)[ind.dimvars])	
			
			## replace NAs in rawData by 0 (required for aggregation)	
			cols <- colnames(rawData)[(length(dimObj)+1):ncol(rawData)]
			ind.na <- list(); length(ind.na) <- length(cols)
			for ( i in 1:length(cols) ) {
				ind.na[[i]] <- which(is.na(rawData[,cols[i],with=FALSE]))
				if ( length(ind.na[[i]]) > 0 ) {
					rawData[ind.na[[i]], cols[i]:=0]
				}		
			}
			
			## merge minDat to fullDat
			fullTabObj <- merge(fullTabObj, rawData, all.x=TRUE)
			
			## set missing combinations of lowest levels to 0
			## problematic are all levels that should exist, but do not exist
			## they are filled with 0 so that we can aggregate
			dim.vars <- colnames(fullTabObj)[ind.dimvars]
			strID <- apply(fullTabObj[,dim.vars,with=FALSE],1,str_c, collapse="")
			
			if ( length(missing.codes) > 0 ) {
				index <- which(strID==missing.codes)
				for ( i in 1:length(cols) ) {
					fullTabObj[index, cols[i]:=0]
				}		
			}
			
			
			## fill up missing dimensions
			not.finished <- TRUE	
			while ( not.finished ) {
				cols <- (nrIndexvars+1):ncol(fullTabObj)
				col.names <- colnames(fullTabObj)[cols]
				for ( i in 1:nrIndexvars ) {
					setkeyv(fullTabObj, dim.vars[-i])	
					dat <- copy(fullTabObj) # we need to copy!
					
					cur.dim <- dimObj[[i]]@dims
					for ( j in length(cur.dim):1 ) {
						cur.levs <-  cur.dim[[j]]
						out <- dat[dat[[ind.dimvars[i]]] %in% cur.levs[-1],]
						out <- out[,lapply(.SD,sum), .SDcols=col.names, by=key(out)]
						
						row.ind <- which(fullTabObj[[ind.dimvars[i]]] == cur.levs[1])
						for ( z in 1:length(col.names) ) {
							v <- out[,col.names[z], with=FALSE]
							fullTabObj[row.ind, col.names[z]:=v]
						}
					}
				}
				if ( !is.na(fullTabObj[1,ind.freq,with=FALSE]) ) {
					not.finished <- FALSE
				}
			}
			
			nrV <- nrow(fullTabObj)	
			f <- fullTabObj[[ind.freq]]
			strID <- apply(fullTabObj[,dim.vars,with=FALSE],1,str_c, collapse="")
			w <- numVarsList <- NULL
			w.ind <- get.dataObj(x, type="weightVarInd")
			if ( !is.null(w.ind) ) {
				w <- fullTabObj[[w.ind]]
			}	
			n.ind <- get.dataObj(x, type="numVarInd")
			if ( !is.null(n.ind) ) {
				numVarsList <- list(); length(numVarsList) <- length(n.ind)
				for ( n in 1:length(n.ind) ) {
					numVarsList[[n]] <- fullTabObj[[n.ind[n]]]
				}
			}		
			
			## replace 0 in rawData by NA if they have been replaced earlier
			for ( i in 1:length(ind.na) ) {
				if ( length(ind.na[[i]]) > 0 ) {
					rawData[ind.na[[i]], cols[i]:=NA]
				}		
			}
			x <- set.dataObj(x, type="rawData", input=list(rawData))	
			
			problemInstance <- new("problemInstance",			
				strID=strID,
				Freq=f,
				w=w,
				numVars=numVarsList,
				lb=rep(0, nrV),
				ub=sapply(f, function(x) { max(2*x, 5)}),
				LPL=rep(1, nrV),
				UPL=rep(1, nrV),
				SPL=rep(0, nrV),
				sdcStatus=rep("s", nrV)
			)		
			partition <- calc.multiple(type='makePartitions', input=list(objectA=problemInstance, objectB=y))
			sdcProblem <- new("sdcProblem",
				dataObj=x,
				dimInfo=y,
				problemInstance=problemInstance,
				partition=partition,
				startI=1,
				startJ=1,
				indicesDealtWith=NULL,
				elapsedTime=(proc.time()-time.start)[3]
			)
			return(sdcProblem)	
		}
	}
)
