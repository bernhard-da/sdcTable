#' @aliases calc.multiple,character,list-method
#' @rdname calc.multiple-method
setMethod(f='calc.multiple', signature=c('character', 'list'),
	definition=function(type, input) {
		if (!type %in% c('makePartitions', 'genMatMFull', 
				'makeAttackerProblem', 'calcFullProblem') ) {
			stop("calc.multiple:: argument 'type' is not valid!\n")
		}	

		if ( type == 'makePartitions' ) {
			pI <- input$objectA
			dimInfoObj <- input$objectB
			dimInfo <- get.dimInfo(dimInfoObj, type='dimInfo')
			strIDs <- get.problemInstance(pI, type='strID')
			
			### create classes and groups
			tmpDat <- expand.grid(lapply(1:length(dimInfo), function(x) { 1:get.dimVar(dimInfo[[x]], type='nrLevels') } ))
			groups <- apply(tmpDat, 1, function(x) { paste(x, collapse="-")})
			classes <- apply(tmpDat, 1, sum)
			sortOrder <- order(classes)
			classes <- classes[sortOrder]
			classesUnique <- unique(classes)
			groups <- groups[sortOrder]
			splitGroups <- split(groups, classes)
			
			### create tables for all classes and groups 
			final <- list()
			final$groups <- as.list(groups)
			final$indices <- list()
			for ( i in 1:length(groups) ) {
				final$indices[[i]] <- list()
				levs <- as.integer(unlist(sapply(groups[[i]], strsplit, "-")))	
				
				###
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
						} # end for(k)-loop				
					} # end else
				} # end for(z)-loop
				final$indices[[i]] <- list()
				final$strIDs[[i]] <- list()
				combs <- expand.grid(lapply(1:length(res), function(x) { 1:length(res[[x]]) } ))
				for ( m in 1:nrow(combs) ) {
					final$strIDs[[i]][[m]] <- pasteStrVec(expand(lapply(1:ncol(combs), function(x) { res[[x]][[combs[m,x]]] })), ncol(combs))
					final$indices[[i]][[m]] <- which(strIDs %in% final$strIDs[[i]][[m]])
				}			
			} # end for(i)-loop
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
			f <- pasteStrVec(expand(lapply(levelObj, function(x) { get.dimVar(x, type='codesDefault')})), length(levelObj))	
			order <- match(strID,f)
			f <- f[order]	
			
			constraintM <- init.simpleTriplet(type='simpleTriplet', input=list(mat=matrix(0, nrow=0, ncol=nrCells)))
			for ( i in 1:nrVars ) {
				lO <- levelObj[[i]]
				keepList <- lapply(get.dimInfo(y, type='strInfo')[-i], function(k) { seq(k[1], k[2]) } )
				keepList2 <- lapply(get.dimInfo(y, type='strInfo')[i], function(k) { seq(k[1], k[2]) } )
				#f2 <- sapply(strID, function(x) { mySplitIndicesList(x, keepList2) } )
				f2 <- mySplitIndicesList(strID, keepList2)
				
				if ( nrVars > 1 ) { 
					#f1 <- as.vector(sapply(strID, function(x) { mySplitIndicesList(x, keepList) } ))
					f1 <- as.vector(mySplitIndicesList(strID, keepList))
				} else {
					f1 <- f2
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
								if ( sum(freqs[spl[[z]]]*ind) != 0 ) {
									stop("something went wrong!\n")
								}
								constraintM <- calc.simpleTriplet(constraintM, type='addRow', input=list(index=spl[[z]], values=ind))
							}
						}				
					} else {
						splitInd <- which(f2 %in% get.dimVar(lO, type='dims')[[which(dimInd==TRUE)]])
						### only 1 dimension?
						if ( nrVars > 1 ) {
							spl <- split(splitInd, f1[splitInd])	
						} else {
							spl <- split(splitInd, rep(1, length(splitInd)))
						}
						
						for ( z in 1:length(spl) ) {
							ind <- rep(1,length(spl[[z]])) 
							ind[which.max(freqs[spl[[z]]])] <- -1
							if ( sum(freqs[spl[[z]]]*ind) != 0 ) {
								stop("something went wrong! (z=",z,"und names(spl)[z]=",names(spl)[z],")\n")
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
			
			### calculating (logical) constraints for the master problem ###
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
			
			dimList <- numList <- weightList <- freqList <- NULL
			dimList <- rawData[get.dataObj(x, type='dimVarInd')]
			freqList <- rawData[get.dataObj(x, type='freqVarInd')]
			
			if ( !is.null(get.dataObj(x, type='weightVarInd')) ) {
				weightList <- rawData[get.dataObj(x, type='weightVarInd')]	
			}
			if ( !is.null(get.dataObj(x, type='numVarInd')) ) {
				numList <- rawData[get.dataObj(x, type='numVarInd')]
			}	
			### match dimVars in rawData to valid codes ###
			dimVarsData <- get.dataObj(x, type='dimVarInd')
			dimVars <- get.dimInfo(y, type='posIndex')
			
			# not all possible dimensional variables are used
			# we need to aggregate
			if ( length(dimVars) < length(dimVarsData) | get.dataObj(x, type='isMicroData') == TRUE ) {
				if ( length(dimVars) < length(dimVarsData) ) {
					dimList <- dimList[-setdiff(dimVarsData, dimVars)]	
				}
				xx <- as.list(aggregate(freqList, dimList, sum))
				if ( !is.null(numList) ) {
					yy <- as.list(aggregate(numList, dimList, sum))
				}				
				if ( !is.null(weightList) ) {
					zz <- as.list(aggregate(weightList, dimList, sum))	
				}
				dimList <- xx[1:(length(xx)-1)]
				freqList <- xx[length(xx)]
				
				if ( !is.null(numList) ) {
					numList <- yy[(length(dimVars)+1):(length(yy))]
				}
				if ( !is.null(weightList) ) {
					weightList <-zz[(length(dimVars)+1):(length(zz))]
				}			
			}
			
			codes <- dimList
			for ( i in 1:length(codes) ) {
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
					codes[[i]] <- sapply(codes[[i]], function(x) { calc.dimVar(dimObj[[i]], type='matchCodeDefault', input=x) } )
				} else if ( all(codes[[i]] %in% cDefault) ) {
					# cat("no recoding necessary!\n")
				} else {
					stop("calc.multiple (type==calcFullProblem):: recoding not possible!\n")
				}
			}
			strID <- as.character(pasteStrVec(unlist(codes), length(dimVars)))
			
			tableVars <- get.dimInfo(y, type='varName')
			minTabObj <- list()
			numVars <- get.dataObj(x, type='numVarNames')
			dimInfo <- get.dimInfo(y, type='dimInfo') # list with objects of class 'dimVar'
			
			exDims <- pasteStrVec(unlist(codes), length(codes))
			possDims <- sort(pasteStrVec(as.character(expand(lapply(dimInfo, function(x) { get.dimVar(x, type='minimalCodesDefault') }), vector=TRUE)), length(dimInfo)))
			
			# problematic are all levels that should exist, but do not exist
			# they are filled with NA/0
			probl <- setdiff(possDims, exDims)
			if ( length(probl) > 0 ) {
				strID <- c(strID, probl)
				freqList[[1]] <- c(freqList[[1]], rep(0, length(probl)))
				if ( !is.null(get.dataObj(x, type='weightVarInd')) ) {
					weightList[[1]] <- c(weightList[[1]], rep(0, length(probl)))
				}
				if ( !is.null(get.dataObj(x, type='numVarInd')) ) {
					lapply(1:length(numList), function(x) { numList[[x]] <<- c(numList[[x]], rep(0, length(probl))) } )
				}
			}	
			
			### fill the table
			nrIndexvars <- length(dimInfo)
			fullDims <- lapply(dimInfo, get.dimVar, type='dims')
			
			# complete
			allDims <- pasteStrVec(expand(lapply(dimInfo, get.dimVar, type='codesDefault')), nrIndexvars)
			
			# the subtotals that need to be calculated
			subTotals <- setdiff(allDims, strID)
			
			fullTabObj 			<- minTabObj
			fullTabObj$strID 	<- c(strID, subTotals)
			fullTabObj$Freq 	<- c(freqList[[1]], rep(NA, length(subTotals)))
			
			#fullTabObj$w <- NULL
			if ( !is.null(get.dataObj(x, type='weightVarInd')) ) {
				fullTabObj$w <- c(unlist(weightList), rep(NA, length(subTotals)))
				attributes(fullTabObj) <- NULL
			}
			
			numVarsList <- NULL
			if ( !is.null(get.dataObj(x, type='numVarInd')) ) {
				for ( k in 1:length(numList) ) {
					numVarsList[[k]] <- c(as.vector(numList[[k]]), rep(NA, length(subTotals)))
					names(numVarsList)[k] <- numVars[k]
				}
			}
			
			info <- get.dimInfo(y, type='strInfo')
			minI <- 1
			maxI <- max(unlist(info))
			
			splitFactors <- keepIndices <- list()
			codeInfo <- numCodes <- numSplitFactors <- list()
			for( i in 1:nrIndexvars ) {	
				codeInfo[[i]] <- list()
				keepIndices[[i]] <- seq(info[[i]][1], info[[i]][2]) 
				splitFactors[[i]] <- mySplit(fullTabObj$strID, keepIndices[[i]])
				
				# prepare all possible codes only once!
				splitF <- unique(splitFactors[[i]])
				for ( j in 1:length(splitF) ) {
					codeInfo[[i]][[j]] <- list()
					codeInfo[[i]][[j]]$upper <- splitF[j]
					codeInfo[[i]][[j]]$lower <- calc.dimVar(dimObj[[i]], type='requiredMinimalCodes', input=splitF[j])
				}
				numCodes[[i]] <- as.numeric(sapply(codeInfo[[i]], function(x) { x$upper}))
				numSplitFactors[[i]] <- as.numeric(splitFactors[[i]])
			}
			facsRawData <- pasteStrVec(unlist(rawData[get.dataObj(x, type='dimVarInd')]), length(get.dataObj(x, type='dimVarInd')), coll="-")
			
			runIndex <- which(is.na(fullTabObj$Freq))
			freqVarInd <- get.dataObj(x, type='freqVarInd')
			numVarInd <- get.dataObj(x, type='numVarInd')
			weightVarInd <- get.dataObj(x, type='weightVarInd')
			
			for ( z in runIndex ) {
				codes <- lapply(1:nrIndexvars, function(x) {
					codeInfo[[x]][[which(numCodes[[x]] == numSplitFactors[[x]][z])]]$lower
				})
				
				facsCodes <- unique(pasteStrVec(expand(codes, vector=TRUE), length(codes), coll="-"))
				
				index <- which(facsRawData %in% facsCodes)
				fullTabObj[["Freq"]][z] <- sum(rawData[[freqVarInd]][index])
				
				if ( !is.null(weightVarInd) ) {
					fullTabObj[["w"]][z] <- sum(rawData[[weightVarInd]][index]) 	
				}
				
				if ( !is.null(numVarInd) ) {
					res <- lapply(numVarInd, function(u) { sum(rawData[[u]][index])  } )
					for ( k in 1:length(numVarsList) ) {
						numVarsList[[k]][z] <- res[[k]]
					} 
				}	
			}
			nrV <- length(fullTabObj$strID)
			
			problemInstance <- new("problemInstance",			
				strID=fullTabObj$strID,
				Freq=fullTabObj$Freq,
				w=fullTabObj$w,
				numVars=numVarsList,
				lb=rep(0, nrV),
				ub=sapply(fullTabObj$Freq, function(x) { max(2*x, 5)}),
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
