#' @aliases get.sdcProblem,sdcProblem,character-method
#' @rdname get.sdcProblem-method
setMethod(f='get.sdcProblem', signature=c('sdcProblem','character'),
  definition=function(object, type) { 
    if ( !type %in% c('dataObj', 'problemInstance', 'partition', 'elapsedTime', 'dimInfo', 'indicesDealtWith',
        'startI', 'startJ', 'innerAndMarginalCellInfo') ) {
      stop("get.sdcProblem:: argument 'type' is not valid!\n")
    }
    if ( type == 'dataObj' ) {
      return(object@dataObj)
    }
    if ( type == 'problemInstance' ) {
      return(object@problemInstance)
    }   
    if ( type == 'partition' ) {
      return(object@partition)
    } 
    if ( type == 'elapsedTime' ) {
      return(object@elapsedTime)
    } 
    if ( type == 'dimInfo' ) {
      return(object@dimInfo)
    }     
    if ( type == 'indicesDealtWith' ) {
      return(object@indicesDealtWith)
    }   
    if ( type == 'startI' ) {
      return(object@startI)
    } 
    if ( type == 'startJ' ) {
      return(object@startJ)
    } 
    
    if ( type == 'innerAndMarginalCellInfo' ) {
      pI <- get.sdcProblem(object, type='problemInstance')
      strIDs <- get.problemInstance(pI, type='strID')
      strInfo <- get.dimInfo(get.sdcProblem(object, type='dimInfo'), type='strInfo')
      
      out <- lapply(1:length(strInfo), function(x) { sort(unique(mySplit(strIDs, strInfo[[x]][1]:strInfo[[x]][2])))[-1] } )
      
      # deal with 'tot' levels
      ind <- which(sapply(out, length) ==0)
      out[ind] <- lapply(ind, function(x) { "0" } )
      
      innerCells <- apply(matrix(unlist(expand(out)), ncol=length(out), byrow=FALSE),1,paste, collapse="")
      totCells <- setdiff(strIDs, innerCells)
      indexTotCells <- match(totCells, strIDs)
      indexInnerCells <- match(innerCells, strIDs)
      return(list(innerCells=innerCells, totCells=totCells, indexInnerCells=indexInnerCells, indexTotCells=indexTotCells))
    }
  }
)

#' @aliases set.sdcProblem,sdcProblem,character,list-method
#' @rdname set.sdcProblem-method
setMethod(f='set.sdcProblem', signature=c('sdcProblem', 'character', 'list'),
  definition=function(object, type, input) { 
    if ( !type %in% c('problemInstance', 'partition', 'rule.freq', 'rule.nk', 
        'rule.p', 'rule.pk', 'startI', 'startJ', 'indicesDealtWith', 'elapsedTime') ) {
      stop("set.sdcProblem:: check argument 'type'!\n")
    }
    
    if ( type == 'problemInstance' ) {
      object@problemInstance <- input[[1]]  
    }
    
    if ( type == 'partition' ) {
      object@partition <- input[[1]]  
    }
    
    if ( type == 'startI' ) {
      object@startI <- input[[1]]
    }
    if ( type == 'startJ' ) {
      object@startJ <- input[[1]]
    } 
    if ( type == 'indicesDealtWith' ) {
      object@indicesDealtWith <- input[[1]]
    } 
    if ( type == 'elapsedTime' ) {
      object@elapsedTime <- input[[1]]
    }       
    
    validObject(object)
    return(object) 
  }
)

#' @aliases calc.sdcProblem,sdcProblem,character,list-method
#' @rdname calc.sdcProblem-method
setMethod(f='calc.sdcProblem', signature=c('sdcProblem', 'character', 'list'),
  definition=function(object, type, input) { 
    if ( !type %in% c('rule.freq', 'rule.nk', 'rule.p', 'rule.pq', 'heuristicSolution',
      'cutAndBranch', 'anonWorker', 'ghmiter', 'preprocess', 'cellID', 
      'finalize', 'ghmiter.diagObj', 'ghmiter.calcInformation', 
      'ghmiter.suppressQuader', 'ghmiter.selectQuader', 
      'ghmiter.suppressAdditionalQuader', 'contributingIndices', 
      'reduceProblem', 'genStructuralCuts') ) {
      stop("calc.sdcProblem:: check argument 'type'!\n")
    }
  
    # frequency-rule
    if ( type == 'rule.freq' ) {
      pI <- get.sdcProblem(object, type='problemInstance')
      
      if ( input$allowZeros == TRUE ) {
        suppInd <- which(get.problemInstance(pI, type='freq') <= input$maxN)  
      } else {
        f <- get.problemInstance(pI, type='freq')
        suppInd <- which(f > 0 & f <= input$maxN)         
        zeroInd <- which(get.problemInstance(pI, type='freq') == 0 )
        if ( length(zeroInd) > 0 ) {
          pI <- set.problemInstance(pI, type='sdcStatus', input=list(index=zeroInd, values=rep("z", length(zeroInd))))  
        }
      }
      
      if ( length(suppInd) > 0 ) {
        pI <- set.problemInstance(pI, type='sdcStatus', input=list(index=suppInd, values=rep("u", length(suppInd))))  
      }
      
      object <- set.sdcProblem(object, type='problemInstance', input=list(pI))
      return(object)
    }
    
    # nk-dominance rule
    if ( type == 'rule.nk' ) {
      nkRule <- function(celltot, sumNcont, k) {
        # if TRUE, cell needs to be suppressed
        (sumNcont) > (k/100*celltot) 
      } 
      if ( !get.dataObj(get.sdcProblem(object, type='dataObj'), type='isMicroData') ) {
        stop("nk-dominance rule can only be applied if micro-data are available!\n")
      }     
      pI <- get.sdcProblem(object, type='problemInstance')
      dataObj <- get.sdcProblem(object, type='dataObj')
      strIDs <- get.problemInstance(pI, type='strID')
      numVarInds <- get.dataObj(dataObj, type='numVarInd')
      
      numVal <- get.dataObj(dataObj, type='rawData')[[numVarInds[input$numVarInd]]]
      if ( any(numVal < 0 ) ) {
        stop("dominance rules can only be applied to numeric variables with only positive values!\n")
      }     
   
      # calculate contributing indices
      indices <- lapply(1:get.problemInstance(pI, type='nrVars'), function(x) { calc.sdcProblem(object, type='contributingIndices', input=list(strIDs[x])) } )
      
      minContributingUnits <- min(setdiff(unique(sapply(indices, length)), 0))    
  
      if ( input$n < 1 | input$n > minContributingUnits ) {
        stop("set.sdcProblem:: parameter 'n' must be >= 1 and <",minContributingUnits,"!\n")
      }   
      
      # values of contributing units
      valueList <- lapply(1:get.problemInstance(pI, type='nrVars'), function(x) { sum(rev(tail(sort(numVal[indices[[x]]]), input$n))) } ) 
      cellTotals <- get.problemInstance(pI, type='numVars')[[input$numVarInd]]    
      # suppStatus: TRUE:unsafe, FALSE: safe
      nkState <- sapply(1:get.problemInstance(pI, type='nrVars'), function(x) {  nkRule(cellTotals[x], valueList[[x]], input$k) } ) 
      
      addSupps <- which(sapply(indices, length) %in% 1:input$n)
      if ( length(addSupps) > 0 ) {
        nkState[addSupps] <- TRUE
      }
      
      suppIndex <- which(nkState==TRUE)
      if ( length(suppIndex) > 0 ) {
        pI <- set.problemInstance(pI, type='sdcStatus', input=list(index=suppIndex, values=rep("u", length(suppIndex))))  
      }
      
      if ( input$allowZeros == FALSE ) {
        indZero <- which(get.problemInstance(pI, type='freq')==0)
        if ( length(indZero) > 0 ) {
          pI <- set.problemInstance(pI, type='sdcStatus', input=list(index=indZero, values=rep("z", length(indZero))))
        }     
      }         
      object <- set.sdcProblem(object, type='problemInstance', input=list(pI))
      return(object)
    }   
    
    # p-percent rule
    if ( type %in% c('rule.p', 'rule.pq') ) {
      pPercRule <- function(celltot, cont1, cont2, p) {
        # if TRUE, cell needs to be suppressed
        (celltot - cont1 - cont2) < (p/100*cont1) 
      }     
      pqRule <- function(celltot, cont1, cont2, p, q) {
        # if TRUE, cell needs to be suppressed        
        (celltot - cont1 - cont2) < (p/q)*cont1 
      }         
      if ( !get.dataObj(get.sdcProblem(object, type='dataObj'), type='isMicroData') ) {
        stop("p-percent rule can only be applied if micro-data are available!\n")
      }
      
      pI <- get.sdcProblem(object, type='problemInstance')
      dataObj <- get.sdcProblem(object, type='dataObj')
      numVarInds <- get.dataObj(dataObj, type='numVarInd')
      strIDs <- get.problemInstance(pI, type='strID')
      
      numVal <- get.dataObj(dataObj, type='rawData')[[numVarInds[input$numVarInd]]]
      if ( any(numVal < 0 ) ) {
        stop("dominance rules can only be applied to numeric variables with only positive values!\n")
      }     

      # calculate contributing indices
      indices <- lapply(1:get.problemInstance(pI, type='nrVars'), function(x) { calc.sdcProblem(object, type='contributingIndices', input=list(strIDs[x])) } )
      
      # values of contributing units
      valueList <- lapply(1:get.problemInstance(pI, type='nrVars'), function(x) { rev(tail(sort(numVal[indices[[x]]]),2)) } ) 
      cellTotals <- get.problemInstance(pI, type='numVars')[[input$numVarInd]]    
      
      # suppStatus: TRUE:unsafe, FALSE: safe
      if (type  == 'rule.p' ) {
        pState <- sapply(1:get.problemInstance(pI, type='nrVars'), function(x) { pPercRule(cellTotals[x], valueList[[x]][1], valueList[[x]][2], input$p) } ) 
      }
      if (type  == 'rule.pq' ) {
        pState <- sapply(1:get.problemInstance(pI, type='nrVars'), function(x) { pqRule(cellTotals[x], valueList[[x]][1], valueList[[x]][2], input$pq[1], input$pq[2]) } ) 
      }
      
      suppIndex <- which(pState==TRUE)
      if ( length(suppIndex) > 0 ) {
        pI <- set.problemInstance(pI, type='sdcStatus', input=list(index=suppIndex, values=rep("u", length(suppIndex))))  
      }
      
      if ( input$allowZeros == FALSE ) {
        indZero <- which(get.problemInstance(pI, type='freq')==0)
        if ( length(indZero) > 0 ) {
          pI <- set.problemInstance(pI, type='sdcStatus', input=list(index=indZero, values=rep("z", length(indZero))))
        }     
      }     
      object <- set.sdcProblem(object, type='problemInstance', input=list(pI))
      return(object)
    }   
    
    if ( type == 'heuristicSolution' ) {
      aProb <- input[[1]]
      validCuts <- input[[2]]
      solver <- input[[3]]$solver
      verbose <- input[[3]]$verbose
      
      ### create incremental attacker problem
      pI <- get.sdcProblem(object, type='problemInstance')
      dimInfoObj <- get.sdcProblem(object, type='dimInfo')
      primSupps <- get.problemInstance(pI, type='primSupps')
      secondSupps <- get.problemInstance(pI, type='secondSupps')
      forcedCells <- get.problemInstance(pI, type='forcedCells')
      nrVars <- get.problemInstance(pI, type='nrVars')
      weights <- ci <- get.problemInstance(pI, type='weight')
      ci[primSupps] <- 0
      ci[secondSupps] <- 0
      
      lb <- get.problemInstance(pI, type='lb')
      ub <- get.problemInstance(pI, type='ub')
      
      # required later in the cleanup-procedure
      LB <- LBdefault <- weights - lb
      UB <- UBdefault <- ub - weights
      
      m1 <- calc.multiple(type='genMatMFull', input=list(objectA=pI, objectB=dimInfoObj))
      m2 <- m1
      m2@v <- -1* m2@v
      AInc <- calc.simpleTriplet(object=m1, type='bind', input=list(m2, bindRow=FALSE))
      
      nrConstraints <- nrow(AInc)
      objective <- rep(ci, 2)
      
      direction <- rep("==", get.simpleTriplet(AInc, type='nrRows', input=list()))
      rhs <- rep(0, get.simpleTriplet(AInc, type='nrRows', input=list()))
      
      types <- rep("C", get.simpleTriplet(AInc, type='nrCols', input=list()))
      boundsLower <- list(ind=1:get.simpleTriplet(AInc, type='nrCols', input=list()), val=rep(0,get.simpleTriplet(AInc, type='nrCols', input=list())))
      boundsUpper <- list(ind=1:get.simpleTriplet(AInc, type='nrCols', input=list()), val=c(UB, LB))
      
      aProbInc <- new("linProb",
        objective=objective,
        constraints=AInc,
        direction=direction,
        rhs=rhs,
        boundsLower=boundsLower,
        boundsUpper=boundsUpper,
        types=types)
      
      # make sure that cells that must be published
      # are not part of the heuristic solution
      if ( length(forcedCells) > 0 ) {
        for ( u in 1:length(forcedCells) ) {
          con <- rep(0, get.simpleTriplet(AInc, type='nrCols', input=list()))
          con[c(forcedCells[u], nrVars+forcedCells[u])] <- c(1,-1)
          aCon <- init.cutList(type='singleCut', input=list(vals=con, dir="==", rhs=0))
          aProbInc <- set.linProb(aProbInc, type='addCompleteConstraint', input=list(aCon))
        }
      }
      
      x <- rep(0, get.simpleTriplet(AInc, type='nrCols', input=list()))
      UPL <- get.problemInstance(pI, type='UPL')
      LPL <- get.problemInstance(pI, type='LPL')
      SPL <- get.problemInstance(pI, type='SPL')
      
      SUP <- primSupps
      
      for ( i in seq_along(primSupps) ) {
        cellInd <- primSupps[i]
        if ( verbose ) {
          cat("finding additional cells to protect primSupp",i,"|",length(primSupps),"...\n")
        }
        con1 <- con2 <- x
        con1[cellInd] <- 1
        con2[nrVars+cellInd] <- 1
        con3 <- con1 - con2 # page 1018: fichetti and salazar!! (- is correct!)
        if ( UPL[cellInd] > 0 ) {
          # update and solve: y_ik_minus <- 0 and y_ik_plus <- UPL_ik
          aCon <- init.cutList(type='multipleCuts', input=list(mat=init.simpleTriplet(type='simpleTriplet', input=list(mat=rbind(con1, con2))), dir=rep("==", 2), rhs=c(UPL[cellInd], 0)))
          
          prob <- set.linProb(aProbInc, type='addCompleteConstraint', input=list(aCon))
          sol <- calc.linProb(prob, type='solveProblem', input=list(solver))$solution
          v <- sol[1:nrVars]+sol[(nrVars+1):length(sol)]
          v[which(is.zero(v))] <- 0
          addIndex <- which ( v > 0 )
          if ( length(addIndex) > 0 ) {
            SUP <- unique(c(SUP, addIndex))#
            ci[SUP] <- 0
            aProbInc <- set.linProb(aProbInc, type='objective', input=list(rep(ci, 2)))
            LB <- ci - lb
            UB <- ub - ci
          }
        }
        if ( LPL[cellInd] > 0 ) {
          # update and solve: y_ik_minus <- LPL_ik and y_ik_plus <- 0
          aCon <- init.cutList(type='multipleCuts', input=list(mat=init.simpleTriplet(type='simpleTriplet', input=list(mat=rbind(con1, con2))), dir=rep("==", 2), rhs=c(0, LPL[cellInd])))
          prob <- set.linProb(aProbInc, type='addCompleteConstraint', input=list(aCon))
          sol <- calc.linProb(prob, type='solveProblem', input=list(solver))$solution
          v <- sol[1:nrVars]+sol[(nrVars+1):length(sol)]
          v[which(is.zero(v))] <- 0
          addIndex <- which ( v > 0 )
          if ( length(addIndex) > 0 ) {
            SUP <- unique(c(SUP, addIndex))
            ci[SUP] <- 0
            aProbInc <- set.linProb(aProbInc, type='objective', input=list(rep(ci, 2)))
            LB <- ci - lb
            UB <- ub - ci
          }
        }
        if ( SPL[cellInd] > 0 ) {
          # update and solve: y_ik_plus + y_ik_minus <- SPL_ik
          aCon <- init.cutList(type='singleCut', input=list(vals=con3, dir="==", rhs=SPL[cellInd]))
          prob <- set.linProb(aProbInc, type='addCompleteConstraint', input=list(aCon))
          sol <- calc.linProb(prob, type='solveProblem', input=list(solver))$solution
          v <- sol[1:nrVars]+sol[(nrVars+1):length(sol)]
          v[which(is.zero(v))] <- 0
          addIndex <- which ( v > 0 )
          if ( length(addIndex) > 0 ) {
            SUP <- unique(c(SUP, addIndex))
            ci[SUP] <- 0
            aProbInc <- set.linProb(aProbInc, type='objective', input=list(rep(ci, 2)))
            LB <- ci - lb
            UB <- ub - ci
          }
        }
      }
      if ( verbose ) {
        cat(length(SUP) - length(primSupps),"additional cells have been suppressed in the heuristic solution!\n")
      }
    
      ### cleanup: remove redundant suppressions....
      # FIXME: use constraint pool to search for violations in the constraint pool
      # aProb has already been calculated and is an input parameter of this method!
      # aProb <- calc.multiple(type='makeAttackerProblem', input=list(objectA=pI, objectB=dimInfoObj))
      nrConstraints <- length(get.linProb(aProb, type='objective')) - 2*length(weights)
      addedSupps <- setdiff(SUP, primSupps)
      orderAddedSupps <- order(weights[addedSupps], decreasing=TRUE)
      xi <- rep(0, length(UPL))
      xi[SUP] <- 1 # we need to check xi
      
      counter <- 0
      for ( i in orderAddedSupps ) {
        counter <- counter + 1
        if ( verbose ) {
          cat("checking if removing cell",counter,"|",length(addedSupps),"still yields a valid suppression pattern...\n")
        }
          
        cellInd <- addedSupps[i]
        limits <- c(LPL[cellInd], UPL[cellInd], SPL[cellInd])
        xiWorking <- xi
        xiWorking[cellInd] <- 0 # we need to check if xi without xi[i] is a valid pattern
        UBWorking <- UB
        LBWorking <- LB
        UBWorking[cellInd] <- UBdefault[cellInd]
        LBWorking[cellInd] <- LBdefault[cellInd]
        
        ######################
        # check if any validCuts are not valid with xiWorking!
        # validCuts needs to be supplemented in function call
        ### get a constraint from validCuts
        if ( get.cutList(validCuts, type='nrConstraints') > 0 ) {
          conMat <- get.cutList(validCuts, type='constraints')
          result <- rep(NA, get.simpleTriplet(conMat, type='nrRows', input=list()))
          for ( z in 1:get.simpleTriplet(conMat, type='nrRows', input=list()) ) {
            expr <- paste(sum(xiWorking[get.simpleTriplet(get.simpleTriplet(conMat, type='getRow', input=list(z)), type='colInd', input=list())]),  get.cutList(validCuts, type='direction')[z], get.cutList(validCuts, type='rhs')[z])
            result[z] <- eval(parse(text=expr))
          }
        } else {
          result <- TRUE
        }
        
        if ( any(result==FALSE) ) {
          #cat("additionally suppressed cell cannot be removed (violated constraint in the pool found)!\n")
        } else {
          # no constraint was violated, we have to solve the incremental attacker problems
          # we need to solve the attackers problem for each sensitive cell twice
          if ( limits[3] != 0 ) {
            # solveAttackerProblem (upper bound)
            rhs <- rep(0, length(get.linProb(aProb, type='rhs')))
            rhs[cellInd] <- 1
            aProb <- set.linProb(aProb, type='rhs', input=list(rhs))
            aProb <- set.linProb(aProb, type='objective', input=list(c(weights + UBWorking*xiWorking, -(weights-xiWorking*LBWorking), rep(0, nrConstraints))))
            up <- calc.linProb(aProb, type='solveProblem', input=list(solver))
            
            # solveAttackerProblem (lower bound)
            aProb <- set.linProb(aProb, type='rhs', input=list(-1*rhs))
            down <- calc.linProb(aProb, type='solveProblem', input=list(solver))
            
            calcDown <- -down$optimum
            calcUp <- up$optimum
          } else {
            # solve attackers problem (minimize)
            if ( limits[1] != 0 ) {
              rhs <- rep(0, length(get.linProb(aProb, type='rhs')))
              rhs[cellInd] <- -1
              aProb <- set.linProb(aProb, type='rhs', input=list(rhs))
              aProb <- set.linProb(aProb, type='objective', input=list(c(weights + UBWorking*xiWorking, -(weights-xiWorking*LBWorking), rep(0, nrConstraints))))
              down <- calc.linProb(aProb, type='solveProblem', input=list(solver))
              calcDown <- -down$optimum
            }
            # solve attackers problem (maximize)
            if ( limits[2] != 0 ) {
              rhs <- rep(0, length(get.linProb(aProb, type='rhs')))
              rhs[cellInd] <- 1
              aProb <- set.linProb(aProb, type='rhs', input=list(rhs))
              aProb <- set.linProb(aProb, type='objective', input=list(c(weights + UBWorking*xiWorking, -(weights-xiWorking*LBWorking), rep(0, nrConstraints))))
              up <- calc.linProb(aProb, type='solveProblem', input=list(solver))
              calcUp <- up$optimum
            }
          }
          
          # check for feasibility
          valid <- TRUE
          if ( limits[3] > 0 & calcUp - calcDown < SPL[i] ) {
            valid <- FALSE
          } else {
            if ( limits[1] > 0 & calcDown > weights[i] - LPL[i] ) {
              valid <- FALSE
            }
            if ( limits[2] > 0 & calcUp < weights[i] + UPL[i] ) {
              valid <- FALSE
            }
          }
          if ( valid ) {
            xi[cellInd] <- 0
            SUP <- setdiff(SUP, cellInd)
            if ( verbose ) {
              cat("redundant suppression found! --> removing cell!\n")
            }
          }
          # else: additionally suppressed cell cannot be removed!
        }
      }
      return(xi)      
    }
    
    if ( type == 'cutAndBranch' ) {
      time.start <- proc.time()
      
      timeLimit <- input$timeLimit
      fixVariables <- input$fixVariables
      maxVars <- input$maxVars
      fastSolution <- input$fastSolution
      approxPerc <- input$approxPerc
      verbose <- input$verbose
      solver <- input$solver
      
      problemInstance <- get.sdcProblem(object, type='problemInstance')
      dimInfo <- get.sdcProblem(object, type='dimInfo')
      nrVars <- get.problemInstance(problemInstance, type='nrVars')
      freqs <- get.problemInstance(problemInstance, type='freq')
      primSupps <- get.problemInstance(problemInstance, type='primSupps')
      publishVars <- which(get.problemInstance(problemInstance, type='sdcStatus')=="z")
      noBranchVars <- unique(c(primSupps, publishVars))
      
      # Nothing to protect here
      if ( !get.problemInstance(problemInstance, type='hasPrimSupps') ) {
        return(object)
      }
      
      # returning heuristic solution
      # only if problem size is too large
      if ( is.null(maxVars) ) {
        maxVars <- nrVars + 1
      }
      if ( fastSolution ) {
        maxVars <- 0
      }
      
      approx <- ifelse(is.null(approxPerc), FALSE, TRUE)
      
      if ( nrVars >= maxVars ) {
        res <- calc.multiple(type='makeAttackerProblem', input=list(objectA=problemInstance, objectB=dimInfo))
        aProb <- res$aProb
        validCuts <- res$newCutsMaster
        
        heuristicSolution <- calc.sdcProblem(object, type='heuristicSolution', input=list(aProb, validCuts, input))
        
        secondSupps <- setdiff(which(heuristicSolution==1), primSupps)
        if (verbose) {
          cat("Result: we are returning a possibly non-optimal solution with",length(secondSupps),"secondary suppressions because of parameter 'fastSolution' or 'maxVars'!\n")
        }
        if ( length(secondSupps) > 0 ) {
          problemInstance <- set.problemInstance(problemInstance, type='sdcStatus', input=list(index=secondSupps, values=rep("x", length(secondSupps))))
        }
        out <- new("sdcProblem",
          dataObj=get.sdcProblem(object, type='dataObj'),
          dimInfo=dimInfo,
          problemInstance=problemInstance,
          partition=get.sdcProblem(object, type='partition'),
          startI=get.sdcProblem(object, type='startI'),
          startJ=get.sdcProblem(object, type='startJ'),
          indicesDealtWith=get.sdcProblem(object, type='indicesDealtWith'),
          elapsedTime=get.sdcProblem(object, type='elapsedTime')+(proc.time()-time.start)[3]
        )
        return(out)
      }
      
      if ( verbose ) {
        cat("running pre-process procedure...\n")
      }
      
      resultPreProcess <- calc.sdcProblem(object, type='preprocess', input=input)
      
      object <- resultPreProcess$sdcProblem
      validCuts <- resultPreProcess$validCuts
      aProb <- resultPreProcess$aProb
      
      # no valid cuts have been generated in preprocessing!
      if ( get.simpleTriplet(get.cutList(validCuts, type='constraints'), type='nrRows', input=list()) == 0 ) {
        return(object)
      }
      
      if ( verbose ) {
        cat("calculating a heuristic solution...\n")
      }
        
      heuristicSolution <- calc.sdcProblem(object, type='heuristicSolution', input=list(aProb, validCuts, input))
      
      ### all solutions found and current best solution
      solutions <- list(); bestSolution <- heuristicSolution
      solutions[[1]] <- heuristicSolution
      
      startTime <- Sys.time()
      timeStop <- FALSE
      
      ### cuts due to hierarchical structure
      if ( verbose ) {
        cat("calculating structural cuts...\n")
      }
        
      structureCuts <- calc.sdcProblem(object, type='genStructuralCuts', input=list())
      # does heuristicSolution violates any cuts from structure-cuts??
      #calc.cutList(structureCuts, type='checkViolation', input=list(heuristicSolution, get.problemInstance(problemInstance, type='weight')))
      validCuts <- calc.cutList(validCuts, type='bindTogether', input=list(structureCuts))
      #######
      
      ### create master problem and add constraints derived in pre-processing
      mProb <- calc.problemInstance(problemInstance, type='makeMasterProblem', input=list())
      mProb <- set.linProb(mProb, type='addCompleteConstraint', input=list(validCuts))
      if ( verbose ) {
        cat("solving the original master problem (no additional constraints)...\n")
      }
      masterSolution <- calc.linProb(mProb, type='solveProblem', input=list(solver))
      masterObj <- masterSolution$optimum
      xi <- masterSolution$solution
      xi[is.zero(xi)] <- 0
      xi[is.one(xi)] <- 1
      
      ### initialize bounds
      currentBestBoundDown <- masterObj
      currentBestBoundUp <- sum(get.linProb(mProb, type='objective') * heuristicSolution)
      branchedVars <- NULL
      
      ### check if we have already the optimum solution (without rounding errors)
      runInd <- TRUE
      if ( abs(masterObj-currentBestBoundUp) < 0.1 ) {
        runInd <- FALSE
      } else {
        ### fixing variables
        if ( fixVariables == TRUE & currentBestBoundUp >= currentBestBoundDown ) {
          if ( verbose ) {
            cat("fixing variables...\n")
          }
          fixedVars <- calc.linProb(mProb, type='fixVariables', input=list(currentBestBoundDown, currentBestBoundUp, primSupps))
          if (length(fixedVars) > 0 ) {
            if ( verbose ) {
              cat("--> setting",length(fixedVars),"variables to 0!\n")
            }
            bounds <- get.linProb(mProb, type='bounds')
            bounds$upper$val[fixedVars] <- 0
            mProb <- set.linProb(mProb, type='bounds', input=bounds)
          }
        }
        
        ### constraint pool initialization
        problemPool <- list()
        problemPool[[1]] <- init.cutList(type='empty', input=list(nrCols=nrVars))
        
        ### solving
        selectFirst <- TRUE
        LPL <- get.problemInstance(problemInstance, type="LPL")
        UPL <- get.problemInstance(problemInstance, type="UPL")
        SPL <- get.problemInstance(problemInstance, type="SPL")
        
        weights <- get.problemInstance(problemInstance, type='weight')
        LB <- weights - get.problemInstance(problemInstance, type="lb")
        UB <- get.problemInstance(problemInstance, type="ub") - weights
        nrConstraints <- length(get.linProb(aProb, type='objective')) - 2*length(weights)
        
        ### initialize constants (probably function-parameters later)
        selectFirst <- FALSE
        
        ### TODO: stop procedure after given time or nr or solutions..
        iter <- 0
      }
      
      while( runInd ) {
        iter <- iter + 1
        selectInd <- ifelse(selectFirst==TRUE, 1, length(problemPool))
        newCuts <- init.cutList(type='empty', input=list(nrCols=nrVars))
        AttProbDown <- AttProbUp <- rep(NA, length(primSupps))
        status <- NULL
        for ( i in 1:length(primSupps) ) {
          cellInd <- primSupps[i]
          limits <- c(LPL[cellInd], UPL[cellInd], SPL[cellInd])
          
          # we need to solve the attackers problem for each sensitive cell twice
          if ( limits[3] != 0 ) {
            # solveAttackerProblem (upper bound)
            rhs <- rep(0, length(get.linProb(aProb, type='rhs')))
            rhs[cellInd] <- 1
            aProb <- set.linProb(aProb, type='rhs', list(rhs))
            aProb <- set.linProb(aProb, type='objective', input=list(c(weights + UB*xi, -(weights-xi*LB), rep(0, nrConstraints))))
            up <- calc.linProb(aProb, type='solveProblem', input=list(solver))
            
            ### solveAttackerProblem (lower bound)
            aProb <- set.linProb(aProb, type='rhs', input=list(-1*rhs))
            down <- calc.linProb(aProb, type='solveProblem', input=list(solver))
            
            AttProbDown[i] <- -down$optimum
            AttProbUp[i] <- up$optimum
            status <- c(status, down$status, up$status)
            #cat('limits ( origValue=',weights[cellInd],') : [',AttProbDown[i],':',AttProbUp[i],']\n')
            
            alpha.down <- down$solution[1:nrVars]
            alpha.up <- up$solution[1:nrVars]
            
            beta.down <- down$solution[(nrVars+1):(2*nrVars)]
            beta.up <- up$solution[(nrVars+1):(2*nrVars)]
          } else {
            # solve attackers problem (minimize)
            if ( limits[1] != 0 ) {
              rhs <- rep(0, length(get.linProb(aProb, type='rhs')))
              rhs[cellInd] <- -1
              aProb <- set.linProb(aProb, type='rhs', input=list(rhs))
              down <- calc.linProb(aProb, type='solveProblem', input=list(solver))
              AttProbDown[i] <- -down$optimum
            }
            # solve attackers problem (maximize)
            if ( limits[2] != 0 ) {
              rhs <- rep(0, length(get.linProb(aProb, type='rhs')))
              rhs[cellInd] <- 1
              aProb <- set.linProb(aProb, type='rhs', input=list(rhs))
              aProb <- set.linProb(aProb, type='objective', input=list(c(weights + UB*xi, -(weights-xi*LB), rep(0, nrConstraints))))
              up <- calc.linProb(aProb, type='solveProblem', input=list(solver))
              AttProbUp[i] <- up$optimum
            }
          }
          # SPL
          if ( limits[3] != 0 & AttProbUp[i] - AttProbDown[i] < limits[3] ) {
            status <- c(status, down$status, up$status)
            alpha.down <- down$solution[1:nrVars]
            alpha.up <- up$solution[1:nrVars]
            beta.down <- down$solution[(nrVars+1):(2*nrVars)]
            beta.up <- up$solution[(nrVars+1):(2*nrVars)]
            
            v <- (alpha.down+alpha.up)*UB + (beta.down+beta.up)*LB
            v[which(is.zero(v))] <- 0
            if ( any(v != 0) )
              newCuts <- set.cutList(newCuts, type='addCompleteConstraint', input=list(init.cutList(type='singleCut', input=list(vals=v, dir=">=", rhs=limits[3]))))
          } else  {
            if ( limits[1] != 0 & freqs[primSupps[i]] - AttProbDown[i] < limits[1] ) { # LPL
              status <- c(status, down$status)
              alpha.down <- down$solution[1:nrVars]
              beta.down <- down$solution[(nrVars+1):(2*nrVars)]
              
              v <- alpha.down*UB + beta.down*LB
              v[which(is.zero(v))] <- 0
              if ( any(v != 0) )
                newCuts <- set.cutList(newCuts, type='addCompleteConstraint', input=list(init.cutList(type='singleCut', input=list(vals=v, dir=">=", rhs=limits[1]))))
            }
            if ( limits[2] != 0 & AttProbUp[i] - freqs[primSupps[i]] < limits[2] ) { # UPL
              status <- c(status, up$status)
              alpha.up <- up$solution[1:nrVars]
              beta.up <- up$solution[(nrVars+1):(2*nrVars)]
              
              v <- alpha.up*UB + beta.up*LB
              v[which(is.zero(v))] <- 0
              if ( any(v != 0) )
                newCuts <- set.cutList(newCuts, type='addCompleteConstraint', input=list(init.cutList(type='singleCut', input=list(vals=v, dir=">=", rhs=limits[2]))))
            }
          }
          #cat('limits ( origValue=',weights[cellInd],') : [',AttProbDown[i],':',AttProbUp[i],']\n')
        }
        
        if ( get.cutList(newCuts, type='nrConstraints') > 0 ) {
          # strengthen cuts
          if ( verbose ) {
            cat("strengthening the cuts and adding",get.cutList(newCuts, type='nrConstraints'),"new derived cuts to master problem...\n")
          }
          newCuts <- calc.cutList(newCuts, type='strengthen', input=list())
          mProb <- set.linProb(mProb, type='addCompleteConstraint', input=list(newCuts))
        }
        ### check for duplicated constraints
        indRem <- which(duplicated(cbind(as.matrix(get.linProb(mProb, type='constraints')), get.linProb(mProb, type='direction'), get.linProb(mProb, type='rhs'))))
        if ( length(indRem) > 0 ) {
          #if ( verbose ) {
          # cat("removing",length(indRem),"duplicated constraints...\n")
          #}
          mProb <- set.linProb(mProb, type='removeCompleteConstraint', input=list(indRem))
        }
        ### bridgeless inequalities only at root-node ####
        bridgelessCandidates <- setdiff(which(xi == 1), primSupps)
        if ( iter == 1 & length(bridgelessCandidates) > 0 ) {
          brCuts <- init.cutList(type='empty', input=list(nrCols=nrVars))
          if ( verbose ) {
            cat("adding",length(bridgelessCandidates),"bridgeless ineqalities!\n")
          }
          for ( z in seq_along(bridgelessCandidates) ) {
            bridgelessInd <- bridgelessCandidates[z]
            ### solveAttackerProblem (upper bound)
            rhs <- rep(0, length(get.linProb(aProb, type='rhs')))
            rhs[bridgelessInd] <- 1
            aProb <- set.linProb(aProb, type='rhs', input=list(rhs))
            aProb <- set.linProb(aProb, type='objective', input=list(c(weights + UB*xi, -(weights-xi*LB), rep(0, nrConstraints))))
            up <- calc.linProb(aProb, type='solveProblem', input=list(solver))
            
            ### solveAttackerProblem (lower bound)
            aProb <- set.linProb(aProb, type='rhs', input=list(-1*rhs))
            down <- calc.linProb(aProb, type='solveProblem', input=list(solver))
            
            alpha.down <- down$solution[1:nrVars]
            alpha.up <- up$solution[1:nrVars]
            
            beta.down <- down$solution[(nrVars+1):(2*nrVars)]
            beta.up <- up$solution[(nrVars+1):(2*nrVars)]
            
            brIneq <- (alpha.down+alpha.up)*UB + (beta.down+beta.up)*LB
            brIneq[is.zero(brIneq)] <- 0
            brIneq[brIneq > 0] <- 1
            brIneq[bridgelessInd] <- -1
            brCuts <- set.cutList(brCuts, type='addCompleteConstraint', input=list(init.cutList(type='singleCut', input=list(vals=brIneq, dir=">=", rhs=0))))
          }
          if ( get.cutList(brCuts, type='nrConstraints') > 0 ) {
            mProb <- set.linProb(mProb, type='addCompleteConstraint', input=list(brCuts))
          }
        }
        
        mProbWorking <- mProb
        # eventually update the lower bound...
        tmpSolution <- calc.linProb(mProbWorking, type='solveProblem', input=list(solver))
        tmpObj <- tmpSolution$optimum
        if ( tmpObj > currentBestBoundDown & tmpObj <= currentBestBoundUp ) {
          currentBestBoundDown <- tmpObj
        }
        if ( abs(currentBestBoundUp - currentBestBoundDown) < 1 ) {
          # optimal solution found!
          break
        }
        
        if ( get.cutList(problemPool[[selectInd]], type='nrConstraints') > 0 ) {
          if ( verbose ) {
            cat("adding",get.cutList(problemPool[[selectInd]], type='nrConstraints'),"constraints to the master problem...\n")
          }
          mProbWorking <- set.linProb(mProb, type='addCompleteConstraint', input=list(problemPool[[selectInd]]))
        }
        
        if ( verbose ) {
          cat("solving the master problem with",length(get.linProb(mProbWorking, type='rhs')),"constraints...\n")
        }
        masterSolution <- calc.linProb(mProbWorking, type='solveProblem', input=list(solver))
        masterObj <- masterSolution$optimum
        xi <- masterSolution$solution
        xi[is.zero(xi)] <- 0
        xi[is.one(xi)] <- 1
        
        if ( verbose ) {
          cat("best-bounds: [",currentBestBoundDown,":",currentBestBoundUp,"] and objVal =",masterObj,"with sum(xi)=",sum(xi),"\n")
        }
        #cat("current boundUp =",currentBestBoundUp,"and objVal =",masterObj,"with sum(xi)=",sum(xi),"\n")
        
        ### again fixing variables
        #if ( fixVariables == TRUE ) {
        # newFixedVars <- calc.linProb(mProb, type='fixVariables', input=list(currentBestBoundDown, currentBestBoundUp, primSupps))
        # if ( !all(newFixedVars) %in% fixedVars ) {
        #   cat("setting",length(newFixedVars),"variables to 0!\n")
        #   bounds <- get.linProb(mProb, type='bounds')
        #   bounds$upper$val[newFixedVars] <- 0
        #   mProb <- set.linProb(mProb, type='bounds', input=bounds)
        # }
        #}
        
        ### checking if we can prune the current node
        prune <- FALSE
        pruneReason <- NULL
        # a) valid (protected) integer solution
        if ( all(is.wholenumber(xi)) && calc.problemInstance(problemInstance, type='isProtectedSolution', input=list(input1=AttProbDown, input2=AttProbUp)) ) {
          prune <- TRUE
          pruneReason <- c(pruneReason, "V") # valid
        }
        # b) infeasibility
        if ( sum(status) != 0 ) {
          prune <- TRUE
          pruneReason <- c(pruneReason, "I") # infeasible
        }
        # c) bounds
        if ( approx == TRUE ) {
          if ( currentBestBoundUp - masterObj <= currentBestBoundUp * (approxPerc/100) ) {
            prune <- TRUE
            pruneReason <- c(pruneReason, "B") # bounds
          }
          if ( masterObj - currentBestBoundDown < 0 ) {
            prune <- TRUE
            pruneReason <- c(pruneReason, "B") # bounds
          }
          
        } else {
          if ( abs(masterObj - currentBestBoundUp) <= 0.01 ) {
            prune <- TRUE
            pruneReason <- c(pruneReason, "B") # bounds
          }
          if ( masterObj - currentBestBoundDown < 0 ) {
            prune <- TRUE
            pruneReason <- c(pruneReason, "B") # bounds
          }
        }
        
        if ( prune == TRUE ) {
          pruneReason <- unique(pruneReason) # remove eventually 2-'Bs'
          if ( length(pruneReason) == 2 ) {
            if ( pruneReason[1] == "V" & pruneReason[2]=="B" ) {
              #cat("Integer-Lösung gefunden, aber es gibt bereits eine bessere -> Pruning by Bounds!\n")
              if ( masterObj < currentBestBoundUp ) {
                pruneReason <- "V"
              } else {
                pruneReason <- "B"
              }
            }
          }
          if ( length(pruneReason) > 1 ) {
            stop("Error: only one pruning reason possible!\n")
          }
          if ( pruneReason == "V") {
            solutions[[length(solutions)+1]] <- as.integer(xi)
            if ( masterObj < currentBestBoundUp ) {
              if ( verbose ) {
                cat("new best integer solution (objval=",masterObj,") found!:\n")
              }
              currentBestBoundUp <- masterObj
              bestSolution <- as.integer(xi)
            }
          }
          #if ( pruneReason == "I") {
          # cat("pruning because of infeasibility!\n")
          #}
          #if ( pruneReason == "B") {
          # cat("pruning because of known bounds!\n")
          #}
          problemPool[[selectInd]] <- NULL
          if ( verbose ) {
            cat("pruning the current node: reason=",pruneReason,"!. Still",length(problemPool),"nodes in the pool!\n")
          }
        } else {
          ## 2) Branching: wir erweitern den ProblemPool und löschen dann den aktuellen Node
          branchedVars <- get.simpleTriplet(get.cutList(problemPool[[selectInd]], type='constraints'), type='colInd', input=list())
          branchVar <- getBranchingVariable(xi, branchedVars, noBranchVars)
          
          if ( length(branchVar) == 1 ) {
            cl <- problemPool[[selectInd]]
            v <- rep(0, nrVars)
            v[branchVar] <- 1
            c1 <- set.cutList(cl, type='addCompleteConstraint', input=list(init.cutList(type='singleCut', input=list(vals=v, dir="==", rhs=0))))
            c2 <- set.cutList(cl, type='addCompleteConstraint', input=list(init.cutList(type='singleCut', input=list(vals=v, dir="==", rhs=1))))
            problemPool[[length(problemPool)+1]] <- c1
            problemPool[[length(problemPool)+1]] <- c2; rm(cl)
            
            # now we can prune the current node
            problemPool[[selectInd]] <- NULL
            rm(c1,c2)
            if ( verbose ) {
              cat("branching was required. Problem pool has now",length(problemPool),"nodes!\n")
            }
          } else {
            if ( verbose ) {
              cat("no further branching possible! all branching variables tried!\n")
            }
            problemPool[[selectInd]] <- NULL
          }
        }
        
        timeSpent <- as.numeric(floor(difftime(Sys.time(), startTime, units="mins")))
        #cat("timeSpent:"); print(timeSpent)
        
        if ( length(problemPool)==0 ) {
          runInd <- FALSE
        } else {
          if ( !is.null(timeLimit) && timeSpent > timeLimit && length(solutions) > 0 ) {
            runInd <- FALSE
            timeStop <- TRUE
          }
          if ( !is.null(timeLimit) && timeSpent > timeLimit && length(solutions) == 0 ) {
            if ( verbose ) {
              cat("Result: the time-limit was reached and no (heuristic) solution could be generated!\n")
            }
            return(object)
          }
        }
      }
      
      secondSupps <- setdiff(which(bestSolution==1), primSupps)
      objVarHeuristic <- sum(get.linProb(mProb, type='objective') * heuristicSolution)
      if ( timeStop==FALSE ) {
        if ( currentBestBoundUp == objVarHeuristic ) {
          if ( verbose ) {
            cat('Result: the heuristic solution was already optimal and has',length(secondSupps),'secondary suppressions!\n')
          }
        } else {
          improvement <- 100 - (100 / objVarHeuristic ) * currentBestBoundUp
          if ( verbose ) {
            cat('Result: the heuristic solution was improved by',format(improvement, digits=2, nsmall=2),'% and has',length(secondSupps),'secondary suppressions!!\n ')
          }
        }
      } else {
        if ( verbose ) {
          cat("Result: we are returning a possibly non-optimal solution with",length(secondSupps),"secondary suppressions because of argument 'timeLimit'!\n")
        }
      }
      if ( length(secondSupps) > 0 ) {
        problemInstance <- set.problemInstance(problemInstance, type='sdcStatus', input=list(index=secondSupps, values=rep("x", length(secondSupps))))
      }
      out <- new("sdcProblem",
        dataObj=get.sdcProblem(object, type='dataObj'),
        dimInfo=dimInfo,
        problemInstance=problemInstance,
        partition=get.sdcProblem(object, type='partition'),
        startI=get.sdcProblem(object, type='startI'),
        startJ=get.sdcProblem(object, type='startJ'),
        indicesDealtWith=get.sdcProblem(object, type='indicesDealtWith'),
        elapsedTime=get.sdcProblem(object, type='elapsedTime')+(proc.time()-time.start)[3]
      )
      return(out)
    }
    
    if ( type == 'anonWorker' ) {
      timeLimit <- input$timeLimit
      verbose <- input$verbose
      save <- input$save
      
      if( save == TRUE ) {
        files <- NULL
      }
      
      start.time <- proc.time()
      pI <- get.sdcProblem(object, type='problemInstance')
      sdcStatusBegin <- get.problemInstance(pI, type='sdcStatus')
      primSupps <- primSuppsOrig <- get.problemInstance(pI, type='primSupps')
      
      indexPool <- numeric()
      allStrIDs <- get.problemInstance(pI, type='strID')
      
      if ( input$method == 'OPT' ) {
        object <- set.sdcProblem(object, type='elapsedTime', input=list(get.sdcProblem(object, type='elapsedTime') + (proc.time()-start.time)[3]))
        
        if ( input$useC == TRUE ) {
          result <- csp_cpp(sdcProblem=object, attackonly=FALSE, verbose=input$verbose)
        } else {
          result <- calc.sdcProblem(object=object, type='cutAndBranch', input=input)
        }
        
        if ( save==TRUE ) {
          fn <- paste(input$method,"_Object-Final.RData", sep="")
          save(object, file=fn)
        }
        return(result)
      }
      
      # HITAS or HYPERCUBE
      # check where we should start (saved file)
      partition <- get.sdcProblem(object, type='partition')
      startI <- get.sdcProblem(object, type='startI')
      startJ <- get.sdcProblem(object, type='startJ')
      
      if ( startI !=1 | startJ != 1 ) {
        maxI <- partition$nrGroups
        if ( startJ < length(partition$indices[[startI]]) ) {
          startJ <- startJ+1
        } else {
          startJ <- 1
          startI <- startI+1
        }
      }     
      
      if ( input$method == 'HITAS' ) {
        for ( i in startI:(partition$nrGroups) ) {
          object <- set.sdcProblem(object, type='startJ', input=list(1)) # reset j before updating i
          object <- set.sdcProblem(object, type='startI', input=list(i))
          #indexPool <- get.sdcProblem(object, type='indicesDealtWith')
          if ( i == 1 ) {
            indexPool <- c()
          } else {
            indexPool <- sort(unique(unlist(partition$indices[1:(i-1)])))
          }
          if ( verbose ) {
            cat("dealing with Group",i,"|",length(partition$groups),"\n")
          }
          ind <- partition$indices[[i]]
          
          beginJ <- ifelse(i==startI, startJ, 1)
          for ( j in beginJ:(length(ind)) ) {
            object <- set.sdcProblem(object, type='startJ', input=list(j))
            currentIndices <- ind[[j]] # within complete table
            
            ### cells with status 'u' or 'x' exist
            pI <- get.sdcProblem(object, type='problemInstance')
            if ( any(get.problemInstance(pI, type='sdcStatus')[currentIndices] %in% c("u","x")) & length(currentIndices) > 1 ) {
              if ( verbose ) {
                cat("starting to solve problem",j,"/",length(ind),"in this group!\n")
              }
              ### if we have cells with "u" or "x" we need to protect
              ### the corresponding subtable
              ### reduce problemInstance
              probNew <- calc.sdcProblem(object, type='reduceProblem', input=list(currentIndices))
              pI.new <- get.sdcProblem(probNew, type='problemInstance')
              
              ### is it necessary to protect the table??
              currentPrimSupps <- primSupps[!is.na(match(primSupps, currentIndices ))]
              
              ### indices that have already been inner cells in
              ### tables dealt earlier
              # FIXME: save indexpool somehow
              indicesDealtWith <- which(currentIndices %in% indexPool) #in current current subproblem
              
              ### fix marginal-cells
              ### --> its-suppression state must not change!
              currentPattern <- get.problemInstance(get.sdcProblem(probNew, type='problemInstance'), type='sdcStatus')
              
              introducedSupps <- indicesDealtWith[which(currentPattern[indicesDealtWith] == "x")]
              if ( length(introducedSupps) > 0 ) {
                ### secondary suppressions from upper tables
                pI.new <- set.problemInstance(pI.new, type='sdcStatus', input=list(index=introducedSupps, values=rep("u", length(introducedSupps))))
                
                ### temporarily change LPL, UPL, SPL for these cells
                LPL.current <- get.problemInstance(pI.new, type='LPL')[introducedSupps]
                UPL.current <- get.problemInstance(pI.new, type='UPL')[introducedSupps]
                SPL.current <- get.problemInstance(pI.new, type='SPL')[introducedSupps]
                
                pI.new <- set.problemInstance(pI.new, type='LPL', input=list(index=introducedSupps, values=rep(0, length(introducedSupps))))
                pI.new <- set.problemInstance(pI.new, type='UPL', input=list(index=introducedSupps, values=rep(0, length(introducedSupps))))
                pI.new <- set.problemInstance(pI.new, type='SPL', input=list(index=introducedSupps, values=rep(0.1, length(introducedSupps))))
              }
              
              ### force non-suppression of cells that have already been dealt with
              indForced <- indicesDealtWith[which(currentPattern[indicesDealtWith] == "s")]
              if ( length(indForced) > 0 ) {
                pI.new <- set.problemInstance(pI.new, type='sdcStatus', input=list(index=indForced, values=rep("z", length(indForced))))
              }
              probNew <- set.sdcProblem(probNew, type='problemInstance', input=list(pI.new))
              
              ### solving the problem
              if ( input$useC == TRUE ) {
                probNew <- csp_cpp(sdcProblem=probNew, attackonly=FALSE, verbose=input$verbose)
              } else {
                probNew <- calc.sdcProblem(object=probNew, type='cutAndBranch', input=input)
              }
              
              ### update sdcStatus
              status <- get.problemInstance(get.sdcProblem(probNew, type='problemInstance'), type='sdcStatus')
              
              pI <- get.sdcProblem(object, type='problemInstance')
              if ( length(indForced) > 0 ) {
                status[indForced] <- "s"
              }
              if ( length(introducedSupps) > 0 ) {
                status[introducedSupps] <- "x"
                pI <- set.problemInstance(pI, type='LPL', input=list(index=currentIndices[introducedSupps], values=LPL.current))
                pI <- set.problemInstance(pI, type='UPL', input=list(index=currentIndices[introducedSupps], values=UPL.current))
                pI <- set.problemInstance(pI, type='SPL', input=list(index=currentIndices[introducedSupps], values=SPL.current))
              }
              
              pI <- set.problemInstance(pI, type='sdcStatus', input=list(index=currentIndices, values=status))
              object <- set.sdcProblem(object, type='problemInstance', input=list(pI))
            }
            if ( save == TRUE ) {
              if ( verbose ) {
                cat("saving object after i=",i,"and j=",j,"\n")
              }
              fn <- paste(input$method,"_Object_",i,"-",j,".RData", sep="")
              files <- c(files, fn)
              save(object, file=fn)
              
              # removing old files
              if ( length(files) > 1 ) {
                sapply(rev(files)[-1], file.remove)
                files <- files[length(files)]
              }
            }
          }
          ### update indices that we have already dealt with
          object <- set.sdcProblem(object, type='indicesDealtWith', input=list(unique(c(indexPool, currentIndices))))
        }       
      }
      
      if ( input$method == 'HYPERCUBE' ) {
        runInd <- TRUE
        nrRuns <- 1
        while ( runInd == TRUE ) { 
          cat("run=",nrRuns,"\n")
          
          tmpSupps <- c(get.problemInstance(get.sdcProblem(object, 'problemInstance'), 'primSupps'), 
              get.problemInstance(get.sdcProblem(object, 'problemInstance'), 'secondSupps'))
          forcedCells <- get.problemInstance(get.sdcProblem(object, 'problemInstance'), 'forcedCells')
      
          for ( i in startI:(partition$nrGroups) ) {
            object <- set.sdcProblem(object, type='startJ', input=list(1)) # reset j before updating i
            object <- set.sdcProblem(object, type='startI', input=list(i))
            if ( verbose ) {
              cat("dealing with Group",i,"|",length(partition$groups),"\n")
            }
            ind <- partition$indices[[i]]
            
            beginJ <- ifelse(i==startI, startJ, 1)
            for ( j in beginJ:(length(ind)) ) {
              object <- set.sdcProblem(object, type='startJ', input=list(j))
              currentIndices <- ind[[j]] # within complete table
              
              ### cells with status 'u' or 'x' exist
              pI <- get.sdcProblem(object, type='problemInstance')
              # when using HYPERCUBE: we only check primary suppressions because we
              # temporarily set secondary suppressions to "u"
              if ( any(get.problemInstance(pI, type='sdcStatus')[currentIndices] == "u") & length(currentIndices) > 1 ) {
                if ( verbose ) {
                  cat("starting to solve problem",j,"(total=",length(ind),") in this group!\n")
                }
                
                ### if we have cells with "u",  we need to protect
                ### the corresponding subtable
                ### reduce problemInstance
                probNew <- calc.sdcProblem(object, type='reduceProblem', input=list(currentIndices))
                pI.new <- get.sdcProblem(probNew, type='problemInstance')

                ### is it necessary to protect the table??
                currentPrimSupps <- primSupps[!is.na(match(primSupps, currentIndices ))]
                
                probNew <- set.sdcProblem(probNew, type='problemInstance', input=list(pI.new))
                
                ### solving the problem
                probNew <- calc.sdcProblem(object=probNew, type='ghmiter', input=input)
                
                ### update sdcStatus
                status <- get.problemInstance(get.sdcProblem(probNew, type='problemInstance'), type='sdcStatus')
                pI <- get.sdcProblem(object, type='problemInstance')
                pI <- set.problemInstance(pI, type='sdcStatus', input=list(index=currentIndices, values=status))
                object <- set.sdcProblem(object, type='problemInstance', input=list(pI))
              }
              if ( save == TRUE ) {
                if ( verbose ) {
                  cat("saving object after i=",i,"and j=",j,"\n")
                }
                fn <- paste(input$method,"_Object_",i,"-",j,".RData", sep="")
                files <- c(files, fn)
                save(object, file=fn)
                
                # removing old files
                if ( length(files) > 1 ) {
                  sapply(rev(files)[-1], file.remove)
                  files <- files[length(files)]
                }
              }
            }
          } 
        
          ### protect secondary suppressions ### 
          pI <- get.sdcProblem(object, 'problemInstance')
          allSupps <- c(
            get.problemInstance(pI, type='primSupps'), 
            get.problemInstance(pI, type='secondSupps')
          )
          newSupps <- setdiff(allSupps, tmpSupps)

          pI <- get.sdcProblem(object, type='problemInstance')
          
          nrVars <- length(get.problemInstance(pI, 'freq'))
          if ( length(newSupps) == 0 ) {
            runInd <- FALSE
            newSdcStatus <- rep('s', length=nrVars)
            newSdcStatus[forcedCells] <- 'z'
            newSdcStatus[tmpSupps] <- 'x'
            newSdcStatus[primSuppsOrig] <- 'u'
            pI <- set.problemInstance(pI, type='sdcStatus', input=list(index=1:nrVars, values=newSdcStatus))
            object <- set.sdcProblem(object, type='problemInstance', input=list(pI))
          } else {
            
            newSdcStatus <- rep('s', length=nrVars)
            newSdcStatus[forcedCells] <- 'z'
            newSdcStatus[tmpSupps] <- 'x'
            newSdcStatus[newSupps] <- 'u'
            pI <- set.problemInstance(pI, type='sdcStatus', input=list(index=1:nrVars, values=newSdcStatus))
            object <- set.sdcProblem(object, type='problemInstance', input=list(pI))
            nrRuns <- nrRuns + 1
          }       
        } # while loop
      }
      return(object)
    }
    
    if ( type == 'ghmiter' ) {
      time.start <- proc.time()
      protectionLevel <- input$protectionLevel
      suppMethod <- input$suppMethod
      suppAdditionalQuader <- input$suppAdditionalQuader
      verbose <- input$verbose
      
      pI <- get.sdcProblem(object, type='problemInstance')
      strIDs <- get.problemInstance(pI, type='strID')
      strInfo <- get.dimInfo(get.sdcProblem(object, type='dimInfo'), type='strInfo')
      
      sdcStatus <- get.problemInstance(pI, type='sdcStatus')
      cellsToProtect <- get.problemInstance(pI, type='primSupps')
      
      freqs <- get.problemInstance(pI, type='freq')
      
      # calc infomation on inner|marginal cells
      cellInfo <- get.sdcProblem(object, type='innerAndMarginalCellInfo')
      
      # replaces f.recodeIndexVars
      indexList <- lapply(1:length(strInfo), function(x) { mySplit(strIDs, strInfo[[x]][1]:strInfo[[x]][2]) } )
      
      for ( i in 1:length(cellsToProtect) ) {
        if ( verbose ) {
          cat("--> Cell",i,"|",length(cellsToProtect)," (ID:",strIDs[cellsToProtect[i]],")...")
        }
        
        diagObj <- calc.sdcProblem(object, type='ghmiter.diagObj', input=list(cellsToProtect[i], indexList, FALSE))   
    
        # calculate required information using diagObj
        infoObj <- calc.sdcProblem(object, type='ghmiter.calcInformation', input=list(diagObj, indexList, protectionLevel, FALSE))
        
        if ( !is.null(infoObj) & length(infoObj) == 0 ) {
          diagObj <- calc.sdcProblem(object, type='ghmiter.diagObj', input=list(cellsToProtect[i], indexList, TRUE))    
          infoObj <- calc.sdcProblem(object, type='ghmiter.calcInformation', list(diagObj, indexList, protectionLevel, TRUE))
          
          if ( !is.null(infoObj) & length(infoObj) == 0 ) {
            stop("ghmiter::: error - could not find sensible cube!\n")
          }
          warning("Cell with Freq=0 has been selected as partner in suppression pattern!\n")
        }
        
        # cellToProtect used from diagObj$cellToProtect
        suppObj <- calc.sdcProblem(object, type='ghmiter.selectQuader', input=list(infoObj, input))
        
        if ( !is.null(suppObj) ) {
          object <- calc.sdcProblem(object, type='ghmiter.suppressQuader', input=suppObj)
          
          # additional quader needs to be found
          # only if it is not a single value in the margins
          # and the cube includes cells with frequency=1
          if ( suppAdditionalQuader==TRUE & suppObj$indikatorSingleItems==TRUE & !(cellsToProtect[i] %in% cellInfo$indexTotCells) ) {
            # find additional cube that does not contain the single cells
            object <- calc.sdcProblem(object, type='ghmiter.suppressAdditionalQuader', input=list(diagObj, infoObj, suppObj, input))
          }
        }
        if ( verbose ) {
          cat("[DONE]\n")
        }
      }
      object <- set.sdcProblem(object, type='elapsedTime', input=list(get.sdcProblem(object, type='elapsedTime')+(proc.time()-time.start)[3]))
      return(object)    
    }
  
    if ( type == 'preprocess') {
      solver <- input$solver
      problemInstance <- get.sdcProblem(object, type='problemInstance')
      if ( !get.problemInstance(problemInstance, type='hasPrimSupps') ) {
        return(object)
      }
      dimInfo <- get.sdcProblem(object, type='dimInfo')
      nrVars <- get.problemInstance(problemInstance, type='nrVars')
      freqs <- get.problemInstance(problemInstance, type='freq')
      primSupps <- get.problemInstance(problemInstance, type='primSupps')
      
      LPL <- get.problemInstance(problemInstance, type="LPL")[primSupps]
      UPL <- get.problemInstance(problemInstance, type="UPL")[primSupps]
      SPL <- get.problemInstance(problemInstance, type="SPL")[primSupps]
      
      weights <- get.problemInstance(problemInstance, type='weight')
      HIGH <- LOW <- weights[primSupps]
      
      # order of calculations
      myOrder <- order(sapply(1:length(primSupps), function(x) { max(SPL[x], (LPL+UPL)[x]) }), decreasing=TRUE)
      LB <- weights - get.problemInstance(problemInstance, type="lb")
      UB <- get.problemInstance(problemInstance, type="ub") - weights
      
      xi <- get.problemInstance(problemInstance, type='suppPattern')
      res <- calc.multiple(type='makeAttackerProblem', input=list(objectA=problemInstance, objectB=dimInfo))
      aProb <- res$aProb
      validCuts <- res$newCutsMaster
      
      nrConstraints <- length(get.linProb(aProb, type='objective')) - 2*length(weights)
      
      #validCuts <- init.cutList(type='empty', input=list(nrCols=nrVars))
      for ( i in myOrder ) {
        if ( i %% 10 == 1 && input$verbose ) {
          cat("preprocessing variable",i,"|",length(myOrder),"...\n")
        }
        cellInd <- primSupps[i]
        limits <- c(LPL[i], UPL[i], SPL[i])
        
        # solveAttackerProblem (upper bound)
        rhs <- rep(0, length(get.linProb(aProb, type='rhs')))
        rhs[cellInd] <- 1
        aProb <- set.linProb(aProb, type='rhs', input=list(rhs))
        aProb <- set.linProb(aProb, type='objective', input=list(c(weights + UB*xi, -(weights-xi*LB), rep(0, nrConstraints))))
        up <- calc.linProb(aProb, type='solveProblem', input=list(solver))
        if ( up$status != 0 ) {
          stop("unsolvable problem (up)!\n")
        }
        calcUp <- up$optimum
        HIGH[i] <- max(HIGH[i], calcUp)
        LOW[i] <- min(LOW[i], calcUp)
        
        # solveAttackerProblem (lower bound)
        aProb <- set.linProb(aProb, type='rhs', input=list(-1*rhs))
        down <- calc.linProb(aProb, type='solveProblem', input=list(solver))
        if ( down$status != 0 ) {
          stop("unsolvable problem (down)!\n")
        }
        calcDown <- -down$optimum
        HIGH[i] <- max(HIGH[i], calcDown)
        LOW[i] <- min(LOW[i], calcDown)
        
        alpha.down <- down$solution[1:nrVars]
        alpha.up <- up$solution[1:nrVars]
        
        beta.down <- down$solution[(nrVars+1):(2*nrVars)]
        beta.up <- up$solution[(nrVars+1):(2*nrVars)]
        
        if ( limits[1] != 0 ) { # LPL
          v <- alpha.down*UB + beta.down*LB
          v[which(is.zero(v))] <- 0
          if ( any(v > 0) )
            validCuts <- set.cutList(validCuts, type='addCompleteConstraint', input=list(init.cutList(type='singleCut', input=list(vals=v, dir=">=", rhs=limits[1]))))
        }
        if ( limits[2] != 0 ) { # UPL
          v <- alpha.up*UB + beta.up*LB
          v[which(is.zero(v))] <- 0
          if ( any(v > 0) )
            validCuts <- set.cutList(validCuts, type='addCompleteConstraint', input=list(init.cutList(type='singleCut', input=list(vals=v, dir=">=", rhs=limits[2]))))
        }
        if ( limits[3] != 0 & calcUp - calcDown < limits[3] ) { # SPL
          v <- (alpha.down+alpha.up)*UB + (beta.down+beta.up)*LB
          v[which(is.zero(v))] <- 0
          if ( any(v > 0) )
            validCuts <- set.cutList(validCuts, type='addCompleteConstraint', input=list(init.cutList(type='singleCut', input=list(vals=v, dir=">=", rhs=limits[3]))))
        }
      }
      
      if ( get.cutList(validCuts, type='nrConstraints') > 0 ) {
        validCuts <- calc.cutList(validCuts, type='strengthen', input=list())
        
        setZeroUPL <- which(freqs[primSupps]+UPL <= HIGH) # -> set UPL = 0
        setZeroLPL <- which(freqs[primSupps]-LPL >= LOW) # -> set LPL = 0
        setZeroSPL <- which(HIGH-LOW >= SPL ) # -> set SPL = 0
        
        if ( length(setZeroUPL) > 0 ) {
          UPL[setZeroUPL] <- 0
          problemInstance <- set.problemInstance(problemInstance, type='UPL', input=list(index=primSupps, values=UPL))
        }
        if ( length(setZeroLPL) > 0 ) {
          LPL[setZeroLPL] <- 0
          problemInstance <- set.problemInstance(problemInstance, type='LPL', input=list(index=primSupps, values=LPL))
        }
        if ( length(setZeroSPL) > 0 ) {
          SPL[setZeroSPL] <- 0
          problemInstance <- set.problemInstance(problemInstance, type='SPL', input=list(index=primSupps, values=SPL))
        }       
      }

      object <- set.sdcProblem(object, type='problemInstance', input=list(problemInstance))
      return(list(sdcProblem=object, aProb=aProb, validCuts=validCuts))     
    }
    
    if ( type == 'cellID' ) {
      para.names <- input[[1]]
      para.codes <- input[[2]]
      para.verbose <- input[[3]]
      
      pI <- get.sdcProblem(object, type='problemInstance')
      dimInfoObj <- get.sdcProblem(object, type='dimInfo')
      
      vNames <- get.dimInfo(dimInfoObj, type='varName')
      vIndex <- get.dimInfo(dimInfoObj, type='posIndex')
      
      indexVar <- match(para.names, vNames)
      strInfo <- get.dimInfo(dimInfoObj, type='strInfo')
      dimInfo <- get.dimInfo(dimInfoObj, type='dimInfo')
      
      ### original Codes berechnen
      codesDefault <- lapply(1:length(strInfo), function(x) { mySplit(get.problemInstance(pI, type='strID'), strInfo[[x]][1]:strInfo[[x]][2]) } )
      codesOrig <- list()
      for ( i in 1:length(codesDefault) ) {
        codesOrig[[i]] <- calc.dimVar(object=dimInfo[[i]], type="matchCodeOrig", input=codesDefault[[i]])
      }
      
      if ( length(input) != 3 ) {
        stop("calc.sdcProblem (type=cellID):: length of argument 'input' must equal 3!\n")
      }
      if ( length(para.names) != length(para.codes) ) {
        stop("calc.sdcProblem (type=cellID):: check argument 'input'!\n")
      }
      if ( !all(para.names %in% vNames) ) {
        stop("calc.sdcProblem (type=cellID):: check variable names in 'input[[1]]'!\n")
      }
      if ( !is.logical(para.verbose) ) {
        stop("calc.sdcProblem (type=cellID):: argument in 'input[[3]]' must be logical!\n")
      }
      
      cellID <- 1:get.problemInstance(pI, type='nrVars')
      for ( i in seq_along(para.names) ) {
        cellID <- intersect(cellID, which(!is.na(match(as.character(codesOrig[[indexVar[i]]]), para.codes[i]))))
      }
      if ( length(cellID) != 1) {
        stop("calc.sdcProblem (type=cellID):: check argument 'input' -> 0 or > 1 cells identified!\n")
      }
      return(cellID)      
    }
  
    if ( type == 'finalize' ) {
      time.start <- proc.time()
      
      pI <- get.sdcProblem(object, type='problemInstance')
      dI <- get.sdcProblem(object, type='dimInfo')
      levelObj <- get.dimInfo(dI, type='dimInfo')
      strInfo <- get.dimInfo(dI, type='strInfo')
      
      sdcStatus <- get.problemInstance(pI, type='sdcStatus')
      
      nrNonDuplicatedCells <- get.problemInstance(pI, type='nrVars')
      nrPrimSupps <- length(which(sdcStatus == 'u'))
      nrSecondSupps <- length(which(sdcStatus == 'x'))
      nrPublishableCells <- length(which(sdcStatus %in% c('z','s')))
      
      #####################################
      ### Merge Bezeichnungen und Codes ###
      #####################################
      codesOriginal <- list()
      strIDs <- get.problemInstance(pI, type='strID')
      for ( i in seq_along(levelObj) ) {
        codesDefault <- mySplit(strIDs, strInfo[[i]][1]:strInfo[[i]][2])
        codesOriginal[[i]] <- calc.dimVar(object=levelObj[[i]], type='matchCodeOrig', input=codesDefault)
      }
      out <- data.frame(codesOriginal);
      colnames(out) <- get.dimInfo(dI, type='varName')
      out$Freq <- get.problemInstance(pI, type='freq')
      
      numVars <- get.problemInstance(pI, type='numVars')
      if ( !is.null(numVars) ) {
        data.obj <- get.sdcProblem(object,"dataObj")
        nV <- as.data.frame(numVars)
        colnames(nV) <- colnames(get.dataObj(data.obj,"rawData"))[get.dataObj(data.obj,"numVarInd")]
        out <- cbind(out, nV)
      }
      out$sdcStatus <- get.problemInstance(pI, type='sdcStatus')
      
      #############################
      ### Duplikate hinzufuegen ###
      #############################
      hasDups <- sapply(1:length(levelObj), function(x) { get.dimVar(levelObj[[x]], type='hasDuplicates') })
      if ( any(hasDups == TRUE) ) {
        for ( i in which(hasDups==TRUE) ) {
          dups <- get.dimVar(levelObj[[i]], type='dups')
          dupsUp <- get.dimVar(levelObj[[i]], type='dupsUp')
          
          runInd <- TRUE
          while ( runInd ) {
            for ( j in 1:length(dups) ) {
              sub <- subset(out, out[,i]==dupsUp[j])
              if ( nrow(sub) > 0 ) {
                sub[,i] <- dups[j]
                out <- rbind(out, sub)
                dups <- dups[-j]
                dupsUp <- dupsUp[-j]
              }
            }
            if ( length(dups) == 0 )
              runInd <- FALSE
          }
        }
      }
      rownames(out) <- NULL
      safeObj <- new("safeObj",
        finalData=out,
        dimInfo=get.sdcProblem(object, type='dimInfo'),
        nrNonDuplicatedCells=nrNonDuplicatedCells,
        nrPrimSupps=nrPrimSupps,
        nrSecondSupps=nrSecondSupps,
        nrPublishableCells=nrPublishableCells,
        suppMethod=input$method,
        elapsedTime=get.sdcProblem(object, type='elapsedTime') + (proc.time()-time.start)[3]
      )
      return(safeObj)     
    }
    
    if ( type == 'ghmiter.diagObj' ) {
      cellToProtect <- input[[1]]
      indexList <- input[[2]]
      allowZeros <- input[[3]]
      
      if ( length(cellToProtect) != 1 ) {
        stop("makeDiagObj:: length of 'cellToProtect' must equal 1!\n")
      }
      pI <- get.sdcProblem(object, type='problemInstance')
      nrVars <- get.problemInstance(pI, type='nrVars')
      if ( !cellToProtect %in% 1:nrVars ) {
        stop("makeDiagObj:: 'cellToProtect' must be >= 1 and <= ",nrVars,"!\n")
      }
      if ( !cellToProtect %in% 1:nrVars ) {
        stop("makeDiagObj:: 'cellToProtect' must be >= 1 and <= ",nrVars,"!\n")
      }
      sdcStatus <- get.problemInstance(pI, type='sdcStatus')
      
      if ( allowZeros == TRUE ) {
        indZeros <- which(sdcStatus=="z" & get.problemInstance(pI, type='freq')==0)
        if ( length(indZeros ) > 0 ) {
          sdcStatus[indZeros] <- "s"
        }
      }
      
      indToProtect <- sapply(indexList, function(x) { x[cellToProtect] } )
      diagIndices <- NULL
      for ( i in seq_along(indexList) ) {
        if ( length(unique(indexList[[i]])) > 1  ) {
          if ( is.null(diagIndices) ) {
            diagIndices <- which(indexList[[i]] != indexList[[i]][cellToProtect])
          } else {
            diagIndices <- intersect(diagIndices, which(indexList[[i]] != indexList[[i]][cellToProtect]))
          }
        }
      }
      diagIndices <- diagIndices[which(sdcStatus[diagIndices]!="z")]
      
      if ( length(diagIndices) == 0 )
        diagIndices <- NULL
      
      supp <- list()
      supp$cellToProtect <- cellToProtect
      supp$indToProtect <- indToProtect
      supp$diagIndices <- diagIndices
      return(supp)      
    }
    
    if ( type == 'ghmiter.calcInformation' ) {
      # calculate some info on a given quader (normalization,...)
      calcQInfo <- function(g, d) {
        if ( length(g) != length(d) ) {
          stop("calcQInfo: 'g' and 'd' must have equal length!\n")
        }
        numberIndexVars <- length(g)
        
        quader <- expand(lapply(1:numberIndexVars, function(x) { c(g[x],d[x]) }  ), vector=FALSE)
        
        quader <- matrix(unlist(quader), length(quader[[1]]), length(quader))
        quader <- quader[!duplicated(quader),,drop=FALSE]
        quader <- lapply(1:ncol(quader), function(x) quader[,x])
        
        ### normquader
        normQ <- list()
        for ( i in seq_along(quader) ) {
          normQ[[i]] <- rep(0, length(quader[[1]]))
          normQ[[i]][which(quader[[i]]==g[i])] <- 1
        }
        
        ### g|u indication?
        indexing <- rep("g", length(quader[[1]]))
        indexing[which(apply(matrix(unlist(normQ),length(quader[[1]]),numberIndexVars),1,sum) %%2 != 0)] <- "u"
        return(list(quader=quader, normQ=normQ, indexing=indexing))     
      }
      
      diagObj <- input[[1]]
      indexList <- input[[2]]
      protectionLevel <- input[[3]]     
      allowZeros <- input[[4]]
      
      # TODO: error checking
      pI <- get.sdcProblem(object, type='problemInstance')
      freqs <- get.problemInstance(pI, type='freq')
      sdcStatus <- get.problemInstance(pI, type='sdcStatus')
      
      if ( allowZeros == TRUE ) {
        indZeros <- which(sdcStatus=="z" & get.problemInstance(pI, type='freq')==0)
        if ( length(indZeros ) > 0 ) {
          sdcStatus[indZeros] <- "s"
        }
      }
      
      ### relevant Indices ### TODO FIXME CHECK!!!!
      relevantIndices <- which(sapply(indexList, function(x) { length(unique(x)) } ) > 1)
      
      resultObj <- list()
      # FIXME: What is with 1-dimensional data?
      limit <- length(diagObj$diagIndices)
      for ( z in 1:limit ) {
        g <- diagObj$indToProtect
        d <- sapply(indexList, function(x) { x[diagObj$diagIndices[z]] } )
        qInfo <- calcQInfo(g, d)
        
        # 2) position (indices==qPosition) of current quader in subTabObj
        valsQ <- pasteStrVec(unlist(qInfo$quader), length(qInfo$quader), coll=NULL)
        qPosition <- match(valsQ, get.problemInstance(pI, type='strID'))
        suppStatus <- sdcStatus[qPosition]
        if ( !any(suppStatus == "z") ) {
          # 3) calculate various information about the selected quader (infoQuader)
          # 3.1) how many values would need to be suppressed for this quader
          indNonSupp <- which(sdcStatus[qPosition] == "s" & sdcStatus[qPosition] != "u")
          nrAdditionalSupps <- length(indNonSupp)
          
          # 3.2) whats the amount of information which needs to be suppressed?
          sumAdditionalSuppsFreq <- sum(freqs[qPosition[indNonSupp]])
          
          # 3.3) does the quader contains other single cells except for
          # the primary suppressed value (diaObj$cellToPretect) to check?
          # subIndices = current quader without primary suppressed cell to check
          indSingleItems <- setdiff(which(freqs[qPosition]==1),1)
          singleItems <- NULL
          indikatorSingleItems <- FALSE
          if( length(indSingleItems) >= 1 ) {
            indikatorSingleItems <- TRUE
            singleItems <- indSingleItems
          }
          
          # 3.5) is the quader protected enough? (protectionLevel)
          # we need to check for interval-protection only if protectionLevel > 0
          schutzInd <- TRUE
          schutz <- protectionLevel
          # FIXME: S|P (what to do with "x" that are temporarily "u"?)
          if( protectionLevel > 0 ) {
            if ( !all(qInfo$indexing =="u") ) {
              range <- min(freqs[qPosition[which(qInfo$indexing =="u")]], na.rm=TRUE) + min(freqs[qPosition[which(qInfo$indexing =="g")]], na.rm=TRUE)
              X <- freqs[diagObj$cellToProtect]
              if( X == 0 ) {
                tmpInd <- which(sdcStatus[qPosition] != "u" & freqs[qPosition] != 0)
                
                if( length(tmpInd) > 0 ) {
                  # TODO: this needs testing !!! (page 60, repsilber)
                  if( range <= min(freqs[tmpInd]) ) {
                    schutzInd <- FALSE
                    protectionLevel <- 0
                  }
                }
              }
              else {
                schutz <- (100*range) / X
                if ( schutz < protectionLevel )
                  schutzInd <- FALSE
              }
            }
          }
          
          # 4) return results
          # in this case, the cell is already protected, so we can stop!
          
          # allowZeros==TRUE: we have not found patterns without zeros
          # so we do not care for an 'optimal' solution that does not exist anyway
          if ( allowZeros == TRUE ) {
            if ( length(resultObj) == 100 ) {
              return(resultObj)
              break
            }
          } else {
            if( nrAdditionalSupps == 0 & schutzInd == TRUE & indikatorSingleItems == FALSE  ) {
              return(erg = NULL)
              break
            }
          }
          
          resultObj[[length(resultObj)+1]] <- list(
            quaderStrID = valsQ,
            indexing = qInfo$indexing,
            qPosition = qPosition,
            nrAdditionalSupps=nrAdditionalSupps,
            sumAdditionalSuppsFreq = sumAdditionalSuppsFreq,
            indikatorSingleItems = indikatorSingleItems,
            singleItems = singleItems,
            schutz = schutz,
            schutzInd = schutzInd
          )
        }
      }
      return(resultObj)   
    }
    
    if ( type == 'ghmiter.suppressQuader' ) {
      pI <- get.sdcProblem(object, type='problemInstance')
      sdcStatus <- get.problemInstance(pI, type='sdcStatus')
      
      suppIndex <- setdiff(input$qPosition, input$qPosition[which(sdcStatus[input$qPosition]=="u")])
      pI <- set.problemInstance(pI, type='sdcStatus', input=list(index=suppIndex, values=rep("x", length(suppIndex))))
      object <- set.sdcProblem(object, type='problemInstance', input=list(pI))
      return(object)      
    }
    
    if ( type == 'ghmiter.selectQuader' ) {
      infoObj <- input[[1]]
      suppMethod <- input[[2]]$suppMethod
      verbose <- input[[2]]$verbose
      
      sdcStatus <- get.problemInstance(get.sdcProblem(object, type='problemInstance'), type='sdcStatus')
      
      relevantIndices <- as.numeric(unlist(lapply(infoObj, '[', 'schutzInd'))[1])
      
      # already protected
      if ( is.null(infoObj) ) {
        suppObj <- NULL
        return(suppObj)
      }
      
      # not protected yet
      # which elements of iqsInfo are NULL?
      nullElements <- which(unlist(lapply(lapply(infoObj, '[[', 'qPosition'), function(x) { length(x) } )) == 0)
      if ( length(nullElements) > 0 ) {
        infoObj <- infoObj[-nullElements]
      }
      
      # put iqs together so that we can choose the optimal suppression scheme
      qIndexNr <- 1:length(infoObj)
      nrAdditionalSupps <- as.numeric(unlist(lapply(infoObj, '[', 'nrAdditionalSupps')))
      sumAdditionalSuppsFreq <- as.numeric(unlist(lapply(infoObj, '[', 'sumAdditionalSuppsFreq')))
      indikatorSingleItems <- as.logical(unlist(lapply(infoObj, '[', 'indikatorSingleItems')))
      schutz <- as.numeric(unlist(lapply(infoObj, '[', 'schutz')))
      schutzInd <- as.logical(unlist(lapply(infoObj, '[', 'schutzInd')))
      schutzInd <- as.logical(unlist(lapply(infoObj, '[', 'schutzInd')))
      
      possQuaders <- data.frame(qIndexNr, nrAdditionalSupps, sumAdditionalSuppsFreq, indikatorSingleItems, schutz, schutzInd)
      
      # are there any suppression schemes satisfying the necessary interval protection?
      indexIntervallOk <- FALSE
      if ( any(possQuaders$schutzInd==TRUE) )
        indexIntervallOk <- TRUE
      
      # do suppression schemes exist that do not contain single values?
      # these are preferred suppression schemes.
      existNonSingles <- FALSE
      if ( any(possQuaders$indikatorSingleItems == FALSE) ) {
        existNonSingles <- TRUE
        possQuaders <- possQuaders[possQuaders$indikatorSingleItems==FALSE,,drop=FALSE]
      }
      
      if( indexIntervallOk ) {
        if ( min(possQuaders$nrAdditionalSupps) > 0 & verbose == TRUE) {
          cat("# additional secondary Supps:", min(possQuaders$nrAdditionalSupps)," ")
        }
        if ( suppMethod == "minSupps" ) {
          possQuaders <- possQuaders[which(possQuaders$nrAdditionalSupps == min(possQuaders$nrAdditionalSupps)),,drop=FALSE]
          possQuaders <- possQuaders[which(possQuaders$sumAdditionalSuppsFreq == min(possQuaders$sumAdditionalSuppsFreq)),,drop=FALSE]
        }
        if ( suppMethod == "minSum" ) {
          possQuaders <- possQuaders[which(possQuaders$sumAdditionalSuppsFreq == min(possQuaders$sumAdditionalSuppsFreq)),,drop=FALSE]
          possQuaders <- possQuaders[which(possQuaders$nrAdditionalSupps == min(possQuaders$nrAdditionalSupps)),,drop=FALSE]
        }
        if ( suppMethod == "minSumLogs" ) {
          possQuaders <- possQuaders[which(possQuaders$sumAdditionalSuppsFreq == min(log(1+possQuaders$sumAdditionalSuppsFreq))),,drop=FALSE]
          possQuaders <- possQuaders[which(possQuaders$nrAdditionalSupps == min(possQuaders$nrAdditionalSupps)),,drop=FALSE]
        }
        
        # finally choose the suppression scheme
        possQuaders <- possQuaders[1,]
        suppObj <- infoObj[[possQuaders$qIndexNr]]
      }
      # problem: no suppression scheme is satisfying the
      # required interval protection
      else {
        # all cells in this subtable are already suppressed
        # -> everything is ok
        if( all(sdcStatus=="u" | sdcStatus == "x") )
          suppObj <- NULL
        
        # no suppression scheme satisfies the required interval protection
        # the suppression pattern with the max. protection level is selected
        else {
          possQuaders <- possQuaders[which(possQuaders$schutz == max(possQuaders$schutz)),]
          possQuaders <- possQuaders[which(possQuaders$nrAdditionalSupps == min(possQuaders$nrAdditionalSupps)),]
          possQuaders <- possQuaders[which(possQuaders$sumAdditionalSuppsFreq == min(possQuaders$sumAdditionalSuppsFreq)),]
          possQuaders <- possQuaders[1,]
          suppObj <- infoObj[[possQuaders$qIndexNr]]
        }
      }
      return(suppObj)     
    }
    
    if ( type == 'ghmiter.suppressAdditionalQuader') {
      diagObj <- input[[1]]
      infoObj <- input[[2]]
      suppObj <- input[[3]]
      suppMethod <- input[[4]]$suppMethod
      verbose <- input[[4]]$suppMethod
      
      ### Task: find quader (from) infoObj with following restrictions
      freqs <- get.problemInstance(get.sdcProblem(object, type='problemInstance'), type='freq')
      
      # - must not be suppObj itself
      cellToProtect <- diagObj$cellToProtect
      suppIndicesOrig <- suppObj$qPosition
      
      # the additional quader must non contain these indices
      # the singletons in the original suppressed pattern
      prohibitedIndices <- setdiff(suppIndicesOrig[which(freqs[suppIndicesOrig] == 1)],cellToProtect)
      
      # the possible indices
      possIndices <- lapply(infoObj, function(x) { x$qPosition } )
      
      # do the indices of the possible patterns contain any of the prohibited cells
      res <- sapply(1:length(possIndices), function(x) { any(prohibitedIndices %in% possIndices[[x]]) } )
      if ( all(res == TRUE ) ) {
        #warning("no additional cube could be found!\n")
        infoObj <- NULL
      } else {
        ind <- which(res==TRUE)
        infoObj <- infoObj[-ind]
      }
      
      suppObjNew <- calc.sdcProblem(object, type='ghmiter.selectQuader', input=list(infoObj, input[[4]]))
      if ( !is.null(suppObjNew) ) {
        object <- calc.sdcProblem(object, type='ghmiter.suppressQuader', input=suppObjNew)
      }
      return(object)      
    }
    
    if ( type == 'contributingIndices' ) {
      strID <- input[[1]]
      dataObj <- get.sdcProblem(object, type='dataObj')
      dimInfoObj <- get.sdcProblem(object, type='dimInfo')
      dimInfo <- get.dimInfo(dimInfoObj, type='dimInfo')
      pI <- get.sdcProblem(object, type='problemInstance')
      
      if ( !strID %in% get.problemInstance(pI, type='strID') ) {
        stop("calc.sdcProblem (type=contributingIndices):: strID not found in the current problem!\n")
      }
      dims <- lapply(dimInfo, function(x) { get.dimVar(x, type='dims') } )
      indexVec <- which(get.dimInfo(dimInfoObj, type='strID')==strID)
      # some (sub)totals need to be considered
      if( length(indexVec) == 0 ) {
        levInfo <- list()
        for ( z in 1:length(dimInfo) ) {
          subLevel <- substr(strID, get.dimInfo(dimInfoObj, type='strInfo')[[z]][1], get.dimInfo(dimInfoObj, type='strInfo')[[z]][2])
          if ( sum(as.numeric(subLevel)) == 0 ) {
            levInfo[[z]] <- sort(unique(unlist(dims[[z]])))
          } else {
            orderInd <- unlist(lapply(dims[[z]], function(x) { match(subLevel, x)}))
            if( min(orderInd, na.rm=TRUE) == 1 ) {
              levInfo[[z]] <- dims[[z]][[which(orderInd==1)]]
            } else {
              levInfo[[z]] <- subLevel
            }
          }
        }
        cellIndex <- pasteStrVec(unlist(expand.grid(levInfo)), length(levInfo))
        indexVec <- which(get.dimInfo(dimInfoObj, type='strID') %in% cellIndex)
      }
      return(indexVec)      
    }
    
    if ( type == 'reduceProblem' ) {
      x <- object
      y <- input[[1]]
      
      pI <- get.sdcProblem(x, type='problemInstance')
      dimInfo <- get.sdcProblem(x, type='dimInfo')
      strInfo <- strInfoOrig <- get.dimInfo(dimInfo, type='strInfo')
      
      if ( length(y) < 1 ) {
        stop("calc.sdcProblem (type=reduceProblem):: length of argument 'y' < 1!\n")
      }
      if ( !all(y %in% 1:get.problemInstance(pI, type='nrVars')) ) {
        stop("reduceProblem (type=reduceProblem):: elements of indices y does not match with problem size!\n")
      }
      
      newDims <- lapply(1:length(strInfo), function(x) { substr(get.problemInstance(pI, type='strID')[y], strInfo[[x]][1], strInfo[[x]][2]) } )
      newDims2 <- lapply(1:length(newDims), function(x) { sort(unique(newDims[[x]])) } )
      newDimsOrigCodes <- lapply(1:length(newDims), function(k) {
        calc.dimVar(object=dimInfo@dimInfo[[k]], type='matchCodeOrig', input=newDims2[[k]])     
      })  
      
      lenNewDims <- sapply(newDims2, length)-1
      codesNew <- lapply(1:length(newDims), function(x) { c("@", rep("@@", lenNewDims[x])) } )
      
      dimInfoOld <- lapply(1:length(newDims2),  function(x) { init.dimVar(input=list(input=data.frame(codesNew[[x]], newDims2[[x]]), vName=paste('V',x,sep="")) ) } )
      dimInfoNew <- lapply(1:length(newDims2),  function(x) { init.dimVar(input=list(input=data.frame(codesNew[[x]], newDimsOrigCodes[[x]]), vName=paste('V',x,sep="")) ) } )
      
      new.codes <- lapply(1:length(newDims), function(x) {
        dimInfoOld[[x]]@codesDefault[match(newDims[[x]], dimInfoOld[[x]]@codesOriginal)]    
      })
      pI@strID <- pasteStrVec(unlist(new.codes), length(newDims))
      pI@Freq <- get.problemInstance(pI, type='freq')[y]
      if ( !is.null(get.problemInstance(pI, type='w') )) {
        pI@w <- get.problemInstance(pI, type='w')[y]
      }
      numVars <- as.list(get.problemInstance(pI, type='numVars'))
      if ( length(numVars) > 0 ) {
        for ( j in 1:length(numVars) ) {
          pI@numVars[[j]] <- numVars[[j]][y]
        }
      }
      pI@lb <- get.problemInstance(pI, type='lb')[y]
      pI@ub <- get.problemInstance(pI, type='ub')[y]
      pI@LPL <- get.problemInstance(pI, type='LPL')[y]
      pI@UPL <- get.problemInstance(pI, type='UPL')[y]
      pI@SPL <- get.problemInstance(pI, type='SPL')[y]
      pI@sdcStatus <- get.problemInstance(pI, type='sdcStatus')[y]
      x@dimInfo@dimInfo <- dimInfoNew
      
      ## strInfo
      info <- c(0, cumsum(sapply(1:length(codesNew), function(x) { sum(sapply(table(codesNew[[x]]), nchar)) } )))
      for ( i in 2:length(info) ) {
        strInfo[[i-1]] <- c(info[i-1]+1, info[i] )
      }
      x@dimInfo@strInfo <- strInfo
      
      #tmpRes <- lapply(1:length(newDims), function(x) { substr(get.problemInstance(pI, type='strID'), strInfoOrig[[x]][1], strInfoOrig[[x]][2]) }  )
      #tmpRes <- lapply(1:length(tmpRes), function(x) {
      # tmpRes[[x]] <- calc.dimVar(object=dimInfoOld[[x]], type='matchCodeDefault', input=tmpRes[[x]])
      #})
      #pI@strID <- pasteStrVec(unlist(tmpRes), length(tmpRes))

      x <- set.sdcProblem(x, type='problemInstance', input=list(pI))
      validObject(x)
      return(x)     
    }
  
    if ( type == 'genStructuralCuts' ) {
      pI <- get.sdcProblem(object, type='problemInstance')
      dimInfoObj <- get.sdcProblem(object, type='dimInfo')
      partition <- calc.multiple(type='makePartitions', input=list(objectA=pI, objectB=dimInfoObj))
      
      dimInfo <- get.dimInfo(dimInfoObj, type='dimInfo')
      nrLevels <- length(dimInfo)
      nrVars <- get.problemInstance(pI, type='nrVars')
      primSupps <- get.problemInstance(pI, type='primSupps')
      strIDs <- get.problemInstance(pI, type='strID')
      indices <- partition$indices
      weights <- get.problemInstance(pI, type='weight')
      requiredCuts <- init.cutList(type='empty', input=list(nrCols=nrVars))
      strInfo <- get.dimInfo(dimInfoObj, type='strInfo')
      x <- rep(0, nrVars)
      
      for ( z in seq_along(primSupps) ) {
        pSupp <- primSupps[z]
        currentPrimSupp <- strIDs[pSupp]
        matchInd <- unlist(lapply(1:length(indices), function(x) { 
          lapply(1:length(indices[[x]]), function(y) { 
            if ( !all(is.na(match(pSupp, indices[[x]][[y]]))) ) {c(x,y)}
          }) 
        }))   
        if ( any(is.na(matchInd)) ) {
          stop('elements of matchInd must not be NA!\n')
        }     
        
        splitMatchInd <- split(matchInd, rep(1:(length(matchInd)/2), each=2))
        
        for ( u in 1:length(splitMatchInd) ) {
          matchInd <- splitMatchInd[[u]]
          nrPow <- nrLevels - length(which(as.numeric(unlist(strsplit(partition$groups[[matchInd[1]]],"-")))==1))
          v1 <- v2 <- x
          index <- indices[[matchInd[1]]][[matchInd[2]]]
          v1[index] <- weights[index]
          v2[index] <- 1
          lim <- sum(sort(weights[index])[1:(2^nrPow)])
          if ( any(v1 != 0) ) {
            requiredCuts <- set.cutList(requiredCuts, type='addCompleteConstraint', input=list(init.cutList(type='singleCut', input=list(vals=v1, dir=">=", rhs=lim))))
          }
          if ( any(v2 != 0) ) {
            requiredCuts <- set.cutList(requiredCuts, type='addCompleteConstraint', input=list(init.cutList(type='singleCut', input=list(vals=v2, dir=">=", rhs=(2^nrPow)))))
          }
          
          ### Todo: at least 2 suppressions in each dimension
          ### there is some error here! -> TODO: CHECK: FIXME!
          #for ( i in 1:length(dimInfo) ) {
          # lO <- dimInfo[[i]]
          # splitList <- lapply(strInfo[i], function(k) { seq(k[1], k[2]) } )
          # subStringToFix <- mySplitIndicesList(currentPrimSupp, splitList)
          # f <- mySplitIndicesList(strIDs, splitList)
          # 
          # index <- which(f == subStringToFix)
          # v3 <- x
          # v4 <- x
          # v3[index] <- weights[index]
          # v4[index] <- 1  
          # lim <- sum(sort(weights[index])[1:2])
          # if ( any(v3 != 0) ) {
          #   if ( !is.na(lim) ) {
          #     requiredCuts <- set.cutList(requiredCuts, type='addCompleteConstraint', input=list(init.cutList(type='singleCut', input=list(vals=v3, dir=">=", rhs=lim))))
          #   }
          # }
          # if ( any(v4 != 0) ) {
          #   requiredCuts <- set.cutList(requiredCuts, type='addCompleteConstraint', input=list(init.cutList(type='singleCut', input=list(vals=v4, dir=">=", rhs=2))))
          # }
          #}        
        }
      } 
      
      dupRows <- get.simpleTriplet(get.cutList(requiredCuts, type='constraints'), type='duplicatedRows', input=list())
      if ( length(dupRows) > 0 ) {
        requiredCuts <- set.cutList(requiredCuts, type='removeCompleteConstraint', input=list(dupRows))
      }
      return(requiredCuts)    
    }
  }
)

#' summarize \code{\link{sdcProblem-class}} objects
#'
#' extract and show relevant information stored in \code{\link{sdcProblem-class}} objects
#'
#' @aliases summary,sdcProblem-method
#' @rdname summary-method
#' @export
#' @docType methods
setMethod(f='summary', signature='sdcProblem',
  definition=function(object, ...) {
    pI <- get.sdcProblem(object, type="problemInstance")
    dO <- get.sdcProblem(object, type="dataObj")
    dI <- get.sdcProblem(object, type="dimInfo")
    if ( get.dataObj(dO, type="isMicroData") ) {
      cat("The raw data contains micro data!")
      if ( length(pI@numVars) > 0 ) {
        cat("--> the use of dominance rules for primary suppressions is possible!")
      }
      cat("\n")
    } else {
      cat("The raw data contain pre-aggregated (tabular) data!\n")
    } 
    
    nrcells <- get.problemInstance(pI, type="nrVars")
    dim_names <- get.dimInfo(dI, type="varName")
 
    cat("\nThe complete table to protect consists of",nrcells,"cells and has",length(dim_names),"spanning variables.") 
    
    cat("\nThe distribution of\n")
    cat("- primary unsafe (u)\n")
    cat("- secondary suppressed (x)\n")
    cat("- forced to publish (z) and\n")
    cat("- selectable for secondary suppression (s) cells is shown below:\n")
    print(table(get.problemInstance(pI, type="sdcStatus")))
    
    nr_tables <- get.sdcProblem(object, type="partition")$nrTables
    cat("\nIf this table is protected with heuristic methods, a total of",nr_tables,"has (sub)tables must be considered!\n")    

  }
)

