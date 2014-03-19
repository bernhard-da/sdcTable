#' @aliases get.linProb,linProb,character-method
#' @rdname get.linProb-method
setMethod(f='get.linProb', signature=c('linProb', 'character'),
  definition=function(object, type) { 
    if ( !type %in% c('constraints', 'direction', 'rhs', 'objective', 'types', 'bounds') ) {
      stop("get.cutList:: argument 'type' is not valid!\n")
    }
    if ( type == 'constraints' ) {
      return(object@constraints)
    }
    if ( type == 'direction' ) {
      return(object@direction)
    }
    if ( type == 'rhs' ) {
      return(object@rhs)
    }       
    if ( type == 'objective' ) {
      return(object@objective)
    } 
    if ( type == 'types' ) {
      return(object@types)
    }     
    if ( type == 'bounds' ) {
      return(list(upper=object@boundsUpper, lower=object@boundsLower))
    }     
  }
)

#' @aliases set.linProb,linProb,character,list-method
#' @rdname set.linProb-method
setMethod(f='set.linProb', signature=c('linProb', 'character', 'list'),
  definition=function(object, type, input) {
    if ( !type %in% c('objective', 'direction', 'rhs', 'types', 
        'removeCompleteConstraint', 'addCompleteConstraint',
        'bounds', 'constraints') ) {
      stop("set.linProb:: check argument 'type'!\n")
    } 
    if ( type == 'objective' ) {      
      object@objective <- input[[1]]    
    }
    if ( type == 'direction' ) {
      object@direction <- input[[1]]      
    }
    if ( type == 'rhs' ) {
      object@rhs <- input[[1]]      
    } 
    if ( type == 'types' ) {
      object@types <- input[[1]]      
    } 
    if  ( type == 'removeCompleteConstraint' ) {
      input <- input[[1]]
      if ( !all(input %in% 1:length(get.linProb(object, type='rhs'))) ) {
        stop("set.linProb:: elements of argument 'input' must be >=1 and <=",length(get.linProb(object, type='rhs')),"!\n")
      }
      object@constraints <- calc.simpleTriplet(get.linProb(object, type='constraints'), type='removeRow', input=list(input))      
      object@direction <- get.linProb(object, type='direction')[-input]
      object@rhs <- get.linProb(object, type='rhs')[-input] 
    }
    if ( type == 'addCompleteConstraint' ) {
      input <- input[[1]]
      if ( get.simpleTriplet(get.linProb(object, type='constraints'), type='nrCols', input=list()) != get.simpleTriplet(get.cutList(input, type='constraints'), type='nrCols', input=list()) ) {
        stop("set.linProb:: nrCols of 'object' and 'input' differ!\n")
      }
      if ( get.cutList(input, type='nrConstraints') > 0 ) {
        con <- get.cutList(input, type='constraints')
        for ( k in 1:get.simpleTriplet(con, type='nrRows', input=list()) ) {
          x <- get.simpleTriplet(con, type='getRow', input=list(k))
          object@constraints <- calc.simpleTriplet(get.linProb(object, type='constraints'), type='addRow', input=list(index=get.simpleTriplet(x, type='colInd', input=list()), values=get.simpleTriplet(x, type='values', input=list())))     
        }
        object@direction <- c(get.linProb(object, type='direction'), get.cutList(input, type='direction'))
        object@rhs <- c(get.linProb(object, type='rhs'), get.cutList(input, type='rhs'))      
      }
    }     

    if ( type == 'bounds' ) {
      # FIXME: check bounds input (lower|upper,...)
      object@boundsLower <- input$lower   
      object@boundsUpper <- input$upper 
    } 
    
    if ( type == 'constraints' ) {
      object@constraints <- input[[1]]  
    } 
    
    
    validObject(object)
    return(object)        
  }
)

#' @aliases calc.linProb,linProb,character,list-method
#' @rdname calc.linProb-method
setMethod(f='calc.linProb', signature=c('linProb', 'character', 'list'),
  definition=function(object, type, input) {
    if ( !type %in% c('solveProblem', 'fixVariables') ) {
      stop("calc.linProb:: check argument 'type'!\n")
    } 
    
    if ( type == 'solveProblem' ) {
      solver <- input[[1]]
    
      if ( !solver %in% c("glpk", "symphony", "lpSolve") ) {
        stop("'solver' needs to be eiter 'glpk', 'lpSolve' or 'symphony!\n'")
      }   
      if ( solver == "glpk" ) {
        sol <- my.Rglpk_solve_LP(
          get.linProb(object, type='objective'), 
          get.linProb(object, type='constraints'),
          get.linProb(object, type='direction'), 
          get.linProb(object, type='rhs'), 
          get.linProb(object, type='types'), 
          max = FALSE, 
          bounds=get.linProb(object, type='bounds'),
          verbose = FALSE)
      }
      if ( solver == "lpSolve" ) {
        stop("solving with 'lpSolve' not yet available!\n")
      }
      if ( solver == "symphony" ) {
        #sol <- Rsymphony_solve_LP(
        # get.linProb(object, type='objective'),
        # get.linProb(object, type='constraints'),
        # get.linProb(object, type='direction'), 
        # get.linProb(object, type='rhs'), 
        # bounds=get.linProb(object, type='bounds'), #bounds
        # get.linProb(object, type='types'),  
        # max = FALSE)
        stop("solving with 'symphony' not yet available!\n")
      } 
      if ( solver == "cplex" ) {
        #directionOrig <- get.linProb(object, type='direction'), 
        #sense <- rep(NA, length(directionOrig))      
        #sense[directionOrig=="=="] <- "E"
        #sense[directionOrig=="<="] <- "L"
        #sense[directionOrig==">="] <- "G"
        #sol <- Rcplex(
        # get.linProb(object, type='objective'), 
        # get.linProb(object, type='constraints'), 
        # get.linProb(object, type='rhs'), 
        # Qmat = NULL,
        # lb = 0, 
        # ub = 1, 
        # control = list(),
        # objsense = "min", 
        # sense = sense, 
        # vtype = get.linProb(object, type='types'), 
        # n = 1)  
        stop("solving with 'cplex' not yet available!\n")
      }   
      return(sol)
    }
    
    if ( type == 'fixVariables' ) {
      lb <- input[[1]]
      ub <- input[[2]]
      primSupps <- input[[3]]

      if ( length(lb) != 1 | length(ub) != 1 ) {
        stop("calc.linProb (type==fixVariables):: length of arguments 'lb' and 'ub' must equal 1!\n")
      }
      if ( !ub > lb ) {
        stop("calc.linProb (type==fixVariables):: arguments 'ub' must be >= argument 'lb'!\n")
      } 
      
      con <- get.linProb(object, type='constraints')
      rhs <- get.linProb(object, type='rhs')
      dir <- get.linProb(object, type='direction')
      obj <- get.linProb(object, type='objective')
      bounds <- get.linProb(object, type='bounds')
      nrVars <- get.simpleTriplet(con, type='nrCols', input=list())
      
      my.lp <- make.lp(0, nrVars)
      
      set.objfn(my.lp, obj)
      set.bounds(my.lp, upper = bounds$upper$val)
      set.bounds(my.lp, lower = bounds$lower$val)
      
      for ( i in 1:get.simpleTriplet(con, type='nrRows', input=list()) ) {
        r <- get.simpleTriplet(con, type='getRow', input=list(i))
        cols <- get.simpleTriplet(r, type='colInd', input=list())
        vals <- get.simpleTriplet(r, type='values', input=list())
        dd <- ifelse(dir[i]=="==", "=", dir[i])
        add.constraint(my.lp, vals, dd, rhs[i], indices=cols)
      }   
      solve(my.lp)
      get.objective(my.lp)
      
      dual <- get.dual.solution(my.lp)
      dual <- dual[(2+length(get.rhs(my.lp))):length(dual)]
      if ( length(dual) != nrVars ) {
        stop("calc.linProb (type==fixVariables):: length of arguments does not match!\n")
      }
      freqs <- obj
      sol <- get.variables(my.lp)
      sol[is.zero(sol)] <- 0
      sol[is.one(sol)] <- 1
      
      ### calculate reduced costs from dual solution 
      reducedCosts <- freqs
      dualInd <- which(dual!=0)
      if ( length(dualInd) > 0 ) {
        reducedCosts[dualInd] <- sapply(dualInd, function(x) { min(freqs[x], dual[x])} )
      }
      
      bas <- which(sol != 0)
      reducedCosts[bas] <- 0
      ### end calculation of reduced costs
      
      ### which variables could be set to zero?
      indSetToZero <- which(lb + reducedCosts >= ub)
      
      # do not fix primary suppressions to zero!
      indSetToZero <- indSetToZero[-which(indSetToZero %in% primSupps)]
      return(indSetToZero)
    } 
  }
)

