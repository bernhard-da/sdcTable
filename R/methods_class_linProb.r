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
      object@constraints <- c_remove_row(get.linProb(object, type='constraints'), input=list(input))
      object@direction <- get.linProb(object, type='direction')[-input]
      object@rhs <- get.linProb(object, type='rhs')[-input]
    }
    if ( type == 'addCompleteConstraint' ) {
      input <- input[[1]]
      if ( g_nr_cols(get.linProb(object, type='constraints')) != g_nr_cols(g_constraints(input)) ) {
        stop("set.linProb:: nrCols of 'object' and 'input' differ!\n")
      }
      if ( g_nr_constraints(input) > 0 ) {
        con <- g_constraints(input)
        for ( k in 1:g_nr_rows(con) ) {
          x <- g_row(con, input=list(k))
          object@constraints <- c_add_row(get.linProb(object, type='constraints'), input=list(index=g_col_ind(x), values=g_values(x)))
        }
        object@direction <- c(get.linProb(object, type='direction'), g_direction(input))
        object@rhs <- c(get.linProb(object, type='rhs'), g_rhs(input))
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
      nrVars <- g_nr_cols(con)

      my.lp <- make.lp(0, nrVars)

      set.objfn(my.lp, obj)
      set.bounds(my.lp, upper = bounds$upper$val)
      set.bounds(my.lp, lower = bounds$lower$val)

      for ( i in 1:g_nr_rows(con) ) {
        r <- g_row(con, input=list(i))
        cols <- g_col_ind(r)
        vals <- g_values(r)
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

