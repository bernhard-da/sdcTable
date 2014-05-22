## Wrapper function for pasting key-Variables
pasteStrVec <- function(strVec, nrVars, coll=NULL) {
  if(length(strVec) %% nrVars != 0)
    stop("Wrong Dimensions!\n")

  if ( is.null(coll) ) {
    .Call( "myPaste", as.character(strVec), nrVars, PACKAGE = "sdcTable")
  } else {
    .Call( "myPasteWithSep", as.character(strVec), nrVars, coll, PACKAGE = "sdcTable")
  }
}

# alternative to expand.grid (used for pasteStrVec!)
expand <- function(inputList, vector=TRUE) {
  uniques <- sapply(inputList, length)
  nrPoss <- prod(uniques)
  if ( vector == TRUE ) {
    out <- NULL
    for ( i in 1:length(inputList) ) {
      if ( i == 1 )
        out <- rep(inputList[[i]], nrPoss/length(inputList[[i]]))
      else
        out <- c(out, rep(inputList[[i]], each=prod(uniques[1:(i-1)]), nrPoss/length(rep(inputList[[i]], each=prod(uniques[1:(i-1)])))))
    }
  }
  else {
    out <- list()
    for ( i in 1:length(inputList) ) {
      if ( i == 1 )
        out[[i]] <- rep(inputList[[i]], nrPoss/length(inputList[[i]]))
      else
        out[[i]] <- rep(inputList[[i]], each=prod(uniques[1:(i-1)]), nrPoss/length(rep(inputList[[i]], each=prod(uniques[1:(i-1)]))))
    }
  }
  out
}

# returns a vector original size or str
mySplit <- function(strVec, keepIndices) {
  if ( min(keepIndices) < 1 | max(keepIndices) > nchar(strVec[1]) ) {
    stop("indices must be in 1:",nchar(strVec[1]),"!\n")
  }
  keepIndices <- unique(keepIndices)-1 # required because of indexing in c++
  return(.Call( "mySplitFn", as.character(strVec), as.numeric(keepIndices), PACKAGE = "sdcTable"))
}

#strs <- rep(paste(LETTERS[1:6],collapse=""), 10000)
#system.time({
#	sapply(strs, mySplit, c(1,6))
#})

mySplitIndicesList <- function(strVec, keepList, coll="-") {
  u <- unlist(keepList)
  if ( min(u) < 1 | max(u) > nchar(strVec[1]) ) {
    stop("indices must be in 1:",nchar(strVec[1]),"!\n")
  }
  out <- list()
  for ( i in 1:length(keepList) ) {
    out[[i]] <- mySplit(strVec, keepList[[i]])
  }
  out <- .Call( "myPasteWithSep", as.character(unlist(out)), length(out), coll, PACKAGE = "sdcTable")
}
# mySplitIndicesList("112233444", list(1:3, 5:6, 7:8))

# check ob 'x' ganzzahlig ist
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  {
  abs(x - round(x)) < tol
}

is.zero <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - 0) < tol
}

is.one <- function(x, tol = .Machine$double.eps^0.5) {
  abs(x - 1) < tol
}

# welche Variable soll als Branching_Variable verwendet werden?
getBranchingVariable <- function(sol, alreadyBranched, primSupps) {
  ind <- setdiff(1:length(sol), c(alreadyBranched, primSupps))
  branchVar <- ind[which.min(0.5 - sol[ind])]
  branchVar
}

my.Rglpk_solve_LP <- function(obj, mat, dir, rhs, types = NULL, max = FALSE, bounds = NULL, verbose = FALSE) {
  if (!identical(max, TRUE) && !identical(max, FALSE))
    stop("'Argument 'max' must be either TRUE or FALSE.")
  direction_of_optimization <- as.integer(max)
  if (!identical(verbose, TRUE) && !identical(verbose, FALSE))
    stop("'Argument 'verbose' must be either TRUE or FALSE.")
  if ( !class(mat) == "simpleTriplet" ) {
    stop("mat must be of class 'simpleTriplet'")
  }

  verb <- as.integer(verbose)
  n_of_constraints <- length(dir)
  direction_of_constraints <- match(dir, c("<", "<=", ">", ">=", "=="))
  if (any(is.na(direction_of_constraints)))
    stop("Argument 'dir' must be either '<', '<=', '>', '>=' or '=='.")
  n_of_objective_vars <- length(obj)
  if (is.null(types))
    types <- "C"
  if (any(is.na(match(types, c("I", "B", "C"), nomatch = NA))))
    stop("'types' must be either 'B', 'C' or 'I'.")
  types <- rep(types, length.out = n_of_objective_vars)
  integers <- types == "I"
  binaries <- types == "B"
  is_integer <- any(binaries | integers)
  bounds <- as.glp_bounds(as.list(bounds), n_of_objective_vars)
  x <- glp_call_interface(
    obj,
    n_of_objective_vars,
    get.simpleTriplet(mat, type='rowInd', input=list()),
    get.simpleTriplet(mat, type='colInd', input=list()),
    get.simpleTriplet(mat, type='values', input=list()),
    length(get.simpleTriplet(mat, type='values', input=list())),

    rhs, direction_of_constraints, n_of_constraints, is_integer,
    integers, binaries,
    direction_of_optimization,
    bounds[,1L],
    bounds[,2L],
    bounds[,3L],
    verb
  )
  solution <- x$lp_objective_vars_values
  solution[integers | binaries] <- round(solution[integers | binaries])
  status <- as.integer(x$lp_status != 5L)
  list(optimum = sum(solution * obj), solution = solution, status = status)
}
environment(my.Rglpk_solve_LP) <- environment(Rglpk_solve_LP)

# calculates years, weeks, days, hours, minutes and seconds from integer number
# secs: count of elapsed seconds (proc.time()[3])
# returns also a formatted string
formatTime <- function(secs){
  time.vec <- rep(NA, 6)
  names(time.vec) <- c('seconds', 'minutes','hours', 'days','weeks','years')

  secs <- ceiling(secs)

  time.vec['years'] <- floor(secs / (60*60*24*7*52))
  if ( time.vec['years'] > 0 ) {
    secs <- secs - (time.vec['years'] * (60*60*24*7*52))
  }

  time.vec['weeks'] <- floor(secs / (60*60*24*7))
  if ( time.vec['weeks'] > 0 ) {
    secs <- secs - (time.vec['weeks'] * (60*60*24*7))
  }

  time.vec['days'] <- floor(secs / (60*60*24))
  if ( time.vec['days'] > 0 ) {
    secs <- secs - time.vec['days']*(60*60*24)
  }

  time.vec['hours'] <- floor(secs / (60*60))
  if ( time.vec['hours'] > 0 ) {
    secs <- secs - time.vec['hours']*(60*60)
  }

  time.vec['minutes'] <- floor(secs / (60))
  if ( time.vec['minutes'] > 0 ) {
    secs <- secs - time.vec['minutes']*(60)
  }

  time.vec['seconds'] <- secs
  time.vec <- rev(time.vec)

  # time str #
  x <- time.vec[time.vec!=0]
  shortNames <- sapply(1:length(x), function(y) { substr(names(x)[y], 1, nchar(names(x)[y])-1)  } )

  time.str <- NULL
  for ( i in seq_along(names(x))) {

    if ( length(x) == 1 ) {
      if ( x[i] > 1 ) {
        time.str <- paste(time.str, x[i], " ", names(x[i]), sep="")
      } else {
        time.str <- paste(time.str, x[i], " ", shortNames[i], sep="")
      }
    }
    else {

      if ( names(x)[i]=="seconds") {
        if ( x[i] > 1 ) {
          time.str <- paste(time.str, "and", x[i], names(x[i]), sep=" ")
        } else {
          time.str <- paste(time.str, "and", x[i], shortNames[i], sep=" ")
        }

      } else {
        if ( x[i] > 1 ) {
          time.str <- paste(time.str, x[i], " ", names(x[i]), sep="")
        } else {
          time.str <- paste(time.str, x[i], " ", shortNames[i], sep="")
        }

        if ( i != length(x)-1 ) {
          time.str <- paste(time.str,", ", sep="")
        }
      }
    }
  }
  return(list(time.vec=time.vec, time.str=time.str))
}

# create default parameter objects suitable for primary|secondary suppression
# if selection == 'control.primary': set arguments suitable for primary suppression
# if selection == 'control.secondary': set arguments suitable for secondary suppression
genParaObj <- function(selection, ...) {
  controlPrimary <- function(...) {
    ### setDefaults ###
    paraObj <- list()

    # freq.rule
    paraObj$maxN <- 3
    paraObj$allowZeros <- FALSE

    # p-percent rule
    paraObj$p <- 80

    # n,k rule
    paraObj$n <- 2
    paraObj$k <- 85

    # pq-rule
    paraObj$pq <- c(25, 50)

    paraObj$numVarInd <- NA

    newPara <- list(...)

    indexNumVarIndices <- which(names(newPara) == "numVarIndices")
    if ( length(indexNumVarIndices) == 0 ) {
      stop("genPara (type=='control.primary'): parameter 'numVarIndices' must be specified\n")
    } else {
      numVarIndices <- newPara[[indexNumVarIndices]]
    }

    for ( i in seq_along(newPara) ) {
      m <- match(names(newPara)[i], names(paraObj))
      if ( !is.na(m) ) {
        paraObj[[m]] <- newPara[[i]]
      }
    }

    #if ( any(sapply(paraObj, length)!=1) ) {
    #	stop("genPara (type=='control.primary'): arguments for primary suppression are not valid!\n")
    #}
    if ( !is.logical(paraObj$allowZeros) ) {
      stop("genPara (type=='control.primary'): argument 'allowZeros' must be logical!\n")
    }
    if ( !all(c(is.numeric(paraObj$maxN), is.numeric(paraObj$p), is.numeric(paraObj$n), is.numeric(paraObj$k))) ) {
      stop("genPara (type=='control.primary'): arguments 'maxN', 'p', 'n' and 'k' must be numeric!\n")
    }
    if ( length(paraObj$pq) != 2 ) {
      stop("genPara (type=='control.primary'): length of argument 'pq' must equal 2!\n")
    }
    if ( paraObj$k < 1 | paraObj$k >= 100) {
      stop("genPara (type=='control.primary'): argument 'k' must be >= 1 and < 100!\n")
    }
    if ( paraObj$p < 1 | paraObj$p >= 100) {
      stop("genPara (type=='control.primary'): argument p must be >= 1 and < 100!\n")
    }
    if ( paraObj$pq[1] < 1 | paraObj$pq[1] >= 100) {
      stop("genPara (type=='control.primary'): argument 'p' of 'pq' must be >= 1 and < 100!\n")
    }
    if ( paraObj$pq[2] < 1 | paraObj$pq[2] >= 100) {
      stop("genPara (type=='control.primary'): argument 'q' of 'pq' must be >= 1 and < 100!\n")
    }
    if ( paraObj$pq[1] >= paraObj$pq[2] ) {
      stop("genPara (type=='control.primary'): argument 'p' of 'pq' must be < argument 'q' of 'pq'\n")
    }
    if ( !is.na(paraObj$numVarInd) ) {
      if ( !paraObj$numVarInd %in% 1:length(numVarIndices) ) {
        stop("genPara (type=='control.primary'): argument 'numVarInd' must be >= 1 and <=",length(numVarIndices),"!\n")
      }
    }
    return(paraObj)
  }

  ### create a parameter list with (...) changing the default-values -> used in protectTable()
  controlSecondary <- function(...) {
    ### setDefaults ###
    paraObj <- list()

    # general parameter
    paraObj$method <- NA
    paraObj$verbose <- FALSE
    paraObj$save <- FALSE
    paraObj$solver <- "glpk"

    # HITAS|OPT - parameter
    paraObj$maxIter <- 10
    paraObj$timeLimit <- NULL
    paraObj$maxVars <- NULL
    paraObj$fastSolution <- FALSE
    paraObj$fixVariables <- TRUE
    paraObj$approxPerc <- 10
    paraObj$useC <- FALSE

    # HYPERCUBE - parameter
    paraObj$protectionLevel <- 80
    paraObj$suppMethod <- "minSupps"
    paraObj$suppAdditionalQuader <- FALSE

    # protectLinkedTables
    paraObj$maxIter <- 5

    newPara <- list(...)
    for ( i in seq_along(newPara) ) {
      m <- match(names(newPara)[i], names(paraObj))
      if ( !is.na(m) ) {
        paraObj[[m]] <- newPara[[i]]
      }
    }

    ### checks
    if ( any(sapply(paraObj, length)!=1) ) {
      stop("genPara (type=='control.secondary'): arguments controlObj for sdc-procedure are not valid!\n")
    }
    if ( !all(c(is.numeric(paraObj$maxIter), is.numeric(paraObj$approxPerc), is.numeric(paraObj$protectionLevel), is.numeric(paraObj$maxIter))) ) {
      stop("genPara (type=='control.secondary'): arguments 'maxIter', 'maxIter', 'protectionLevel' and 'maxIter' must be numeric!\n")
    }
    if ( !all(c(is.logical(paraObj$verbose), is.logical(paraObj$save), is.logical(paraObj$fastSolution), is.logical(paraObj$fixVariables), is.logical(paraObj$suppAdditionalQuader))) ) {
      stop("genPara (type=='control.secondary'): arguments 'verbose', 'save', 'fastSolution' 'fixVariables' and 'suppAdditionalQuader' must be numeric!\n")
    }
    if ( !is.null(paraObj$timeLimit) && !paraObj$timeLimit %in% 1:3000 ) {
      stop("genPara (type=='control.secondary'): argument 'timeLimit' must be >= 1 and <= 3000 minutes!\n")
    }
    if ( !length(paraObj$approxPerc) & !paraObj$approxPerc %in% 1:100 ) {
      stop("genPara (type=='control.secondary'): argument 'approxPerc' must be >= 1 and <= 100!\n")
    }
    if ( !paraObj$method %in% c('SIMPLEHEURISTIC', 'HITAS', 'HYPERCUBE', 'OPT') ) {
      stop("genPara (type=='control.secondary'): 'method' must be either 'SIMPLEHEURISTIC', 'HITAS', 'HYPERCUBE' or 'OPT'!\n")
    }
    if ( !paraObj$suppMethod %in% c('minSupps', 'minSum', 'minSumLogs') ) {
      stop("genPara (type=='control.secondary'): 'suppMethod' must be either 'minSupps', 'minSum' or 'minSumLogs'!\n")
    }
    return(paraObj)
  }

  if ( !selection %in% c('control.primary', 'control.secondary') ) {
    stop("genPara:: argument 'selection' must be either 'control.primary' or 'control.secondary'!\n")
  }

  if ( selection == 'control.primary' ) {
    paraObj <- controlPrimary(...)
  }
  if ( selection == 'control.secondary' ) {
    paraObj <- controlSecondary(...)
  }
  return(paraObj)
}

# convert simple triplet to matrix
st_to_mat <- function(x) {
  n.rows <- get.simpleTriplet(x, type='nrRows', input=list())
  n.cols <- get.simpleTriplet(x, type='nrCols', input=list())
  M <- matrix(0, nrow=n.rows, ncol=n.cols)

  i.x <- get.simpleTriplet(x, type='rowInd', input=list())
  j.x <- get.simpleTriplet(x, type='colInd', input=list())
  v.x <- get.simpleTriplet(x, type='values', input=list())
  for ( i in 1:get.simpleTriplet(x, type='nrCells', input=list()) ) {
    M[i.x[i], j.x[i]] <- v.x[i]
  }
  # matrizen from attackers problem are transposed -> switch!
  return(t(M))
}

csp_cpp <- function(sdcProblem, attackonly=FALSE, verbose) {
  pI <- g_problemInstance(sdcProblem)
  dimInfo <- g_dimInfo(sdcProblem)
  aProb <- c_make_att_prob(input=list(objectA=pI, objectB=dimInfo))$aProb

  # already suppressed cells
  ind_prim <- as.integer(sort(c(g_primSupps(pI), g_secondSupps(pI))))
  len_prim <- as.integer(length(ind_prim))
  bounds_min <- bounds_max <- rep(0, len_prim)

  ind_fixed <- as.integer(g_forcedCells(pI))
  len_fixed <- as.integer(length(ind_fixed))

  aProb <- c_make_att_prob(input=list(objectA=pI, objectB=dimInfo))$aProb
  attProbM <- init.simpleTriplet("simpleTriplet", input=list(mat=st_to_mat(aProb@constraints)))

  ia <- as.integer(c(0, get.simpleTriplet(attProbM, type="rowInd", list())))
  ja <- as.integer(c(0, get.simpleTriplet(attProbM, type="colInd", list())))
  ar <- as.double(c(0, get.simpleTriplet(attProbM, type="values", list())))

  cells_mat <- as.integer(length(ia))
  nr_vars <- as.integer(get.simpleTriplet(attProbM, type="nrCols", list()))
  nr_rows <- as.integer(get.simpleTriplet(attProbM, type="nrRows", list()))

  vals <- as.integer(g_freq(pI))

  lb <- as.double(g_lb(pI))
  ub <- as.double(g_ub(pI))

  LPL <- as.integer(g_LPL(pI))
  UPL <- as.integer(g_UPL(pI))
  SPL <- as.integer(g_SPL(pI))

  if ( attackonly == TRUE ) {
    attackonly <- as.integer(1)
  } else {
    attackonly <- as.integer(0)
  }
  final_pattern <- as.integer(rep(0, length(vals)))
  time.start <- proc.time()
  res <- .C("csp",
    ind_prim=ind_prim,
    len_prim=len_prim,
    bounds_min=bounds_min,
    bounds_max=bounds_max,
    ind_fixed=ind_fixed,
    len_fixed=len_fixed,
    ia=ia,
    ja=ja,
    ar=ar,
    cells_mat=cells_mat,
    nr_vars=nr_vars,
    nr_rows=nr_rows,
    vals=vals,
    lb=lb, ub=ub,
    LPL=LPL,
    UPL=UPL,
    SPL=SPL,
    final_pattern=final_pattern,
    attackonly=attackonly,
    verbose=as.integer(verbose),
    is_ok=0L
  )

  if ( attackonly ) {
    df <- data.frame(prim_supps=res$ind_prim, val=res$vals[res$ind_prim], bounds_low=res$bounds_min, bounds_up=res$bounds_max)
    df$protected <- df$bounds_low <= df$val - LPL[df$prim_supps]  &
    df$bounds_up >=  df$val + UPL[df$prim_supps] &
      df$bounds_up - df$bounds_low >= SPL[df$prim_supps]

    if ( length(g_secondSupps(pI)) > 0 ) {
      index <- g_primSupps(pI)
      df <- df[which(df$prim_supps %in% index),]
    }
    return(df)
  } else {
    if  ( res$is_ok != 0 ) {
      warning("no valid solution was obtained!\n")
      return(NULL)
    } else {
      nr_vars <- g_nrVars(g_problemInstance(sdcProblem))
      status_new <- rep("s", nr_vars)
      status_new[res$final_pattern!=0] <- "x"
      status_new[ind_prim] <- "u"
      if ( length(g_secondSupps(pI)) > 0 ) {
        status_new[g_secondSupps(pI)] <- "x"
      }
      if ( length(ind_fixed) > 0 ) {
        status_new[ind_fixed] <- "z"
      }
      
      pI <- g_problemInstance(sdcProblem)
      s_sdcStatus(pI) <- list(index=1:nr_vars, vals=status_new)
      s_problemInstance(sdcProblem) <- pI

      time.el <- g_elapsedTime(sdcProblem)+(proc.time()-time.start)[3]
      s_elapsedTime(sdcProblem) <- time.el
      s_indicesDealtWith(sdcProblem) <- 1:nr_vars
      return(sdcProblem)
    }
  }
}

### Primaersperrungen ###
# object (class=sdcProblem)
performQuickSuppression <- function(object, input) {
  suppMultDimTable <- function(dat, dimVars, freqInd) {
    # protect n-dimensional table
    simpleSupp <- function(splList, freqInd) {
      runInd <- TRUE
      counter <- 0
      override <- FALSE
      while(runInd) {
        runInd <- FALSE
        counter <- counter + 1
        #cat("run:", counter,"\n")
        for ( i in 1:length(splList)) {
          allOk <- all(splList[[i]]$sdcStatus %in% c("z"))

          if ( !allOk & nrow(spl[[i]]) > 1 & length(which(splList[[i]]$sdcStatus %in% c("u", "x")))==1 ) {
            #cat("we need to doe something: i=",i,"\n")
            runInd <- TRUE
            ind.x <- which(splList[[i]]$sdcStatus=='s')
            f <- splList[[i]][,freqInd]
            toSupp <- ind.x[order(f[ind.x], decreasing=FALSE)[1]]
            #cat("toSupp:", toSupp,"\n")
            if ( is.na(toSupp) ) {
              cat("Problem bei i=",i,"\n")
              ind.x <- which(splList[[i]]$sdcStatus %in% c('s','z') & splList[[i]]$freq!=0)
              f <- splList[[i]][,freqInd]
              toSupp <- ind.x[order(f[ind.x], decreasing=FALSE)[1]]
              override <- TRUE

              if ( splList[[i]]$freq[toSupp]==0) {
                stop("Fehler!\n")
              }
            }
            splList[[i]]$sdcStatus[toSupp] <- 'x'
          }
        }
      }

      s <- do.call("rbind", splList)
      rownames(s) <- NULL
      suppsAdded <- TRUE
      if ( counter == 1 ) {
        suppsAdded <- FALSE
      }
      return(list(s=s, suppsAdded=suppsAdded, override=override))
    }

    nDims <- length(dimVars)
    combs <- combn(nDims, nDims-1)

    runInd <- TRUE
    counter <- 0
    override <- FALSE
    patternOrig <- dat$sdcStatus
    while ( runInd ) {
      counter <- counter + 1
      suppsAdded <- rep(NA, ncol(combs))
      for ( i in 1:ncol(combs)) {
        f <- apply(dat, 1, function(x) { paste(x[combs[,i]], collapse="-") } )
        spl <- split(dat, f)
        res <- simpleSupp(spl, freqInd)

        dat <- res$s
        if ( override == FALSE & res$override == TRUE) {
          override <- TRUE
        }
        suppsAdded[i] <- res$suppsAdded
      }
      #cat("counter:", counter, "\n")
      #cat("suppsAdded:\n"); print(suppsAdded)
      if ( all(suppsAdded == FALSE) ) {
        #cat("finished! (counter=",counter,")\n")
        runInd <- FALSE
      }
    }
    pattern <- dat$sdcStatus
    pattern[which(dat$sdcStatus =="s")] <- "z"
    return(list(pattern=pattern, ids=dat$id, override=override))
  }

  verbose <- input$verbose
  pI <- g_problemInstance(object)

  strIDs <- g_strID(pI)

  dat <- data.frame(
    id=1:length(strIDs),
    strID=strIDs,
    freq=g_freq(pI),
    sdcStatus=g_sdcStatus(pI), stringsAsFactors=F
  )

  indices <- g_partition(object)$indices
  dimInfo <- g_dimInfo(object)
  strInfo <- get.dimInfo(dimInfo, type="strInfo")
  vNames <- get.dimInfo(dimInfo, type="varName")

  for ( i in seq_along(vNames) ) {
    dat[,vNames[i]] <- str_sub(dat$strID, strInfo[[i]][1], strInfo[[i]][2])
  }

  dimVars <- 1:length(vNames)
  dat <- cbind(dat[,5:ncol(dat)], dat[,1:4])
  #freqInd <- length(vNames)+3
  freqInd <- match("freq", colnames(dat))

  runInd <- TRUE
  while( runInd ) {
    override <- FALSE
    for ( i in 1:length(indices) ) {
      for ( j in 1:length(indices[[i]]) ) {
        curIndices <- indices[[i]][[j]]
        subDat <- dat[curIndices,]

        nrSupps <- length(which(subDat$sdcStatus%in%c("u","x")))

        if ( nrSupps > 0 ) {
          if ( verbose ) {
            cat("group:",i,"| ")
            cat("table",j,"/",length(indices[[i]]),"| ")
            cat("nrCells:",length(curIndices),"| ")
            cat("nrPrimSupps:",nrSupps,"| ")
            cat("override:",override,"\n")
          }

          res <- suppMultDimTable(subDat, dimVars, freqInd)
          matchInd <- match(curIndices, res$ids)
          dat$sdcStatus[curIndices] <- res$pattern[matchInd]
          if ( override == FALSE & res$override==TRUE ) {
            override <- TRUE
          }
        } else {
          ind <- which(subDat$sdcStatus=="s")
          dat$sdcStatus[curIndices[ind]] <- "z"
        }
      }
      if ( length(which(dat$freq==0 & dat$sdcStatus%in%c("u","x"))) > 0 ) {
        stop("fehler2!\n")
      }
    }
    if ( override == TRUE ) {
      # alle zellen %in% c("s", "z") muessen auf "s" gesetzt werden
      dat[dat$sdcStatus %in% c("s","z"),"sdcStatus"] <- "s"
      # new
      dat$sdcStatus[dat$freq==0] <- "z"
    } else {
      if ( verbose ) {
        cat("finished!\n")
      }
      runInd <- FALSE
    }
  }

  matchID <- match(dat$strID,strIDs )
  s_sdcStatus(pI) <- list(index=matchID, vals=dat$sdcStatus)
  s_problemInstance(object) <- pI
  object
}

# using optimally solving a problem using c++ implementation (used in protectTable())
opt_cpp <- function(object, input) {
  if ( !class(object) == "sdcProblem" ) {
    stop("check input 'object' (not of class 'sdcProblem')!\n")
  }

  timeLimit <- input$timeLimit
  verbose <- input$verbose
  save <- input$save

  start.time <- proc.time()
  pI <- g_problemInstance(object)
  sdcStatusBegin <- g_sdcStatus(pI)
  primSupps <- primSuppsOrig <- g_primSupps(pI)

  indexPool <- numeric()
  allStrIDs <- g_strID(pI)

  s_elapsedTime(object) <- g_elapsedTime(object) + (proc.time()-start.time)[3]
  invisible(csp_cpp(sdcProblem=object, attackonly=FALSE, verbose=input$verbose))
}

# heuristically solving a problem using c++ implementation (used in protectTable())
hitas_cpp <- function(object, input) {
  if ( !class(object) == "sdcProblem" ) {
    stop("check input 'object' (not of class 'sdcProblem')!\n")
  }

  timeLimit <- input$timeLimit
  verbose <- input$verbose
  save <- input$save

  start.time <- proc.time()
  pI <- g_problemInstance(object)
  sdcStatusBegin <- g_sdcStatus(pI)
  primSupps <- primSuppsOrig <- g_primSupps(pI)

  indexPool <- numeric()
  allStrIDs <- g_strID(pI)

  partition <- g_partition(object)
  run_ind <- TRUE
  while ( run_ind ) {
    for ( i in 1:(partition$nrGroups) ) {
      s_startJ(object) <- 1 # reset j before updating i
      s_startI(object) <- i

      indexPool <- NULL
      if ( i > 1 ) {
        indexPool <- sort(unique(unlist(partition$indices[1:(i-1)])))
      }
      ind <- partition$indices[[i]]

      for ( j in 1:(length(ind)) ) {
        is_ok <- TRUE
        s_startJ(object) <- j
        currentIndices <- ind[[j]] # within complete table

        ### cells with status 'u' or 'x' exist
        pI <- g_problemInstance(object)
        if ( any(g_sdcStatus(pI)[currentIndices] %in% c("u","x")) & length(currentIndices) > 1 ) {
          if ( verbose ) {
            cat("starting to solve problem",j,"/",length(ind),"in group",i,"/",partition$nrGroups,"!\n")
          }
          # if we have cells with "u" or "x" we need to protect
          # the corresponding subtable --> reduce problemInstance
          probNew <- c_reduce_problem(object, input=list(currentIndices))
          pI.new <- g_problemInstance(probNew)

          # is it necessary to protect the table??
          currentPrimSupps <- primSupps[!is.na(match(primSupps, currentIndices ))]

          # indices that have already been inner cells in tables dealt earlier
          # FIXME: save indexpool somehow
          indicesDealtWith <- which(currentIndices %in% indexPool) #in current current subproblem

          ### fix marginal-cells
          ### --> its-suppression state must not change!
          currentPattern <- g_sdcStatus(g_problemInstance(probNew))

          introducedSupps <- indicesDealtWith[which(currentPattern[indicesDealtWith] == "x")]
          if ( length(introducedSupps) > 0 ) {
            ### secondary suppressions from upper tables
            s_sdcStatus(pI.new) <- list(index=introducedSupps, vals=rep("u", length(introducedSupps)))

            ### temporarily change LPL, UPL, SPL for these cells
            LPL.current <- g_LPL(pI.new)[introducedSupps]
            UPL.current <- g_UPL(pI.new)[introducedSupps]
            SPL.current <- g_SPL(pI.new)[introducedSupps]

            s_LPL(pI.new) <- list(index=introducedSupps, vals=rep(0, length(introducedSupps)))
            s_UPL(pI.new) <- list(index=introducedSupps, vals=rep(0, length(introducedSupps)))
            s_SPL(pI.new) <- list(index=introducedSupps, vals=rep(0.1, length(introducedSupps)))            
          }

          ### force non-suppression of cells that have already been dealt with
          indForced <- indicesDealtWith[which(currentPattern[indicesDealtWith] == "s")]
          if ( length(indForced) > 0 ) {
            s_sdcStatus(pI.new) <- list(index=indForced, vals=rep("z", length(indForced)))
          }
          s_problemInstance(probNew) <- pI.new
          
          # solve the problem using c++ implementation
          res <- csp_cpp(sdcProblem=probNew, attackonly=FALSE, verbose=input$verbose)
          if ( is.null(res) ) {
            cat("\nWe got a problem and need to relax some conditions!\n\n")
            old.status <- probNew@problemInstance@sdcStatus
            ii <- which(old.status %in% c("z") & probNew@problemInstance@Freq > 0)
            if ( length(ii) == 0 ) {
              stop("This is a really nasty problem. No solution can be computed. Please contact the package maintainer.\n")
            }
            probNew@problemInstance@sdcStatus[ii] <- "s"
            res <- csp_cpp(sdcProblem=probNew, attackonly=FALSE, verbose=input$verbose)
            if ( is.null(res) ) {
              stop("This is a really nasty problem. No solution can be computed. Please contact the package maintainer.\n")
            }
            probNew <- res
            new.status <- probNew@problemInstance@sdcStatus
            xx <- data.frame(currentIndices=currentIndices,old=old.status, new=new.status, freq=probNew@problemInstance@Freq)

            ii <- which(new.status=="x" & old.status=="z")

            updated_status <- rep("s", length(old.status))
            updated_status[which(old.status=="u")] <- "u"
            updated_status[probNew@problemInstance@Freq==0] <- "z"
            updated_status[ii] <- "u"# previously "z", now "u"

            xx <- sdcStatusBegin
            xx[currentIndices] <- updated_status
            pI <- g_problemInstance(object)
            s_LPL(pI) <- list(index=currentIndices[ii], vals=rep(0, length(ii)))
            s_UPL(pI) <- list(index=currentIndices[ii], vals=rep(0, length(ii)))
            s_SPL(pI) <- list(index=currentIndices[ii], vals=rep(0.1, length(ii)))
            
            s_sdcStatus(pI) <- list(index=1:length(xx), vals=xx)
            s_problemInstance(object) <- pI
            is_ok <- FALSE
          } else {
            probNew <- res
          }
          # break j-loop
          if ( !is_ok ) {
            break
          }

          ### update sdcStatus
          status <- g_sdcStatus(g_problemInstance(probNew))

          pI <- g_problemInstance(object)
          if ( length(indForced) > 0 ) {
            status[indForced] <- "s"
          }
          if ( length(introducedSupps) > 0 ) {
            status[introducedSupps] <- "x"
            s_LPL(pI) <- list(index=currentIndices[introducedSupps], vals=LPL.current)
            s_UPL(pI) <- list(index=currentIndices[introducedSupps], vals=UPL.current)
            s_SPL(pI) <- list(index=currentIndices[introducedSupps], vals=SPL.current)        
          }
          s_sdcStatus(pI) <- list(index=currentIndices, vals=status)
          s_problemInstance(object) <- pI
        }
      }
      # break i-loop
      if ( !is_ok ) {
        break
      }
      ### update indices that we have already dealt with
      s_indicesDealtWith(object) <- unique(c(indexPool, currentIndices))
    }
    if ( is_ok ) {
      run_ind <- FALSE
    }
  }
  invisible(object)
}
