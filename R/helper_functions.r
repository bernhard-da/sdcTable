## Wrapper function for pasting key-Variables
pasteStrVec <- function(strVec, nrVars, coll=NULL) {
  if(length(strVec) %% nrVars != 0)
    stop("Wrong Dimensions!\n")
  if ( is.null(coll) ) {
    return(cpp_myPaste(as.character(strVec), as.integer(nrVars)[1], NA))
  } else {
    return(cpp_myPaste(as.character(strVec), as.integer(nrVars)[1], as.character(coll[1])))
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
  keepIndices <- unique(keepIndices)
  return(cpp_mySplit(as.character(strVec), as.integer(keepIndices)))
}

mySplitIndicesList <- function(strVec, keepList, coll="-") {
  u <- unlist(keepList)
  if ( min(u) < 1 | max(u) > nchar(strVec[1]) ) {
    stop("indices must be in 1:",nchar(strVec[1]),"!\n")
  }
  out <- list()
  for ( i in 1:length(keepList) ) {
    out[[i]] <- mySplit(strVec, keepList[[i]])
  }
  out <- cpp_myPaste(as.character(unlist(out)), as.integer(length(out)), coll)
}
# mySplitIndicesList("112233444", list(1:3, 5:6, 7:8))

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
    if (length(indexNumVarIndices) == 0) {
      stop("genPara (type=='control.primary'): parameter 'numVarIndices' must be specified\n")
    } else {
      numVarIndices <- newPara[[indexNumVarIndices]]
    }

    for (i in seq_along(newPara)) {
      m <- match(names(newPara)[i], names(paraObj))
      if (!is.na(m)) {
        paraObj[[m]] <- newPara[[i]]
      }
    }

    if (!is.logical(paraObj$allowZeros) ) {
      stop("genPara (type=='control.primary'): argument 'allowZeros' must be logical!\n")
    }
    if (!all(c(is.numeric(paraObj$maxN), is.numeric(paraObj$p), is.numeric(paraObj$n), is.numeric(paraObj$k))) ) {
      stop("genPara (type=='control.primary'): arguments 'maxN', 'p', 'n' and 'k' must be numeric!\n")
    }
    if (length(paraObj$pq) != 2) {
      stop("genPara (type=='control.primary'): length of argument 'pq' must equal 2!\n")
    }
    if (paraObj$k < 1 | paraObj$k >= 100) {
      stop("genPara (type=='control.primary'): argument 'k' must be >= 1 and < 100!\n")
    }
    if (paraObj$p < 1 | paraObj$p >= 100) {
      stop("genPara (type=='control.primary'): argument p must be >= 1 and < 100!\n")
    }
    if (paraObj$pq[1] < 1 | paraObj$pq[1] >= 100) {
      stop("genPara (type=='control.primary'): argument 'p' of 'pq' must be >= 1 and < 100!\n")
    }
    if (paraObj$pq[2] < 1 | paraObj$pq[2] >= 100) {
      stop("genPara (type=='control.primary'): argument 'q' of 'pq' must be >= 1 and < 100!\n")
    }
    if (paraObj$pq[1] >= paraObj$pq[2]) {
      stop("genPara (type=='control.primary'): argument 'p' of 'pq' must be < argument 'q' of 'pq'\n")
    }
    if (!is.na(paraObj$numVarInd)) {
      if (!paraObj$numVarInd %in% 1:length(numVarIndices) ) {
        stop("genPara (type=='control.primary'): argument 'numVarInd' must be >= 1 and <=",length(numVarIndices),"!\n")
      }
    }
    return(paraObj)
  }

  ### create a parameter list with (...) changing the default-values -> used in protectTable()
  controlSecondary <- function(...) {
    ### setDefaults ###
    paraObj <- list()

    paraObj$method <- "SIMPLEHEURISTIC"
    # general parameter
    paraObj$verbose <- FALSE
    paraObj$detectSingletons <- FALSE

    # protectLinkedTables
    paraObj$maxIter <- 10

    newPara <- list(...)
    for (i in seq_along(newPara)) {
      m <- match(names(newPara)[i], names(paraObj))
      if (!is.na(m)) {
        paraObj[[m]] <- newPara[[i]]
      }
    }

    ### checks
    if (any(sapply(paraObj, length)!=1)) {
      stop("genPara (type=='control.secondary'): arguments controlObj for sdc-procedure are not valid!\n")
    }
    if (!is.numeric(paraObj$maxIter)) {
      stop("genPara (type=='control.secondary'): argument 'maxIter'  must be numeric!\n")
    }
    if (!all(c(is.logical(paraObj$verbose), is.logical(paraObj$detectSingletons)))) {
      stop("genPara (type=='control.secondary'): arguments 'verbose' and 'detectSingletons' must be logical!\n")
    }
    if (!is.null(paraObj$timeLimit) && !paraObj$timeLimit %in% 1:3000) {
      stop("genPara (type=='control.secondary'): argument 'timeLimit' must be >= 1 and <= 3000 minutes!\n")
    }
    return(paraObj)
  }

  if (!selection %in% c('control.primary', 'control.secondary')) {
    stop("genPara:: argument 'selection' must be either 'control.primary' or 'control.secondary'!\n")
  }

  if (selection == "control.primary") {
    paraObj <- controlPrimary(...)
  }
  if (selection == "control.secondary") {
    paraObj <- controlSecondary(...)
  }
  return(paraObj)
}

singletonDetectionProcedure <- function(dat, indices, subIndices) {
  id <- freq <- sdcStatus <- NULL
  nrAddSupps <- 0
  suppIds <- c()

  # temporarily recode primary suppressions and check, if they are really singletons
  id_changed <- dat[sdcStatus=="u" & freq>1, id]
  if (length(id_changed)>0) {
    dat[id_changed, sdcStatus:="x"]
  }
  for (i in 1:length(indices)) {
    sI <- subIndices[[i]]
    for (j in 1:length(sI)) {
      sJ <- sI[[j]]
      for (z in 1:length(sJ)) {
        poss <- sJ[[z]]
        mm <- max(poss)
        for (k in 1:mm) {
          ii <- indices[[i]][[j]][which(poss==k)]
          # only if we have a real subtable
          if (length(ii) > 1) {
            # tau-argus strategy
            ind_u <- which(dat$sdcStatus[ii]=="u")

            if (length(ind_u)==2) {
              ff <- dat$freq[ii[ind_u]]
              # at least one cell of two supps is a singleton
              # 1. If on a row or column of a subtable there are only two singletons and no other
              # primary suppressions.
              # 2. If there is only one singleton and one multiple primary unsafe cell.
              # one or two singletons
              if (any(ff==1) & sum(dat$sdcStatus[ii]=="x")==0 & sum(dat$freq[ii]>0)>2) {
                # we have two singletons, we need to add one additional suppression
                ss <- dat[ii]
                ss <- ss[sdcStatus=="s"]
                suppId <- ss$id[which.min(ss$freq)]
                if (length(suppId)==0) {
                  stop("error finding an additional primary suppression (1)\n")
                }
                if (dat[suppId, freq]==1) {
                  dat[suppId, sdcStatus:="u"]
                } else {
                  dat[suppId, sdcStatus:="x"]
                }
                nrAddSupps <- nrAddSupps + 1
                suppIds <- c(suppIds, suppId)
              }
            }
            # 3. If a frequency rule is used, it could happen that two cells on a row/column are
            # primary unsafe, but the sum of the two cells could still be unsafe. In that case
            # it should be prevented that these two cells protect each other.
            if (length(ind_u)==3) {
              # the sum is primary suppressed, thus the other two primary suppressions are within the row/col
              if (dat$sdcStatus[ii[1]]=="u" & sum(dat$freq[ii]>0)>3) {
                # we need to find an additional suppression
                ss <- dat[ii]
                ss <- ss[sdcStatus=="s"]
                suppId <- ss$id[which.min(ss$freq)]
                if (length(suppId)==0) {
                  stop("error finding an additional primary suppression (2)\n")
                }
                if (dat[suppId, freq]==1) {
                  dat[suppId, sdcStatus:="u"]
                } else {
                  dat[suppId, sdcStatus:="x"]
                }
                nrAddSupps <- nrAddSupps + 1
                suppIds <- c(suppIds, suppId)
              }
            }
          }
        }
      }
    }
  }

  # reset primary suppressions
  if (length(id_changed)>0) {
    dat[id_changed, sdcStatus:="u"]
  }
  #if (length(nrAddSupps)>0) {
  #  dat[suppIds, sdcStatus:="u"]
  #}
  invisible(list(dat=dat, nrAddSupps=nrAddSupps, suppIds=suppIds))
}
