#' @aliases get.sdcProblem,sdcProblem,character-method
#' @rdname get.sdcProblem-method
setMethod(f="get.sdcProblem", signature=c("sdcProblem", "character"),
  definition=function(object, type) {
    if ( !type %in% c("dataObj", "problemInstance", "partition", "elapsedTime", "dimInfo", "indicesDealtWith",
        "startI", "startJ") ) {
      stop("get.sdcProblem:: argument 'type' is not valid!\n")
    }
    if ( type == "dataObj" ) {
      return(g_dataObj(object))
    }
    if ( type == "problemInstance" ) {
      return(g_problemInstance(object))
    }
    if ( type == "partition" ) {
      return(g_partition(object))
    }
    if ( type == "elapsedTime" ) {
      return(g_elapsedTime(object))
    }
    if ( type == "dimInfo" ) {
      return(g_dimInfo(object))
    }
    if ( type == "indicesDealtWith" ) {
      return(g_indicesDealtWith(object))
    }
    if ( type == "startI" ) {
      return(g_startI(object))
    }
    if ( type == "startJ" ) {
      return(g_startJ(object))
    }
  }
)

#' @aliases set.sdcProblem,sdcProblem,character,list-method
#' @rdname set.sdcProblem-method
setMethod(f="set.sdcProblem", signature=c("sdcProblem", "character", "list"),
  definition=function(object, type, input) {
    if (!type %in% c("problemInstance", "partition", "rule.freq", "rule.nk",
        "rule.p", "rule.pk", "startI", "startJ", "indicesDealtWith", "elapsedTime") ) {
      stop("set.sdcProblem:: check argument 'type'!\n")
    }
    if (type == "problemInstance") {
      s_problemInstance(object) <- input[[1]]
    }
    if (type == "partition") {
      s_partition(object) <- input[[1]]
    }
    if (type == "startI") {
      s_startI(object) <- input[[1]]
    }
    if (type == "startJ") {
      s_startJ(object) <- input[[1]]
    }
    if (type == "indicesDealtWith") {
      s_indicesDealtWith(object) <- input[[1]]
    }
    if (type == "elapsedTime") {
      s_elapsedTime(object) <- input[[1]]
    }
    validObject(object)
    return(object)
  }
)

#' @aliases calc.sdcProblem,sdcProblem,character,list-method
#' @rdname calc.sdcProblem-method
setMethod(f="calc.sdcProblem", signature=c("sdcProblem", "character", "list"),
  definition=function(object, type, input) {
    if (!type %in% c("rule.freq", "rule.nk", "rule.p", "rule.pq",
      "cellID", "finalize", "contributingIndices")) {
      stop("calc.sdcProblem:: check argument 'type'!\n")
    }
    # frequency-rule
    if (type == "rule.freq" ) {
      return(c_rule_freq(object, input))
    }
    # nk-dominance rule
    if (type == "rule.nk") {
      return(c_rule_nk(object, input))
    }
    # p%-rule
    if (type == "rule.p") {
      return(c_rule_p(object, input))
    }
    # pq-rule
    if (type == "rule.pq") {
      return(c_rule_pq(object, input))
    }
    if (type == "cellID") {
      return(c_cellID(object, input))
    }
    if (type == "finalize") {
      return(c_finalize(object, input))
    }
    if (type == "contributingIndices") {
      return(c_contributing_indices(object, input))
    }
  }
)

#' summarize object of class \code{\link{sdcProblem-class}} or \code{\link{safeObj-class}}.
#'
#' extract and show relevant information stored in object ofs class \code{\link{sdcProblem-class}} or \code{\link{safeObj-class}}.
#'
#' @aliases summary,sdcProblem-method
#' @rdname summary.sdcProblem-method
#' @param object Objects of either class \code{\link{sdcProblem-class}} or \code{\link{safeObj-class}}.
#' @param ... currently not used.
#' @export
#' @docType methods
setMethod(f="summary", signature="sdcProblem",
  definition=function(object, ...) {
    pI <- g_problemInstance(object)
    dO <- g_dataObj(object)
    dI <- g_dimInfo(object)
    if ( g_is_microdata(dO) ) {
      cat("The raw data contains micro data!")
      if ( length(pI@numVars) > 0 ) {
        cat("--> the use of dominance rules for primary suppressions is possible!")
      }
      cat("\n")
    } else {
      cat("The raw data contain pre-aggregated (tabular) data!\n")
    }

    nrcells <- g_nrVars(pI)
    dim_names <- g_varname(dI)

    cat("\nThe complete table to protect consists of",nrcells,"cells and has",length(dim_names),"spanning variables.")

    cat("\nThe distribution of\n")
    cat("- primary unsafe (u)\n")
    cat("- secondary suppressed (x)\n")
    cat("- forced to publish (z) and\n")
    cat("- selectable for secondary suppression (s) cells is shown below:\n")
    print(table(g_sdcStatus(pI)))

    nr_tables <- g_partition(object)$nrTables
    cat("\nIf this table is protected with heuristic methods, a total of",nr_tables,"has (sub)tables must be considered!\n")
  }
)


#' print objects of class \code{\link{sdcProblem-class}}.
#'
#' print some useful information instead of just displaying the entire object (which may be large)
#'
#' @aliases print,sdcProblem-method
#' @rdname print.sdcProblem-method
#' @param x an objects of class \code{\link{sdcProblem-class}}
#' @param ... currently not used.
#' @export
#' @docType methods
setMethod("print", signature="sdcProblem",
  definition=function(x, ...) {
    dims <- x@dimInfo@dimInfo
    nrDims <- length(dims)
    nrCells <- length(x@problemInstance@strID)
    cat(paste("The object is an 'sdcProblem' with",nrCells,"cells in",nrDims, "dimension(s)!\n"))
    cat("\nThe dimensions are:\n")
    for ( i in 1:nrDims ) {
      nrCodes <- length(dims[[i]]@codesOriginal)
      nrAggregates <- sum(dims[[i]]@codesMinimal==FALSE)
      maxHier <- length(dims[[i]]@structure)
      cat(paste0("\t- ",names(dims)[i]," (",maxHier," levels; ",nrCodes," codes; of these being ",nrAggregates," aggregates)\n"))
    }
    cat("\nCurrent suppression pattern:\n")
    cat("\t- Primary suppressions:",sum(x@problemInstance@sdcStatus=="u"),"\n")
    cat("\t- Secondary suppressions:",sum(x@problemInstance@sdcStatus=="x"),"\n")
    cat("\t- Publishable cells:",sum(x@problemInstance@sdcStatus%in% c("s","z")),"\n")
  }
)

#' show objects of class \code{\link{sdcProblem-class}}.
#'
#' just calls the corresponding print-method
#'
#' @aliases show,sdcProblem-method
#' @rdname show.sdcProblem-method
#' @param object an objects of class \code{\link{sdcProblem-class}}
#' @export
#' @docType methods
setMethod("show", signature="sdcProblem",
  definition=function(object) {
    print(object)
  }
)

setMethod("g_problemInstance", signature="sdcProblem", definition=function(object) {
  object@problemInstance
})

setMethod("g_dimInfo", signature="sdcProblem", definition=function(object) {
  object@dimInfo
})

setMethod("g_partition", signature="sdcProblem", definition=function(object) {
  object@partition
})

setMethod("g_elapsedTime", signature="sdcProblem", definition=function(object) {
  object@elapsedTime
})

setMethod("g_dataObj", signature="sdcProblem", definition=function(object) {
  object@dataObj
})

setMethod("g_startI", signature="sdcProblem", definition=function(object) {
  object@startI
})

setMethod("g_startJ", signature="sdcProblem", definition=function(object) {
  object@startJ
})

setMethod("g_indicesDealtWith", signature="sdcProblem", definition=function(object) {
  object@indicesDealtWith
})

setMethod("g_df", signature="sdcProblem", definition=function(object, addDups=FALSE, addNumVars=FALSE) {
  xx <- strID <- NULL
  pI <- g_problemInstance(object)
  dt <- data.table(
    strID=g_strID(pI),
    freq=g_freq(pI),
    sdcStatus=g_sdcStatus(pI)
  )
  if (addNumVars & !is.null(pI@numVars) ) {
    dt <- cbind(dt, as.data.table(pI@numVars))
  }
  dI <- g_dimInfo(object)
  strInfo <- g_str_info(dI)
  dimObj <- g_dim_info(dI)
  vNames <- g_varname(dI)
  res <- as.data.table(cpp_splitByIndices(g_strID(pI), strInfo))
  setnames(res, vNames)
  dt <- cbind(dt, res)
  for ( i in 1:length(strInfo) ) {
    v <- paste0(vNames[i],"_o",sep="")
    dt[[v]] <- c_match_orig_codes(object=dimObj[[i]], input=dt[[vNames[i]]])
  }

  if ( addDups ) {
    dims <- g_dim_info(dI)
    for ( i in seq_along(dims) ) {
      if ( g_has_dups(dims[[i]]) ) {
        dU <- dims[[i]]@dupsUp
        dL <- dims[[i]]@dups
        vName <- paste0(dims[[i]]@vName,"_o")
        for ( j in 1:length(dL) ) {
          cmd <- paste0("xx <- dt[",vName,"=='",dU[j],"']")
          eval(parse(text=cmd))
          if ( !is.numeric(dt[[vName]]) ) {
            cmd <- paste0("xx[,",vName,":='",dL[j],"']")
          } else {
            cmd <- paste0("xx[,",vName,":=",dL[j],"]")
          }
          eval(parse(text=cmd))
          dt <- rbind(dt, xx); rm(xx)
        }
      }
    }
  }
  setkey(dt, strID)
  return(dt)
})

setReplaceMethod("s_problemInstance", signature=c("sdcProblem", "problemInstance"), definition=function(object, value) {
  object@problemInstance <- value
  validObject(object)
  object
})

setReplaceMethod("s_partition", signature=c("sdcProblem"), definition=function(object, value) {
  object@partition <- value
  validObject(object)
  object
})

setReplaceMethod("s_startI", signature=c("sdcProblem", "numeric"), definition=function(object, value) {
  object@startI <- value
  validObject(object)
  object
})

setReplaceMethod("s_startJ", signature=c("sdcProblem", "numeric"), definition=function(object, value) {
  object@startJ <- value
  validObject(object)
  object
})

setReplaceMethod("s_indicesDealtWith", signature=c("sdcProblem"), definition=function(object, value) {
  object@indicesDealtWith <- value
  validObject(object)
  object
})

setReplaceMethod("s_elapsedTime", signature=c("sdcProblem"), definition=function(object, value) {
  object@elapsedTime <- value
  validObject(object)
  object
})

setMethod("c_rule_freq", signature=c("sdcProblem", "list"), definition=function(object, input) {
  pI <- g_problemInstance(object)
  if ( input$allowZeros == TRUE ) {
    suppInd <- which(g_freq(pI) <= input$maxN)
  } else {
    f <- g_freq(pI)
    suppInd <- which(f > 0 & f <= input$maxN)
    zeroInd <- which(g_freq(pI) == 0 )
    if ( length(zeroInd) > 0 ) {
      s_sdcStatus(pI) <- list(index=zeroInd, vals=rep("z", length(zeroInd)))
    }
  }
  if ( length(suppInd) > 0 ) {
    s_sdcStatus(pI) <- list(index=suppInd, vals=rep("u", length(suppInd)))
  }
  s_problemInstance(object) <- pI
  validObject(object)
  return(object)
})

setMethod("c_rule_nk", signature=c("sdcProblem", "list"), definition=function(object, input) {
  nkRule <- function(celltot, sumNcont, k) {
    # if TRUE, cell needs to be suppressed
    (sumNcont) > (k/100*celltot)
  }
  if ( !g_is_microdata(g_dataObj(object)) ) {
    stop("nk-dominance rule can only be applied if micro-data are available!\n")
  }
  pI <- g_problemInstance(object)
  dataObj <- g_dataObj(object)
  strIDs <- g_strID(pI)
  numVarInds <- g_numvar_ind(dataObj)

  numVal <- g_raw_data(dataObj)[[numVarInds[input$numVarInd]]]
  if ( any(numVal < 0 ) ) {
    stop("dominance rules can only be applied to numeric variables with only positive values!\n")
  }

  # calculate contributing indices
  indices <- lapply(1:g_nrVars(pI), function(x) {
    c_contributing_indices(object, input=list(strIDs[x]))
  })

  minContributingUnits <- min(setdiff(unique(sapply(indices, length)), 0))
  if ( input$n < 1 | input$n > minContributingUnits ) {
    stop("set.sdcProblem:: parameter 'n' must be >= 1 and <",minContributingUnits,"!\n")
  }

  # values of contributing units
  valueList <- lapply(1:g_nrVars(pI), function(x) {
    sum(rev(tail(sort(numVal[indices[[x]]]), input$n)))
  })

  cellTotals <- g_numVars(pI)[[input$numVarInd]]
  # suppStatus: TRUE:unsafe, FALSE: safe
  nkState <- sapply(1:g_nrVars(pI), function(x) {
    nkRule(cellTotals[x], valueList[[x]], input$k)
  })

  addSupps <- which(sapply(indices, length) %in% 1:input$n)
  if ( length(addSupps) > 0 ) {
    nkState[addSupps] <- TRUE
  }

  suppIndex <- which(nkState==TRUE)
  if ( length(suppIndex) > 0 ) {
    s_sdcStatus(pI) <- list(index=suppIndex, vals=rep("u", length(suppIndex)))
  }

  if ( input$allowZeros == FALSE ) {
    indZero <- which(g_freq(pI)==0)
    if ( length(indZero) > 0 ) {
      s_sdcStatus(pI) <- list(index=indZero, vals=rep("z", length(indZero)))
    }
  }
  s_problemInstance(object) <- pI
  validObject(object)
  return(object)
})

setMethod("c_rule_p", signature=c("sdcProblem", "list"), definition=function(object, input) {
  pPercRule <- function(celltot, cont1, cont2, p) {
    # if TRUE, cell needs to be suppressed
    (celltot - cont1 - cont2) < (p/100*cont1)
  }

  if ( !g_is_microdata(g_dataObj(object)) ) {
    stop("p-percent rule can only be applied if micro-data are available!\n")
  }

  pI <- g_problemInstance(object)
  dataObj <- g_dataObj(object)
  numVarInds <- g_numvar_ind(dataObj)
  strIDs <- g_strID(pI)

  numVal <- g_raw_data(dataObj)[[numVarInds[input$numVarInd]]]
  if ( any(numVal < 0 ) ) {
    stop("dominance rules can only be applied to numeric variables with only positive values!\n")
  }

  # calculate contributing indices
  indices <- lapply(1:g_nrVars(pI), function(x) {
    c_contributing_indices(object, input=list(strIDs[x]))
  })

  # values of contributing units
  valueList <- lapply(1:g_nrVars(pI), function(x) {
    rev(tail(sort(numVal[indices[[x]]]),2))
  })
  cellTotals <- g_numVars(pI)[[input$numVarInd]]

  # suppStatus: TRUE:unsafe, FALSE: safe
  pState <- sapply(1:g_nrVars(pI), function(x) {
    pPercRule(cellTotals[x], valueList[[x]][1], valueList[[x]][2], input$p)
  })

  suppIndex <- which(pState==TRUE)
  if ( length(suppIndex) > 0 ) {
    s_sdcStatus(pI) <- list(index=suppIndex, vals=rep("u", length(suppIndex)))
  }

  if ( input$allowZeros == FALSE ) {
    indZero <- which(g_freq(pI)==0)
    if ( length(indZero) > 0 ) {
      s_sdcStatus(pI) <- list(index=indZero, vals=rep("u", length(indZero)))
    }
  }
  s_problemInstance(object) <- pI
  validObject(object)
  return(object)
})

setMethod("c_rule_pq", signature=c("sdcProblem", "list"), definition=function(object, input) {
  pqRule <- function(celltot, cont1, cont2, p, q) {
    # if TRUE, cell needs to be suppressed
    (celltot - cont1 - cont2) < (p/q)*cont1
  }
  if ( !g_is_microdata(g_dataObj(object)) ) {
    stop("p-percent rule can only be applied if micro-data are available!\n")
  }

  pI <- g_problemInstance(object)
  dataObj <- g_dataObj(object)
  numVarInds <- g_numvar_ind(dataObj)
  strIDs <- g_strID(pI)

  numVal <- g_raw_data(dataObj)[[numVarInds[input$numVarInd]]]
  if ( any(numVal < 0 ) ) {
    stop("dominance rules can only be applied to numeric variables with only positive values!\n")
  }

  # calculate contributing indices
  indices <- lapply(1:g_nrVars(pI), function(x) {
    c_contributing_indices(object, input=list(strIDs[x]))
  })

  # values of contributing units
  valueList <- lapply(1:g_nrVars(pI), function(x) {
    rev(tail(sort(numVal[indices[[x]]]),2))
  })
  cellTotals <- g_numVars(pI)[[input$numVarInd]]

  # suppStatus: TRUE:unsafe, FALSE: safe
  pState <- sapply(1:g_nrVars(pI), function(x) {
    pqRule(cellTotals[x], valueList[[x]][1], valueList[[x]][2], input$pq[1], input$pq[2])
  })

  suppIndex <- which(pState==TRUE)
  if ( length(suppIndex) > 0 ) {
    s_sdcStatus(pI) <- list(index=suppIndex, vals=rep("u", length(suppIndex)))
  }

  if ( input$allowZeros == FALSE ) {
    indZero <- which(g_freq(pI)==0)
    if ( length(indZero) > 0 ) {
      s_sdcStatus(pI) <- list(index=indZero, vals=rep("u", length(indZero)))
    }
  }
  s_problemInstance(object) <- pI
  validObject(object)
  return(object)
})

setMethod("c_quick_suppression", signature=c("sdcProblem", "list"), definition=function(object, input) {
  freq <- id <- sdcStatus <- weights <- NULL
  verbose <- input$verbose
  detectSingletons <- input$detectSingletons
  pI <- g_problemInstance(object)
  indices <- g_partition(object)$indices
  dimInfo <- g_dimInfo(object)
  strInfo <- g_str_info(dimInfo)
  vNames <- g_varname(dimInfo)

  if (verbose) {
    cat("calculating subIndices (this may take a while) ...")
  }

  dat <- as.data.table(cpp_splitByIndices(g_strID(pI), strInfo))
  setnames(dat, vNames)
  dat[,id:=1:nrow(dat)]
  dat[,freq:=g_freq(pI)]
  dat[,weights:=g_weight(pI)]
  dat[,sdcStatus:=g_sdcStatus(pI)]
  dimVars <- match(vNames, names(dat))
  nDims <- length(dimVars)
  freqInd <- match("freq", colnames(dat))
  if (length(vNames)==1) {
    combs <- combn(vNames, 1)
  } else {
    combs <- combn(vNames, length(vNames)-1)
  }

  tmpIndices <- rep(NA, length(vNames))

  nrGroups <- length(indices)
  subIndices <- list(); length(subIndices) <- nrGroups

  for (group in 1:nrGroups) {
    nrTabs <- length(indices[[group]])
    subIndices[[group]] <- list()
    length(subIndices[[group]]) <- nrTabs
    for (tab in 1:nrTabs) {
      subDat <- dat[indices[[group]][[tab]],]
      # only one dimension!
      if (ncol(combs) == 1) {
        subDat$ind_1_tmp <- 1
        tmpIndices[1] <- ncol(subDat)
      } else {
        for (i in 1:ncol(combs)) {
          setkeyv(subDat, combs[,i])
          cn <- paste0("ind_",i,"_tmp")
          expr <- parse(text = paste0(cn, ":=.GRP"))
          subDat[,eval(expr), by=key(subDat)]
          tmpIndices[i] <- ncol(subDat)
        }
      }
      setkeyv(subDat, vNames)
      subIndices[[group]][[tab]] <- as.list(subDat[,tmpIndices, with=F])
    }
  }
  if (verbose) {
    cat("[done]\n");
  }

  if (detectSingletons==TRUE) {
    if (verbose) {
      cat("start singleton detection procedure!\n")
    }
    res <- singletonDetectionProcedure(dat=dat, indices=indices, subIndices=subIndices)
    if (verbose) {
      cat("singleton-detection procedure finished with",res$nrAddSupps,"additional suppressions!\n")
    }
    dat <- res$dat; rm(res)
  }

  res <- greedyMultDimSuppression(dat, indices, subIndices, dimVars, verbose=verbose)
  if (verbose) {
    cat("finishing output...")
  }
  s_sdcStatus(pI) <- list(index=res$id, vals=res$sdcStatus)
  s_problemInstance(object) <- pI
  if (verbose) {
    cat("[done]\n")
  }
  invisible(list(object=object, zstatus=res$status_z))
})

setMethod("c_cellID", signature=c("sdcProblem", "list"), definition=function(object, input) {
  para.names <- input[[1]]
  para.codes <- input[[2]]
  para.verbose <- input[[3]]

  pI <- g_problemInstance(object)
  dimInfoObj <- g_dimInfo(object)

  vNames <- g_varname(dimInfoObj)
  vIndex <- g_pos_index(dimInfoObj)

  indexVar <- match(para.names, vNames)
  strInfo <- g_str_info(dimInfoObj)
  dimInfo <- g_dim_info(dimInfoObj)

  # calculate original codes
  codesDefault <- lapply(1:length(strInfo), function(x) {
    mySplit(g_strID(pI), strInfo[[x]][1]:strInfo[[x]][2])
  })
  codesOrig <- list()
  for ( i in 1:length(codesDefault) ) {
    codesOrig[[i]] <- c_match_orig_codes(object=dimInfo[[i]], input=codesDefault[[i]])
  }

  if ( length(input) != 3 ) {
    stop("c_cellID:: length of argument 'input' must equal 3!\n")
  }
  if ( length(para.names) != length(para.codes) ) {
    stop("c_cellID:: check argument 'input'!\n")
  }
  if ( !all(para.names %in% vNames) ) {
    stop("c_cellID:: check variable names in 'input[[1]]'!\n")
  }
  if ( !is.logical(para.verbose) ) {
    stop("c_cellID:: argument in 'input[[3]]' must be logical!\n")
  }

  cellID <- 1:g_nrVars(pI)
  for ( i in seq_along(para.names) ) {
    cellID <- intersect(cellID, which(!is.na(match(as.character(codesOrig[[indexVar[i]]]), para.codes[i]))))
  }
  if ( length(cellID) != 1) {
    stop("c_cellID:: check argument 'input' -> 0 or > 1 cells identified!\n")
  }
  return(cellID)
})

setMethod("c_finalize", signature=c("sdcProblem", "list"), definition=function(object, input) {
  Freq <- NULL
  time.start <- proc.time()

  pI <- g_problemInstance(object)
  dI <- g_dimInfo(object)
  levelObj <- g_dim_info(dI)
  strInfo <- g_str_info(dI)

  sdcStatus <- g_sdcStatus(pI)

  nrNonDuplicatedCells <- g_nrVars(pI)
  nrPrimSupps <- length(which(sdcStatus == 'u'))
  nrSecondSupps <- length(which(sdcStatus == 'x'))
  nrPublishableCells <- length(which(sdcStatus %in% c('z','s')))

  # merge codes and labels
  codesOriginal <- list()
  strIDs <- g_strID(pI)
  for ( i in seq_along(levelObj) ) {
    codesDefault <- mySplit(strIDs, strInfo[[i]][1]:strInfo[[i]][2])
    codesOriginal[[i]] <- c_match_orig_codes(object=levelObj[[i]], input=codesDefault)
  }
  out <- as.data.table(codesOriginal)
  setnames(out, g_varname(dI))
  out[,Freq:=g_freq(pI)]

  numVars <- g_numVars(pI)
  if (!is.null(numVars) ) {
    data.obj <- g_dataObj(object)
    nV <- as.data.table(numVars)
    setnames(nV, colnames(g_raw_data(data.obj))[g_numvar_ind(data.obj)])
    out <- cbind(out, nV)
  }
  out[,sdcStatus:=g_sdcStatus(pI)]
  out[sdcStatus=="z", sdcStatus:="s"]

  # add duplicates
  hasDups <- sapply(1:length(levelObj), function(x) {
    g_has_dups(levelObj[[x]])
  })
  if (any(hasDups)) {
    for (i in which(hasDups==TRUE)) {
      dups <- g_dups(levelObj[[i]])
      dupsUp <- g_dups_up(levelObj[[i]])
      runInd <- TRUE
      while (runInd) {
        add <- list(); length(add) <- length(dups)
        for (j in 1:length(dups)) {
          if (!is.na(dups[j])) {
            cmd <- paste0("sub <- out[",names(out)[i],"=='",dupsUp[j],"',]")
            eval(parse(text=cmd))

            if (nrow(sub) > 0) {
              sub[[i]] <- dups[j]
              add[[j]] <- sub
              dups[j] <- dupsUp[j] <- NA
            }
          }
        }
        add <- rbindlist(add)
        out <- rbind(out, add); rm(add)
        if (all(is.na(dups))) {
          runInd <- FALSE
        }
      }
    }
  }
  safeObj <- new("safeObj",
    finalData=out,
    dimInfo=dI,
    nrNonDuplicatedCells=nrNonDuplicatedCells,
    nrPrimSupps=nrPrimSupps,
    nrSecondSupps=nrSecondSupps,
    nrPublishableCells=nrPublishableCells,
    suppMethod=input$method,
    elapsedTime=g_elapsedTime(object) + (proc.time()-time.start)[3]
  )
  return(safeObj)
})

setMethod("c_contributing_indices", signature=c("sdcProblem", "list"), definition=function(object, input) {
  strID <- input[[1]]
  dataObj <- g_dataObj(object)
  dimInfoObj <- g_dimInfo(object)
  dimInfo <- g_dim_info(dimInfoObj)
  pI <- g_problemInstance(object)

  if (!strID %in% g_strID(pI)) {
    stop("c_contributing_indices:: strID not found in the current problem!\n")
  }
  dims <- lapply(dimInfo, function(x) {
    g_dims(x)
  })
  indexVec <- which(g_str_id(dimInfoObj)==strID)
  # some (sub)totals need to be considered
  if (length(indexVec) == 0) {
    levInfo <- list()
    for ( z in 1:length(dimInfo) ) {
      subLevel <- substr(strID, g_str_info(dimInfoObj)[[z]][1], g_str_info(dimInfoObj)[[z]][2])
      if (sum(as.numeric(subLevel)) == 0) {
        levInfo[[z]] <- sort(unique(unlist(dims[[z]])))
      } else {
        orderInd <- unlist(lapply(dims[[z]], function(x) { match(subLevel, x)}))
        if (min(orderInd, na.rm=TRUE) == 1) {
          levInfo[[z]] <- dims[[z]][[which(orderInd==1)]]
        } else {
          levInfo[[z]] <- subLevel
        }
      }
    }
    cellIndex <- pasteStrVec(unlist(expand.grid(levInfo)), length(levInfo))
    indexVec <- which(g_str_id(dimInfoObj) %in% cellIndex)
  }
  return(indexVec)
})
