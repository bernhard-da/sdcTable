##############################################
##### some general helper/check function #####
##############################################
# check if vector is integer
check_int <- function(x) {
  if (!is.numeric(x)) {
    return(FALSE)
  }
  all(na.omit(x)%%1==0)
}

# different checks
check_numeric <- function(inp) {
  is.numeric(inp)
}
check_character <- function(inp) {
  is.character(inp)
}
check_len <- function(inp, len=1) {
  length(inp)==len
}
check_range <- function(inp, v_min, v_max) {
  between(inp, v_min, v_max)
}

check_pos_number <- function(inp, v_min=1, v_max=100, default=1) {
  res <- check_numeric(inp) && check_int(inp) && check_len(inp) && check_range(inp, v_min, v_max)
  if (res) {
    return(inp)
  } else {
    return(default)
  }
}

# check variable input for batch-file
check_varinput <- function(obj, type, responsevar, shadowvar, costvar, requestvar, holdingvar) {
  if (!class(obj)=="sdcProblem") {
    stop("argument 'obj' must be of class 'sdcProblem'!\n")
  }
  stopifnot(type %in% c("microdata", "tabular"))

  if (is.null(responsevar)) {
    if (type=="microdata") {
      responsevar <- "<freq>"
    }
    if (type == "tabular") {
      responsevar <- "freq"
    }
  }

  df <- g_df(obj, addNumVars=TRUE)
  vNames <- get.sdcProblem(obj, "dimInfo")@vNames
  vNames_o <- paste0(vNames,"_o")
  # "allowed" variables
  cn <- setdiff(names(df), c("strID", vNames, vNames_o,"sdcStatus"))
  if (responsevar!="<freq>") {
    if (!check_character(responsevar) | !check_len(responsevar, 1)) {
      stop("argument 'responsevar' must be a character vector with 1 element!\n")
    }
    if (!responsevar %in% cn) {
      stop("non-valid variable selected for choice 'responsevar'.\n")
    }
  }

  if (!is.null(shadowvar)) {
    if (!check_character(shadowvar) | !check_len(shadowvar, 1)) {
      stop("argument 'shadowvar' must be a character vector with 1 element!\n")
    }
    if (!shadowvar %in% cn) {
      stop("non-valid variable selected for choice 'shadowvar'.\n")
    }
  } else {
    shadowvar <- ""
  }
  if (!is.null(costvar)) {
    if (!check_character(costvar) | !check_len(costvar, 1)) {
      stop("argument 'costvar' must be a character vector with 1 element!\n")
    }
    if (!costvar %in% cn) {
      stop("non-valid variable selected for choice 'costvar'.\n")
    }
  } else {
    costvar <- ""
  }

  if (!is.null(requestvar)) {
    if (!check_character(requestvar) | !check_len(requestvar, 1)) {
      stop("argument 'requestvar' must be a character vector with 1 element!\n")
    }
    if (!all(obj@dataObj@rawData[[requestvar]]  %in% c(0,1))) {
      stop("variable 'requestvar' must contain 0 and 1 only!\n")
    }
    if (!requestvar %in% cn) {
      stop("non-valid variable selected for choice 'requestvar'.\n")
    }
  } else {
    requestvar <- ""
  }

  if (!is.null(holdingvar)) {
    if (!check_character(holdingvar) | !check_len(holdingvar, 1)) {
      stop("argument 'holdingvar' must be a character vector with 1 element!\n")
    }
    if (!holdingvar %in% cn) {
      stop("non-valid variable selected for choice 'holdingvar'.\n")
    }
  } else {
    holdingvar <- ""
  }

  list(responsevar=responsevar, shadowvar=shadowvar, costvar=costvar,
    requestvar=requestvar, holdingvar=holdingvar)
}

# check suppression-method
check_suppmethod <- function(method) {
  if (!check_len(method, len=1) | !method %in% c("OPT","GH","MOD")) {
    stop("argument 'method' must be either 'OPT', 'MOD' or 'GH'!\n")
  }
  method
}


###############################################################################################################
##### helper-functions that are used to create hierarchy-files for tau-argus using an sdcProblem-instance #####
###############################################################################################################
## note: we do not need to create cdl (codelist-files) as they are only relevant for the GUI
## create required inputs
hrc <- function(obj) {
  convert_dim <- function(dim) {
    if (!class(dim)=="dimVar") {
      stop("convert_dim(): wrong input!\n")
    }
    codes <- get.dimVar(dim, "codesOriginal")
    levs <- get.dimVar(dim, "levels")-1
    ll <- sapply(1:length(codes), function(x) {
      str_sub(paste(rep("@", levs[x]), collapse=""), 2, -1)
    })

    hrc <- data.table(l=ll, codes=codes)[-1]
    hrc
  }

  dI <- get.sdcProblem(obj, "dimInfo")
  dims <- get.dimInfo(dI, "dimInfo")
  dim <- dims[[1]]
  res <- lapply(dims, function(x) {
    convert_dim(x)
  })
  res
}

## write to file (microdata)
write_hrc <- function(inp, fOut, nr_digits) {
  out <- l <- codes <- NULL
  unique_vals <- unique(inp$l)
  inp[,out:=""]
  for (vv in unique_vals) {
    if (vv=="") {
      inp[l==vv, out:=str_pad(codes, width=nr_digits, side="left")]
    } else {
      inp[l==vv, out:=paste0(l, str_pad(codes, nr_digits))]
    }
  }
  inp[,codes:=NULL]
  write.table(inp[,out], file=fOut, sep=" ", row.names=FALSE, col.names=FALSE, quote=FALSE,  eol="\r\n")
}

###################################################################################################################
#### helper-function to write data-files and metadata files required for tau-argus based on sdcProblem-objects ####
###################################################################################################################
## using microdata
create_microdata_and_metadata <- function(obj, verbose, digits=2, path=getwd(), ID, requestvar=NULL, holdingvar=NULL) {
  f1 <- generateStandardizedNames(path=NULL, lab=paste0("metadata_", ID), ext=".rda")
  f2 <- generateStandardizedNames(path=NULL, lab=paste0("microdata_", ID), ext=".asc")

  f_metadata <- file.path(path, f1)
  f_microdata <- file.path(path, f2)

  # get data
  dataObj <- get.sdcProblem(obj, "dataObj")
  mdat <- get.dataObj(dataObj, "rawData")
  vNames <- names(mdat)

  # convert num to integer
  ii <- which(sapply(mdat, class)%in% c("numeric","integer"))
  if (length(ii)>0) {
    cn <- vNames[ii]

    for (v in cn) {
      if (check_int(mdat[[v]])) {
        cmd <- paste0("mdat[,",v,":=as.integer(",v,")]")
        eval(parse(text=cmd))
      }
    }
  }

  if (!get.dataObj(dataObj, "isMicroData")) {
    stop("using complete tables (isMicroData==FALSE) is currently unsupported!\n")
  }

  bl <- "  "
  cmds <- list()
  cmds <- append(cmds, paste("<SEPARATOR>", dQuote(",")))

  required_digits <- c()

  # define output-matrix (required for fixed-file format)
  mat <- matrix("", nrow=nrow(mdat), ncol=1)

  # calculate hierarchiy-files
  hiercodes <- hrc(obj)

  # write microdata
  # 1: dim-vars
  # --> all are <RECODEABLE>
  # --> no codelist, we use the values in the hrc-file
  # --> all are <HIERARCHICAL>
  # --> no <HIERLEVELS> because we are using <HIERCODELIST>
  # --> all have <HIERCODELIST> (hrc-file) generated with hrc()
  # --> all have <HIERLEADSTRING>
  # --> none have <REQUEST> or <HOLDING>: currently not supported
  dim_vars <- vNames[get.dataObj(dataObj, "dimVarInd")]
  f_hrcs <- c()
  for (i in seq_along(dim_vars)) {
    vv <- dim_vars[i]

    startpos <- sum(required_digits)+i
    cur_dig <- max(nchar(as.character(mdat[[vv]])))
    required_digits <- c(required_digits, cur_dig)

    # write hierachy-files
    f_hrc <- generateStandardizedNames(path=path, lab=paste0("hier_",ID,"_",vv), ext=".hrc")
    f_hrcs <- c(f_hrcs, f_hrc)
    write_hrc(inp=hiercodes[[vv]], fOut=f_hrc, nr_digits=cur_dig)
    f_hrc <- normalizePath(f_hrc, winslash="/", mustWork=TRUE)
    cmds <- append(cmds, paste(vv, cur_dig)) # no missings allowed in sdcProblem-objects
    cmds <- append(cmds, paste(bl, "<RECODEABLE>"))
    cmds <- append(cmds, paste(bl, "<HIERCODELIST>", dQuote(f_hrc)))
    cmds <- append(cmds, paste(bl, "<HIERLEADSTRING>", dQuote("@")))
    cmds <- append(cmds, paste(bl, "<HIERARCHICAL>"))

    if (i==1) {
      mat[,1] <- str_pad(mdat[[vv]], width=cur_dig, side="left")
    } else {
      mat <- cbind(mat, str_pad(mdat[[vv]], width=cur_dig, side="left"))
    }
  }

  # 2: optionally (weight)
  # 3: numeric variables

  num_vars <- setdiff(colnames(mdat), c("freq", dim_vars))
  for (i in seq_along(num_vars)) {
    vv <- num_vars[i]

    has_decimals <- FALSE
    if (!check_int(mdat[[vv]])) {
      digits <- digits
      mdat[[vv]] <- formatC(mdat[[vv]], format="f", digits=digits)
      has_decimals <- TRUE
    }

    startpos <- sum(required_digits)+i
    cur_dig <- max(nchar(as.character(mdat[[vv]])))
    required_digits <- c(required_digits, cur_dig)
    cmds <- append(cmds, paste(vv, cur_dig, dQuote(str_pad("", width=cur_dig, pad="9"))))

    if (!is.null(requestvar) && requestvar==vv) {
      cmds <- append(cmds, paste(bl, "<REQUEST>", dQuote(1)))
    } else if (!is.null(holdingvar) && holdingvar==vv) {
      cmds <- append(cmds, paste(bl, "<HOLDING>"))
    } else {
      cmds <- append(cmds, paste(bl, "<NUMERIC>"))
      if (has_decimals) {
        cmds <- append(cmds, paste(bl, "<DECIMALS>", digits))
      }
      wV <- vNames[get.dataObj(dataObj, "sampWeightInd")]
      if (length(wV)>0) {
        if (vv==wV) {
          cmds <- append(cmds, paste(bl, "<WEIGHT>"))
        }
      }
    }
    mat <- cbind(mat, str_pad(mdat[[vv]], width=cur_dig, side="left"))
  }

  if (verbose) {
    cat("writing metadatafile to",shQuote(f_metadata),"\n")
  }
  cmds <- unlist(cmds)
  cmds[length(cmds)] <- paste0(cmds[length(cmds)],"\r")
  cat(cmds, sep="\r\n", file=f_metadata)

  if (verbose) {
    cat("writing microdatafile to",shQuote(f_microdata),"\n")
  }
  write.table(mat, file=f_microdata, sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE, eol="\r\n")
  invisible(list(f_md=f1, f_dat=f2, f_hrcs=f_hrcs))
}

## using tabular data
create_tabdata_and_metadata <- function(obj, verbose, responsevar, digits=2, path=getwd(), ID) {
  lpl <- upl <- sdcStatus <- NULL
  f1 <- generateStandardizedNames(path=NULL, lab=paste0("metadata_", ID), ext=".rda")
  f2 <- generateStandardizedNames(path=NULL, lab=paste0("tabdata_", ID), ext=".tab")

  f_metadata <- file.path(path, f1)
  f_tabdata <- file.path(path, f2)

  # get data
  dataObj <- get.sdcProblem(obj, "dataObj")
  pI <- get.sdcProblem(obj, "problemInstance")
  mdat <- copy(sdcProb2df(obj, addNumVars=TRUE, dimCodes="original"))
  mdat[,lpl:=get.problemInstance(pI, "LPL")]
  mdat[,upl:=get.problemInstance(pI, "UPL")]

  vNames <- names(mdat)

  # convert num to integer
  ii <- which(sapply(mdat, class)%in% c("numeric","integer"))
  if (length(ii)>0) {
    cn <- vNames[ii]

    for (v in cn) {
      if (check_int(mdat[[v]])) {
        cmd <- paste0("mdat[,",v,":=as.integer(",v,")]")
        eval(parse(text=cmd))
      }
    }
  }

  bl <- "  "
  cmds <- list()
  cmds <- append(cmds, paste("<SEPARATOR>", dQuote(",")))
  cmds <- append(cmds, paste("<SAFE>", dQuote("S")))
  cmds <- append(cmds, paste("<UNSAFE>", dQuote("U")))
  cmds <- append(cmds, paste("<PROTECT>", dQuote("Z")))

  required_digits <- c()

  # define output-matrix (required for fixed-file format)
  mat <- matrix("", nrow=nrow(mdat), ncol=1)

  # calculate hierarchiy-files
  hiercodes <- hrc(obj)

  # write microdata
  # 1: dim-vars
  # --> all are <RECODEABLE>
  # --> no codelist, we use the values in the hrc-file
  # --> all are <HIERARCHICAL>
  # --> no <HIERLEVELS> because we are using <HIERCODELIST>
  # --> all have <HIERCODELIST> (hrc-file) generated with hrc()
  # --> all have <HIERLEADSTRING>
  # --> none have <REQUEST> or <HOLDING>: currently not supported
  dim_vars <- names(hiercodes)
  f_hrcs <- c()
  for (i in seq_along(dim_vars)) {
    vv <- dim_vars[i]
    startpos <- sum(required_digits)+i

    tot_code <- get.sdcProblem(obj, "dimInfo")@dimInfo[[i]]@codesOriginal[1]
    nch_tot <- nchar(tot_code)

    cur_dig <- max(nchar(as.character(mdat[[vv]])))
    cur_dig <- max(cur_dig, nch_tot)
    required_digits <- c(required_digits, cur_dig)

    # write hiercode-file
    f_hrc <- generateStandardizedNames(path=path, lab=paste0("hier_",ID,"_",vv), ext=".hrc")
    f_hrcs <- c(f_hrcs, f_hrc)
    write_hrc(inp=hiercodes[[vv]], fOut=f_hrc, nr_digits=0)
    f_hrc <- normalizePath(f_hrc, winslash="/", mustWork=TRUE)

    cmds <- append(cmds, paste(vv)) # no missings allowed in sdcProblem-objects
    cmds <- append(cmds, paste(bl, "<RECODEABLE>"))
    cmds <- append(cmds, paste(bl, "<TOTCODE>", dQuote(tot_code)))
    cmds <- append(cmds, paste(bl, "<HIERCODELIST>", dQuote(f_hrc)))
    cmds <- append(cmds, paste(bl, "<HIERLEADSTRING>", dQuote("@")))
    cmds <- append(cmds, paste(bl, "<HIERARCHICAL>"))

    if (i==1) {
      mat[,1] <- str_pad(mdat[[vv]], width=cur_dig, side="left")
    } else {
      mat <- cbind(mat, str_pad(mdat[[vv]], width=cur_dig, side="left"))
    }
  }

  # 2b: numeric variables, frequency, lpl/upl
  num_vars <- c("sdcStatus", responsevar, "lpl","upl")

  for (i in seq_along(num_vars)) {
    vv <- num_vars[i]
    has_decimals <- FALSE
    if (!check_int(mdat[[vv]])) {
      digits <- digits
      mdat[[vv]] <- formatC(mdat[[vv]], format="f", digits=digits)
      has_decimals <- TRUE
    }

    # is non-integer!
    if (vv=="sdcStatus") {
      has_decimals <- FALSE
      mdat[,sdcStatus:=toupper(sdcStatus)] # bug in tau-argus?
    }

    if (any(is.na(mdat[[vv]]))) {
      stop("missing value detected in variable ",dQuote(vv),". This is currently not supported!\n")
    }
    startpos <- sum(required_digits)+i
    cur_dig <- max(nchar(as.character(mdat[[vv]])))
    required_digits <- c(required_digits, cur_dig)

    if (vv=="freq") {
      cmds <- append(cmds, vv)
      cmds <- append(cmds, paste(bl, "<NUMERIC>"))
    } else if (vv=="lpl") {
      cmds <- append(cmds, vv)
      cmds <- append(cmds, paste(bl, "<NUMERIC> <LOWERPL>"))
    } else if (vv=="upl") {
      cmds <- append(cmds, vv)
      cmds <- append(cmds, paste(bl, "<NUMERIC> <UPPERPL>"))
    } else if (vv=="sdcStatus") {
      cmds <- append(cmds, vv)
      cmds <- append(cmds, paste(bl, "<STATUS>"))
    }
    # do not append any other variable since batch-input assumes one numeric variable

    if (has_decimals) {
      cmds <- append(cmds, paste(bl, "<DECIMALS>", digits))
    }
    wV <- vNames[get.dataObj(dataObj, "sampWeightInd")]
    if (length(wV)>0) {
      if (vv==wV) {
        cmds <- append(cmds, paste(bl, "<WEIGHT>"))
      }
    }
    mat <- cbind(mat, str_pad(mdat[[vv]], width=cur_dig, side="left"))
  }

  if (verbose) {
    cat("writing metadatafile to",shQuote(f_metadata),"\n")
  }
  cmds <- unlist(cmds)
  cmds[length(cmds)] <- paste0(cmds[length(cmds)],"\r")
  cat(cmds, sep="\r\n", file=f_metadata)

  if (verbose) {
    cat("writing tabular data file to",shQuote(f_tabdata),"\n")
  }
  write.table(mat, file=f_tabdata, sep=",", row.names=FALSE, col.names=FALSE, quote=FALSE, eol="\r\n")
  invisible(list(f_md=f1, f_dat=f2, f_hrcs=f_hrcs))
}


##################################################################################
##### helper-functions to deal with primary suppression rules (safety-rules) #####
##################################################################################
# FREQ(n, rg)
sf_freq <- function(n, rg) {
  if (is.null(n) || !check_int(n) || !check_range(n, 1, Inf)) {
    stop("sf_freq(): check parameter 'n'\n")
  }
  if (is.null(rg) || !check_int(rg) || !check_range(rg, 1, 99)) {
    stop("sf_freq(): check parameter 'rg'\n")
  }
  paste0("FREQ(",n,",",rg,")")
}

# P(p, n)
sf_p <- function(p, n) {
  if (is.null(p) || !check_int(p) || !check_range(p, 1, 99)) {
    stop("sf_p(): check parameter 'p'\n")
  }
  if (is.null(n) || !check_int(n) || !check_range(n, 1, Inf)) {
    stop("sf_p(): check parameter 'n'\n")
  }
  paste0("P(",p,",",n,")")
}

# NK(n, k)
sf_nk <- function(n, k) {
  if (is.null(n) || !check_int(n) || !check_range(n, 1, Inf)) {
    stop("sf_nk(): check parameter 'n'\n")
  }
  if (is.null(k) || !check_int(k) || !check_range(k, 1, 99)) {
    stop("sf_nk(): check parameter 'k'\n")
  }
  paste0("NK(",n,",",k,")")
}

# ZERO(rg)
sf_zero <- function(rg) {
  if (is.null(rg) || !check_int(rg) || !check_range(rg, 1, 99)) {
    stop("sf_zero(): check parameter 'rg'\n")
  }
  paste0("ZERO(",rg,")")
}

# MIS(val)
sf_mis <- function(val) {
  if (is.null(val) || !check_int(val) || !check_len(val, 1) || !val %in% c(0,1)) {
    stop("sf_mis(): check parameter 'val'\n")
  }
  paste0("MIS(",val,")")
}

# WGT(val)
sf_wgt <- function(val) {
  if (is.null(val) || !check_int(val) || !check_len(val, 1) || !val %in% c(0,1)) {
    stop("sf_wgt(): check parameter 'val'\n")
  }
  paste0("WGT(",val,")")
}

# MAN(val)
sf_man <- function(val) {
  if (is.null(val) || !check_int(val) || !check_len(val, 1) || !check_range(val, 1, 99)) {
    stop("sf_man(): check parameter 'val'\n")
  }
  paste0("MAN(",val,")")
}


# REQ(Percent1, Percent2, SafetyMargin)
#sf_req <- function(p1, p2, margin) {
#  if (is.null(rg) || !check_int(rg) || !check_range(rg, 1, 99)) {
#    stop("sf_req(): check parameter 'rg'\n")
#  }
#  paste0("REQ(",p1,",",p2,",",margin,")")
#}


## create safety-rules
# todo: REQ()
srule <- function(type, ...) {
  args <- list(...)
  type <- tolower(type)
  if (type=="p") {
    return(sf_p(p=args$p, n=args$n))
  } else if (type=="nk") {
    return(sf_nk(n=args$n, k=args$k))
  } else if (type=="freq") {
    return(sf_freq(n=args$n, rg=args$rg))
  } else if (type=="zero") {
    return(sf_zero(rg=args$rg))
  } else if (type=="mis") {
    return(sf_mis(val=args$val))
  } else if (type=="wgt") {
    return(sf_wgt(val=args$val))
  } else if (type=="man") {
    return(sf_man(val=args$val))
  } else {
    stop("invalid choice of 'type'!\n")
  }
}

## check for invalid inputs
check_primrules <- function(primSuppRules, responsevar) {
  out <- list(); length(out) <- length(primSuppRules)
  types <- sapply(primSuppRules, function(x) x$type)

  # at least one suppression-rule is required
  if (length(types)<1) {
    return(NA)
  }

  if (responsevar!="<freq>") {
    # max 2 NK() or P()-rules for magnitute-tables possible
    if (sum(types=="nk")>=2) {
      return(NA)
    }
    if (sum(types=="p")>=2) {
      return(NA)
    }
  } else {
    # in case of frequency tables, only one FREQ()-rule is possible
    if (length(types)>1 || primSuppRules[[1]]$type!="freq") {
      return(NA)
    }
  }

  for (i in 1:length(primSuppRules)) {
    res <- try(do.call(srule, primSuppRules[[i]]), silent=TRUE)
    if ("try-error" %in% class(res)) {
      return(NA)
    }
    out[[i]] <- res
  }
  out
}

#primSuppRules <- list()
#primSuppRules[[1]] <- list(type="freq", n=5, rg=20)
#primSuppRules[[2]] <- list(type="p", n=5, p=20)
#primSuppRules[[3]] <- list(type="nk", n=5, k=20.5)
#primSuppRules[[4]] <- list(type="mis", val=1)
#primSuppRules[[5]] <- list(type="wgt", val=1)
#primSuppRules[[6]] <- list(type="man", val=25)

#check_primrules(primSuppRules, responsevar="asf")

##############################################################
##### helper-function to generate standardized filenames #####
##############################################################
generateStandardizedNames <- function(path=tempdir(), lab, ext=".log") {
  if (is.null(path)) {
    return(paste0(lab, ext))
  } else {
    return(file.path(path, paste0(lab, ext)))
  }
}


#######################################################
##### helper-function to convert sdc-status codes #####
#######################################################
## convert status-codes from tau-argus (manual page 103) to sdcTable codes or from sdcTable to tau-argus
recode_sdcStati <- function(status_in, is_argus=TRUE) {
  status_out <- NULL
  dt <- data.table(status_in=status_in, status_out="")
  if (is_argus==TRUE) {
    dt[status_in %in% 1:2, status_out:="s"]
    dt[status_in %in% 3:9, status_out:="u"]
    dt[status_in %in% 10, status_out:="s"] # protected -> is could also be "x"?
    dt[status_in %in% 11:12, status_out:="x"]
    dt[status_in %in% 13:14, status_out:="z"]
  } else {
    dt[status_in=="z", status_out:=14]
    dt[status_in=="s", status_out:=1]
    dt[status_in=="u", status_out:=3]
    dt[status_in=="x", status_out:=11]
  }
  dt$status_out
}


#########################################################
##### helper-function to read in tau-argus solution #####
#########################################################
## reads an solution written to file from tau-argus into R
read_ArgusSolution <- function(fIn) {
  if (!file.exists(fIn)) {
    stop("Output file",dQuote(fIn),"does not exist!\n")
  }
  fread(fIn)
}


###########################################################
##### helper-functions to create tau_BatchObj-objects #####
###########################################################

## based on microdata input
tauBatchInput_microdata <- function(obj,
    verbose,
    path=getwd(),
    solver="FREE",
    method,
    primSuppRules,
    responsevar="<freq>",
    shadowvar=NULL,
    costvar=NULL,
    requestvar=NULL,
    holdingvar=NULL, ...) {

  args <- list(...)

  # create and check variable-input
  vars <- check_varinput(obj, type="microdata", responsevar, shadowvar, costvar, requestvar, holdingvar)
  responsevar <- vars$responsevar
  shadowvar <- vars$shadowvar
  costvar <- vars$costvar
  requestvar <- vars$requestvar
  holdingvar <- vars$holdingvar

  ## check argument 'method'
  method <- check_suppmethod(method)

  batchObj <- tau_BatchObj()
  batchID <- paste(sample(c(letters, LETTERS, 0:9), 10, replace=TRUE), collapse="")
  batchObj <- setId(batchObj, batchID)
  batchObj <- setObj(batchObj, obj)

  if (!file.exists(path)) {
    res <- dir.create(path)
  }

  batchObj <- setPath(batchObj, path)

  f_log <- generateStandardizedNames(path=NULL, lab=paste0("arguslog_", batchID), ext=".log")
  f_tab <- generateStandardizedNames(path=NULL, lab=paste0("tabout_", batchID), ext=".txt")

  res <- create_microdata_and_metadata(obj, verbose=verbose, digits=2, path=path, ID=batchID,
    requestvar=requestvar, holdingvar=holdingvar)

  ## logfile
  batchObj <- setLogbook(batchObj, f=f_log)

  ## microdata-file
  batchObj <- setMicrodata(batchObj, f=res$f_dat)

  ## metadata-file
  batchObj <- setMetadata(batchObj, f=res$f_md)

  ## specify the table
  p1 <- paste0(dQuote(obj@dimInfo@vNames), collapse="")
  lambda <- 1.0
  if (costvar!="") {
    if (!is.null(args$lambda)) {
      lambda <- args$lambda
    }
    tablestr <- paste0(p1,"|",dQuote(responsevar),"|",dQuote(shadowvar),"|",dQuote(costvar),",",lambda)
  } else {
    tablestr <- paste0(p1,"|",dQuote(responsevar),"|",dQuote(shadowvar),"|",dQuote(costvar),"")
  }
  batchObj <- setTable(batchObj, tablestr)

  ## primary suppression
  rules <- check_primrules(primSuppRules, responsevar)
  if (is.na(rules[1])) {
    stop("Invalid primary-suppression rules specified!\n")
  }
  batchObj <- setSafetyrules(batchObj, paste0(paste(unlist(rules), collapse="|"),"|"))

  ## read the data
  batchObj <- setReadInput(batchObj, "<READMICRODATA>")

  ## solver
  license <- NULL
  if (solver=="CPLEX") {
    license <- args$licensefile
  }
  batchObj <- setSolver(batchObj, list(solver=solver, license=license))

  # modular/hitas
  if (method=="MOD") {
    MaxTimePerSubtable <- check_pos_number(inp=args$MaxTimePerSubtable, default=1)
    SingleSingle <- check_pos_number(inp=args$SingleSingle, default=1, v_min=0, v_max=1)
    SingleMultiple <- check_pos_number(inp=args$SingleMultiple, default=1, v_min=0, v_max=1)
    MinFreq <- check_pos_number(inp=args$MinFreq, default=1, v_min=0, v_max=1)
    suppstr <- paste0("MOD(1,",MaxTimePerSubtable,",",SingleSingle,",",SingleMultiple,",",MinFreq,")")
  }
  # hypercube
  if (method=="GH") {
    # A priori Bounds Percentage
    # ModelSize 0 = normal, 1 indicates the large model.
    # ApplySingleton: 1 = yes, 0 = no; default = yes if the table has frequency-information, no if not
    BoundPercentage <- check_pos_number(inp=args$BoundPercentage, default=1, v_min=1, v_max=100)
    ModelSize <- check_pos_number(inp=args$ModelSize, default=1, v_min=0, v_max=1)
    ApplySingleton <- check_pos_number(inp=args$ApplySingleton, default=1, v_min=0, v_max=1)
    suppstr <- paste0("GH(1,",BoundPercentage,",",ModelSize,",",ApplySingleton,")")
  }
  # optimal
  if (method=="OPT") {
    MaxComputingTime <- check_pos_number(inp=args$MaxComputingTime, default=1, v_min=1, v_max=Inf)
    suppstr <- paste0("OPT(1,",MaxComputingTime,")")
  }
  batchObj <- setSuppress(batchObj, suppstr)

  ## output-table
  f_tab <- normalizePath(file.path(path, f_tab), winslash="/", mustWork=FALSE)
  batchObj <- setWritetable(batchObj, paste0("(1, 3, AS+FL+, ", dQuote(f_tab),")"))
  invisible(batchObj)
}

## based on tabular input
tauBatchInput_table <- function(obj,
    verbose,
    path=getwd(),
    solver="FREE",
    method,
    responsevar="freq",
    shadowvar=NULL,
    costvar=NULL, ...) {

  args <- list(...)

  # create and check variable-input
  vars <- check_varinput(obj, type="tabular", responsevar, shadowvar, costvar,
    requestvar=NULL, holdingvar=NULL)
  responsevar <- vars$responsevar
  shadowvar <- vars$shadowvar
  costvar <- vars$costvar

  # check argument 'method'
  method <- check_suppmethod(method)

  batchObj <- tau_BatchObj()
  batchObj <- setIs_table(batchObj, TRUE)
  batchObj <- setObj(batchObj, obj)

  batchID <- paste(sample(c(letters, LETTERS, 0:9), 10, replace=TRUE), collapse="")
  batchObj <- setId(batchObj, batchID)

  if (!file.exists(path)) {
    res <- dir.create(path)
  }
  batchObj <- setPath(batchObj, path)

  f_log <- generateStandardizedNames(path=NULL, lab=paste0("arguslog_",batchID), ext=".log")
  f_tab <- generateStandardizedNames(path=NULL, lab=paste0("tabout_", batchID), ext=".txt")

  res <- create_tabdata_and_metadata(obj, verbose=verbose, responsevar=responsevar, digits=2, path=path, ID=batchID)

  ## logfile
  batchObj <- setLogbook(batchObj, f=f_log)

  ## microdata-file
  batchObj <- setMicrodata(batchObj, f=res$f_dat)

  ## metadata-file
  batchObj <- setMetadata(batchObj, f=res$f_md)

  ## specify the table
  p1 <- paste0(dQuote(obj@dimInfo@vNames), collapse="")
  lambda <- 1.0
  if (costvar!="") {
    if (!is.null(args$lambda)) {
      lambda <- args$lambda
    }
    tablestr <- paste0(p1,"|",dQuote(responsevar),"|",dQuote(shadowvar),"|",dQuote(costvar),",",lambda)
  } else {
    tablestr <- paste0(p1,"|",dQuote(responsevar),"|",dQuote(shadowvar),"|",dQuote(costvar),"")
  }
  batchObj <- setTable(batchObj, tablestr)

  ## todo, parametrize safety-range
  batchObj <- setSafetyrules(batchObj, "MAN(20)")

  ## read the data
  batchObj <- setReadInput(batchObj, "<READTABLE>")

  ## solver
  license <- NULL
  if (solver=="CPLEX") {
    license <- args$licensefile
  }
  batchObj <- setSolver(batchObj, list(solver=solver, license=license))

  # modular/hitas
  if (method=="MOD") {
    MaxTimePerSubtable <- check_pos_number(inp=args$MaxTimePerSubtable, default=1)
    SingleSingle <- check_pos_number(inp=args$SingleSingle, default=1, v_min=0, v_max=1)
    SingleMultiple <- check_pos_number(inp=args$SingleMultiple, default=1, v_min=0, v_max=1)
    MinFreq <- check_pos_number(inp=args$MinFreq, default=1, v_min=0, v_max=1)
    suppstr <- paste0("MOD(1,",MaxTimePerSubtable,",",SingleSingle,",",SingleMultiple,",",MinFreq,")")
  }
  # hypercube
  if (method=="GH") {
    # A priori Bounds Percentage
    # ModelSize 0 = normal, 1 indicates the large model.
    # ApplySingleton: 1 = yes, 0 = no; default = yes if the table has frequency-information, no if not
    BoundPercentage <- check_pos_number(inp=args$BoundPercentage, default=1, v_min=1, v_max=100)
    ModelSize <- check_pos_number(inp=args$ModelSize, default=1, v_min=0, v_max=1)
    ApplySingleton <- check_pos_number(inp=args$ApplySingleton, default=1, v_min=0, v_max=1)
    suppstr <- paste0("GH(1,",BoundPercentage,",",ModelSize,",",ApplySingleton,")")
  }
  # optimal
  if (method=="OPT") {
    MaxComputingTime <- check_pos_number(inp=args$MaxComputingTime, default=1, v_min=1, v_max=Inf)
    suppstr <- paste0("OPT(1,",MaxComputingTime,")")
  }
  batchObj <- setSuppress(batchObj, suppstr)

  ## output-table
  f_tab <- normalizePath(file.path(path, f_tab), winslash="/", mustWork=FALSE)
  batchObj <- setWritetable(batchObj, paste0("(1, 3, AS+FL+, ", dQuote(f_tab),")"))
  invisible(batchObj)
}


###############################################################################
##### helper-function to merge sdcTable input and solution from tau-argus #####
###############################################################################
## merges sdc-results from argus to sdcTable output
# if obj is NULL, only the solution fro tau-Argus is read in and variable names are extracted from batch-Input
combineInputs <- function(obj=NULL, res_argus, batchF) {
  extractVarnames <- function(batchF) {
    inp <- readLines(batchF, warn=FALSE)
    inp <- inp[grep("SPECIFY", inp)]
    inp <- unlist(strsplit(inp, "> "))[-1]
    inp <- unlist(strsplit(inp, '"'))
    to <- which(inp=="|")[1]
    inp <- inp[1:(to-1)]
    vNames <- inp[inp!=""]
    vNames
  }

  id <- sdcStatus_argus <- cellvalue <- NULL
  if (is.null(obj)) {
    dN <- extractVarnames(batchF)
  } else {
    dN <- obj@dimInfo@vNames
  }
  setnames(res_argus, c(dN, "cellvalue", "sdcStatus_argus"))

  # if obj is NULL, we do not have information on "duplicates"
  if (!is.null(obj)) {
    ## check and re-add dups
    ## due to ordering, the corresponding value is listed in the row above the dups
    orig <- g_df(obj, addDups=TRUE)
    for (i in 1:length(dN)) {
      cmd <- paste0("orig[,",dN[i],":=",dN[i],"_o]")
      eval(parse(text=cmd))
      cmd <- paste0("orig[,", dN[i],"_o:=NULL]")
      eval(parse(text=cmd))

      cmd <- paste0("res_argus[,", dN[i],":=str_trim(",dN[i],")]")
      eval(parse(text=cmd))
    }
    orig$id <- 1:nrow(orig)
    setkeyv(orig, dN)
    setkeyv(res_argus, dN)
    orig <- merge(orig, res_argus, all.x=TRUE)
    setkey(orig, id)
    orig[,c("id","strID"):=NULL]

    ii <- which(is.na(orig[,sdcStatus_argus]))
    if (length(ii)>0) {
      vv <- orig[ii-1,sdcStatus_argus]
      orig[ii,sdcStatus_argus:=vv]
      cv <- orig[ii-1,cellvalue]
      orig[ii,cellvalue:=cv]
    }
  } else {
    orig <- copy(res_argus)
    suppressWarnings(orig[,cellvalue:=as.numeric(cellvalue)])
  }

  orig[,sdcStatus_argus:=recode_sdcStati(sdcStatus_argus, is_argus=TRUE)]
  orig
}

# read file-path from batch-files
infoFromBatch <- function(batchF, typ="LOGBOOK") {
  inp <- readLines(batchF, warn=FALSE)
  inp <- inp[grep(typ, inp)]
  if (length(inp)==0) {
    stop(paste(dQuote(typ),"not found in batch-file",dQuote(batchF),"\n"))
  }
  filepath <- tail(unlist(strsplit(inp, " ")),1)
  filepath <- gsub('\"', "", filepath)

  if (typ=="WRITETABLE") {
    filepath <- gsub(")", "", filepath)
  }
  filepath
}

#' argusVersion
#'
#' returns the version and build number of a given tau-argus executable
#' specified in argument \code{exe}.
#'
#' @param exe a path to a tau-argus executable
#' @param verbose (logical) if \code{TRUE}, the version info and build number
#' of the given tau-argus executable will be printed.
#'
#' @return a list with two elements being the tau-argus version and the build-number.
#' @export
#'
#' @examples
#' \dontrun{
#' argusVersion(exe="C:\\Tau\\TauArgus.exe", verbose=TRUE)
#' }
argusVersion <- function(exe, verbose=FALSE) {
  fInfo <- tempfile(fileext=".txt")
  fBatch <- tempfile(fileext=".arb")
  fLog <- tempfile(fileext=".log")
  cat(paste("<LOGBOOK>", dQuote(fLog),"\n<VERSIONINFO>", dQuote(fInfo),"\n"), file=fBatch)

  cmd <- paste(shQuote(exe), fBatch)
  res <- suppressWarnings(system(cmd, intern=TRUE, ignore.stdout=TRUE, ignore.stderr=FALSE))
  if (!is.null(attributes(res)$status)) {
    stop("Please use an tau-argus version >= 4.1.6\n")
  }

  res <- readLines(fInfo, warn=FALSE)
  res <- unlist(strsplit(res, ";"))
  version <- str_trim(unlist(strsplit(res[1], ":"))[2])
  build <- str_trim(unlist(strsplit(res[2], ":"))[2])
  if (verbose) {
    cat(paste0("Tau-Argus version: ", version, " (Build: ",build,")\n"))
  }

  # cleanup
  xx <- file.remove(c(fInfo, fBatch, fLog))
  return(list(version=version, build=build))
}
