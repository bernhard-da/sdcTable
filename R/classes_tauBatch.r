setClassUnion("characterOrNULL", c("character", "NULL"))
setClassUnion("sdcProblemOrNULL", c("sdcProblem", "NULL"))
setClassUnion("listOrNull", c("list", "NULL"))

tau_BatchObj <- setClass("tau_BatchObj",
  slots = c(
    path = "characterOrNULL",
    id = "characterOrNULL",
    logbook = "characterOrNULL",
    microdata = "characterOrNULL",
    metadata = "characterOrNULL",
    table = "characterOrNULL",
    safetyrules = "characterOrNULL",
    readInput = "characterOrNULL",
    solver = "listOrNull",
    suppress = "characterOrNULL",
    writetable = "characterOrNULL",
    is_table = "logical",
    obj = "sdcProblemOrNULL"
  ),

  # Set the default values for the slots. (optional)
  prototype = list(
    path = NULL,
    id = NULL,
    logbook = NULL,
    microdata = NULL,
    metadata = NULL,
    table = NULL,
    safetyrules = NULL,
    readInput = NULL,
    solver = NULL,
    suppress = NULL,
    writetable = NULL,
    is_table = FALSE,
    obj = NULL
  ),
  validity = function(object) {
    if (length(object@is_table) != 1) {
      stop("length(is_table) != 1\n")
    }
    if (!is.null(object@solver)) {
      if (object@solver$solver == "CPLEX" &&
          !file.exists(object@solver$license)) {
        stop("No valid licensefile given!\n")
      }
    }
    return(TRUE)
  }
)


## define generics
setGeneric(
  name = "setPath",
  def = function(obj, f) {
    standardGeneric("setPath")
  }
)
setGeneric(
  name = "setId",
  def = function(obj, f) {
    standardGeneric("setId")
  }
)
setGeneric(
  name = "setLogbook",
  def = function(obj, f) {
    standardGeneric("setLogbook")
  }
)
setGeneric(
  name = "setMicrodata",
  def = function(obj, f) {
    standardGeneric("setMicrodata")
  }
)
setGeneric(
  name = "setMetadata",
  def = function(obj, f) {
    standardGeneric("setMetadata")
  }
)
setGeneric(
  name = "setTable",
  def = function(obj, f) {
    standardGeneric("setTable")
  }
)
setGeneric(
  name = "setSafetyrules",
  def = function(obj, f) {
    standardGeneric("setSafetyrules")
  }
)
setGeneric(
  name = "setReadInput",
  def = function(obj, f) {
    standardGeneric("setReadInput")
  }
)
setGeneric(
  name = "setSolver",
  def = function(obj, f) {
    standardGeneric("setSolver")
  }
)
setGeneric(
  name = "setSuppress",
  def = function(obj, f) {
    standardGeneric("setSuppress")
  }
)
setGeneric(
  name = "setWritetable",
  def = function(obj, f) {
    standardGeneric("setWritetable")
  }
)
setGeneric(
  name = "setIs_table",
  def = function(obj, f) {
    standardGeneric("setIs_table")
  }
)
setGeneric(
  name = "setObj",
  def = function(obj, f) {
    standardGeneric("setObj")
  }
)

setMethod(
  f = "setPath",
  signature = "tau_BatchObj",
  definition = function(obj, f) {
    obj@path <- f
    validObject(obj)
    return(obj)
  }
)
setMethod(
  f = "setId",
  signature = "tau_BatchObj",
  definition = function(obj, f) {
    obj@id <- f
    validObject(obj)
    return(obj)
  }
)

setMethod(
  f = "setLogbook",
  signature = "tau_BatchObj",
  definition = function(obj, f) {
    obj@logbook <- f
    validObject(obj)
    return(obj)
  }
)

setMethod(
  f = "setMicrodata",
  signature = "tau_BatchObj",
  definition = function(obj, f) {
    obj@microdata <- f
    validObject(obj)
    return(obj)
  }
)

setMethod(
  f = "setMetadata",
  signature = "tau_BatchObj",
  definition = function(obj, f) {
    obj@metadata <- f
    validObject(obj)
    return(obj)
  }
)

setMethod(
  f = "setTable",
  signature = "tau_BatchObj",
  definition = function(obj, f) {
    obj@table <- f
    validObject(obj)
    return(obj)
  }
)

setMethod(
  f = "setSafetyrules",
  signature = "tau_BatchObj",
  definition = function(obj, f) {
    obj@safetyrules <- f
    validObject(obj)
    return(obj)
  }
)

setMethod(
  f = "setReadInput",
  signature = "tau_BatchObj",
  definition = function(obj, f) {
    obj@readInput <- f
    validObject(obj)
    return(obj)
  }
)

setMethod(
  f = "setSolver",
  signature = "tau_BatchObj",
  definition = function(obj, f) {
    obj@solver <- f
    validObject(obj)
    return(obj)
  }
)

setMethod(
  f = "setSuppress",
  signature = "tau_BatchObj",
  definition = function(obj, f) {
    obj@suppress <- f
    validObject(obj)
    return(obj)
  }
)

setMethod(
  f = "setWritetable",
  signature = "tau_BatchObj",
  definition = function(obj, f) {
    obj@writetable <- f
    validObject(obj)
    return(obj)
  }
)

setMethod(
  f = "setIs_table",
  signature = "tau_BatchObj",
  definition = function(obj, f) {
    obj@is_table <- f
    validObject(obj)
    return(obj)
  }
)

setMethod(
  f = "setObj",
  signature = "tau_BatchObj",
  definition = function(obj, f) {
    obj@obj <- f
    validObject(obj)
    return(obj)
  }
)

setGeneric(
  name = "writeBatchFile",
  def = function(obj) {
    standardGeneric("writeBatchFile")
  }
)
setMethod(
  f = "writeBatchFile",
  signature = c("tau_BatchObj"),
  definition = function(obj) {
    is_table <- obj@is_table
    path <- obj@path

    cmds <- list()

    ## metainformation
    cmds <- append(cmds, "//This batch file was generated by sdcTable")
    cmds <- append(cmds, paste("//Date:", Sys.time()))
    cmds <- append(cmds, "//")

    ## tauArgus batch file
    f_log <- normalizePath(
      file.path(path, obj@logbook),
      winslash = "/",
      mustWork = FALSE
    )
    f_data <- normalizePath(
      file.path(path, obj@microdata),
      winslash = "/",
      mustWork = TRUE
    )
    f_metadata <- normalizePath(
      file.path(path, obj@metadata),
      winslash = "/",
      mustWork = TRUE
    )
    cmds <- append(cmds, paste("<LOGBOOK>", dQuote(f_log)))
    if (is_table) {
      cmds <- append(cmds, paste("<OPENTABLEDATA>", dQuote(f_data)))
    } else {
      cmds <- append(cmds, paste("<OPENMICRODATA>", dQuote(f_data)))
    }
    cmds <- append(cmds, paste("<OPENMETADATA>", dQuote(f_metadata)))
    cmds <- append(cmds, paste("<SPECIFYTABLE>", obj@table))
    cmds <- append(cmds, paste("<SAFETYRULE>", obj@safetyrules))
    cmds <- append(cmds, obj@readInput)

    solver <- slot(obj, "solver")$solver
    if (solver == "CPLEX") {
      f_license <- normalizePath(
        slot(obj, "solver")$license,
        winslash = "/",
        mustWork = TRUE
      )
      cmds <- append(cmds, paste0("<SOLVER> ", solver, ",", dQuote(f_license)))
    } else {
      cmds <- append(cmds, paste("<SOLVER>", solver))
    }

    cmds <- append(cmds, paste("<SUPPRESS>", obj@suppress))
    cmds <- append(cmds, paste("<WRITETABLE>", obj@writetable))

    f_batch <- generateStandardizedNames(
      path = obj@path,
      lab = paste0("batch_", obj@id),
      ext = ".arb"
    )
    cmds <- unlist(cmds)
    cmds[length(cmds)] <- paste0(cmds[length(cmds)], "\r")
    cat(cmds, sep = "\r\n", file = f_batch)
    invisible(f_batch)
  }
)
