#' Write a problem in jj-format to a file
#'
#' This function allows to write a problem instance in
#' JJ-Format to a file.
#'
#' @param x an input produced by [createJJFormat()]
#' @param tabvar the name of the variable that will be used when producing the
#' problem in JJ format. It is possible to specify `freqs` (the default) or the
#' name of a numeric variable that was available in the [sdcProblem-class] object
#' used in [makeProblem()].
#' @param path a scalar character defining the name of the file that
#' should be written. This can be an absolute or relative URL; however the file
#' must not exist.
#' @param overwrite logical scalar, if `TRUE` the file specified in `path` will
#' be overwritten if it exists
#' @return invisibly the path to the file that was created.
#' @export
#' @md
#' @examples
# # use test data
#' data("microData1", package = "sdcTable")
#'
#' # create hierarchies
#' dimList <- list(
#'   region = hier_create(
#'     root = "Total",
#'     nodes = LETTERS[1:4]
#'   ),
#'   gender = hier_create(
#'     root = "Total",
#'     nodes = c("male", "female")
#'   )
#' )
#'
#' # create a problem instance
#' prob <- makeProblem(
#'   data = microData1,
#'   dimList = dimList,
#'   numVarInd = "val"
#' )
#'
#' # create suitable input for `writeJJFormat`
#' inp <- createJJFormat(prob); inp
#'
#' # write files to disk
#' # frequency table by default
#' writeJJFormat(inp, path = "prob_freqs.jj", overwrite = TRUE)
#'
#' # or using the numeric variable `val` previously specified
#' writeJJFormat(inp, tabvar = "val", path = "prob_val.jj", overwrite = TRUE)
writeJJFormat <- function(x, tabvar = "freqs", path = "out.jj", overwrite = FALSE) {
  if (!inherits(x, "jjformat")) {
    e <- "Invalid input. Please use `createJJFormat()`."
    stop(e, call. = FALSE)
  }

  if (!is_scalar_character(path)) {
    stop("`path` needs to be a scalar character.", call. = FALSE)
  }

  if (!is_scalar_logical(overwrite)) {
    stop("Argument `overwrite` must be scalar logical.", call. = FALSE)
  }

  if (!overwrite & file.exists(path)) {
    e <- c(
      "File", shQuote(path), "exists.",
      "Please remove it, specify another path or set argument `overwrite` to `TRUE`."
    )
    stop(paste(e, collapse = " "), call. = FALSE)
  }

  if (!is_scalar_character(tabvar)) {
    e <- "You need to specify a single variable that should be summed up."
    stop(e, call. = FALSE)
  }

  poss <- c("freqs", attributes(x)$numvars)

  if (!tabvar %in% poss) {
    stop("Invalid name detected in argument `tabvar`.", call. = FALSE)
  }

  keep <- c("ind", tabvar, "costs", "lbi", "ubi", "LPL", "UPL", "SPL")
  x[[3]] <- x[[3]][, keep, with = FALSE]
  lapply(x, function(y) {
    write.table(
      data.frame(y),
      path,
      append = TRUE,
      sep = " " ,
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )
  })
  message("File ", shQuote(path), " successfully written.")
  return(invisible(path))
}
