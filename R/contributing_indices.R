#' Compute contributing units to table cells
#'
#' This function computes (with respect to the raw input data) the indices of all
#' contributing units to given cells identified by `ids`.
#'
#' @param prob a [sdcProblem-class] object created with [makeProblem()]
#' @param ids a character vector containing default ids (strIDs) that define table
#' cells. Valid inputs can be extracted by using [sdcProb2df()] and looking at
#' column `strID`. If this argument is `NULL`, the correspondig units are computed
#' for all cells in the table.
#'
#' @return a named `list where names correspond to the given `ids` and the values
#' to the row numbers within the raw input data.
#' @export
#' @md
#' @examples
#' # loading test data
#' data("microData1", package="sdcTable")
#'
#' # specify hierarchies for `age` and `region`
#' dim_region <- hier_create(root = "Total", nodes = LETTERS[1:4])
#' dim_gender <- hier_create(root = "Total", nodes = c("male", "female"))
#' dl <- list(region = dim_region, gender = dim_gender)
#'
#' # no variables holding counts, numeric values, weights or sampling
#' # weights are available in the input data
#'
#' # using variable names is also possible
#' prob <- makeProblem(
#'   data = microData1,
#'   dimList = dl
#' )
#'
#' df <- sdcProb2df(prob, dimCodes = "original")
#'
#' # which units contribute to cell region = "A" and gender = "female"?
#'
#' # compute the id ("0101")
#' df[region == "A" & gender == "female", strID]
#'
#' # which indices contribute to the cell?
#' ids <- contributing_indices(prob = prob, ids = "0101")
#'
#' # check
#' dataObj <- get.sdcProblem(prob, "dataObj")
#' rawData <- slot(dataObj, "rawData")
#' rawData[ids[["0101"]]]
#'
#' # compute contributing ids for each cell
#' contributing_indices(prob)
#'
contributing_indices <- function(prob, ids = NULL) {
  pi <- slot(prob, "problemInstance")
  poss_ids <- slot(pi, "strID")

  if (is.null(ids)) {
    ids <- poss_ids
  } else {
    if (!is.character(ids)) {
      stop("Please provide a character vector in argument `ids`.", call. = FALSE)
    }
    if (!all(ids %in% poss_ids)) {
      e <- c(
        "Some values provided in `ids` are not valid. ",
        "See column `strID` in `sdcProb2df()` for valid ids."
      )
      stop(paste(e, collapse = " "), call. = FALSE)
    }
  }

  contr_indices <- lapply(1:length(ids), function(x) {
    c_contributing_indices(
      object = prob,
      input = list(ids[x]))
  })
  names(contr_indices) <- ids
  contr_indices
}
