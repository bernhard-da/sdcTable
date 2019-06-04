#' Create input for jj_format
#'
#' This function transforms a [sdcProblem-class] object into a list that can
#' be used as input for [writeJJFormat()] to write a problem in `"JJ-format"` to
#' disk.
#'
#' @param x a [sdcProblem-class] object
#' @return an input suitable for [writeJJFormat()]
#' @author Bernhard Meindl (bernhard.meindl@@statistik.gv.at) and
#' Sapphire Yu Han (y.han@@cbs.nl)
#' @export
#' @md
#' @inherit writeJJFormat examples
createJJFormat <- function(x) {
  i <- j <- v <- NULL

  if (!inherits(x, "sdcProblem")) {
    stop("`x` must be a `sdcProblem-object`.", call. = FALSE)
  }

  # extracts the "problemInstance"
  pi <- get.sdcProblem(x, type = "problemInstance")

  nr_cells <- get.problemInstance(pi, "nrVars")
  numvars <- as.data.table(get.problemInstance(pi, "numVars"))

  # prepare output
  jj <- vector("list", 5)

  jj[[1]] <- 0 # first element must be 0
  jj[[2]] <- nr_cells # number of cells

  # matrix part containing data values and bounds
  if (ncol(numvars) > 0) {
    m <- data.table(ind = (1:nr_cells) - 1, numvars)
  } else {
    m <- data.table(ind = (1:nr_cells) - 1)
  }

  m$freqs <- get.problemInstance(pi, "freq")
  m$costs <- get.problemInstance(pi, "weight")
  m$status <- get.problemInstance(pi, "sdcStatus")
  m$lbi <- get.problemInstance(pi, "lb")
  m$ubi <- get.problemInstance(pi, "ub")
  m$LPL <- get.problemInstance(pi, "LPL")
  m$UPL <- get.problemInstance(pi, "UPL")
  m$SPL <- get.problemInstance(pi, "SPL")
  jj[[3]] <- m

  # number of linear dependencies
  # all deps as sparse matrix
  st <- c_gen_mat_m(
    input = list(
      objectA = pi,
      objectB = get.sdcProblem(x, type = "dimInfo")
    )
  )

  nr_constraints <- slot(st, "nrRows")
  jj[[4]] <- nr_constraints

  # the fifth element of the list are the
  # linear dependences as vector
  mm <- data.table(v1 = rep("0.0", nr_constraints))
  mm$v2 <- as.character(table(slot(st, "i")))
  mm$v3 <- rep(":", nr_constraints)

  dt <- data.table(
    i = slot(st, "i"),
    j = slot(st, "j") - 1,
    v = slot(st, "v")
  )
  setkey(dt, i)
  con <- dt[, paste(j, "(", v, ")", collapse = " ") , by = key(dt)]

  mm$v4 <- con[["V1"]]
  jj[[5]] <- mm

  setattr(jj, "numvars", names(numvars))
  class(jj) <- "jjformat"
  return(invisible(jj))
}
