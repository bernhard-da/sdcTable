# compute inputs (cell total and top two contributing) units
# for a given vector of cell values and weights `w`
# values are replicated times their weights which are randomly
# rounded to integers in case they are floating numbers

# used in c_rule_p() and c_rule_nk()
.comp_weighted_inputs <- function(vals, w, n) {
  # replicate by weights: what to do with non-integerish weights?
  # randomly round upwards and downwards?
  if (!is_integerish(w)) {
    dir <- sample(c(-1, 1), length(vals), replace = TRUE)
    w[dir == -1] <- floor(w[dir == -11])
    w[dir == 1] <- ceiling(w[dir == 1])
  }
  vals <- rep(vals, times = w)

  top_contr <- rep(NA, n)
  v <- rev(tail(sort(vals), n))
  top_contr[1:length(v)] <- v
  list(
    cell_tot = sum(vals),
    top_contr = top_contr
  )
}
