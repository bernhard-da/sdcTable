domRule <- function(object, params, type) {
  p_rule <- function(params) {
    stopifnot(is_scalar_double(params$cell_tot))

    cont1 <- params$top_contr[1]
    cont2 <- params$top_contr[2]
    stopifnot(is_scalar_double(cont1))
    stopifnot(is_scalar_double(cont2))
    stopifnot(is_scalar_integerish(params$p))
    stopifnot(params$p > 0, params$p <= 100)

    # if TRUE, cell needs to be suppressed
    (params$cell_tot - sum(cont1, cont2)) < (params$p / 100 * cont1)
  }
  pq_rule <- function(params) {
    stopifnot(is_scalar_double(params$cell_tot))
    cont1 <- params$top_contr[1]
    cont2 <- params$top_contr[2]
    stopifnot(is_scalar_double(cont1))
    stopifnot(is_scalar_double(cont2))
    stopifnot(is_integerish(params$pq))
    p <- params$pq[1]
    q <- params$pq[2]

    # if TRUE, cell needs to be suppressed
    (params$cell_tot - sum(cont1, cont2)) < (p / q) * cont1
  }
  nk_rule <- function(params) {
    stopifnot(is_scalar_double(params$cell_tot))
    stopifnot(is_scalar_integerish(params$n))

    # if TRUE, cell needs to be suppressed
    sum(params$top_contr) > (params$k / 100 * params$cell_tot)
  }

  if (type == "p") {
   fun <- p_rule
  }
  if (type == "pq") {
    fun <- pq_rule
  }
  if (type == "nk") {
    fun <- nk_rule
  }

  # compute inputs (cell total and top two contributing) units
  # for a given vector of cell values and weights `w`
  # values are replicated times their weights which are randomly
  # rounded to integers in case they are floating numbers
  .comp_weighted_inputs <- function(vals, w, n) {
    if (length(vals) == 0) {
      return(NULL)
    }
    # replicate by weights: what to do with non-integerish weights?
    # randomly round upwards and downwards?
    if (!is_integerish(w)) {
      dir <- sample(c(-1, 1), length(vals), replace = TRUE)
      w[dir == -1] <- floor(w[dir == -11])
      w[dir == 1] <- ceiling(w[dir == 1])
    }
    vals <- rep(vals, times = w)

    top_contr <- rep(0, n)
    v <- rev(tail(sort(vals), n))
    top_contr[1:length(v)] <- v
    list(
      cell_tot = sum(vals),
      top_contr = top_contr
    )
  }

  if (!is_scalar_character(type)) {
    stop("`type` needs to be a scalar character.", call. = FALSE)
  }
  if (!type %in% c("pq", "p", "nk")) {
    stop("`type` needs to be either `pq`, `p` or `nk`.", call. = FALSE)
  }

  if (!g_is_microdata(g_dataObj(object))) {
    e <- "nk-dominance rule can only be applied if micro-data are available!"
    stop(e, call. = FALSE)
  }

  if (type %in% c("p", "pq")) {
    n <- 2
  } else {
    n <- params$n
    if (n < 2) {
      stop("Parameter `n` must be >= 2 for nk-dominance rule.", call. = FALSE)
    }
  }

  pI <- g_problemInstance(object)
  dataObj <- g_dataObj(object)
  numVarInds <- g_numvar_ind(dataObj)
  strIDs <- g_strID(pI)
  numVal <- g_raw_data(dataObj)[[numVarInds[params$numVarInd]]]

  if (any(na.omit(numVal) < 0)) {
    e <- c(
      "dominance rules can only be applied to numeric variables",
      "with only positive values!"
    )
    stop(paste(e, collapse = " "), call. = FALSE)
  }

  # sampweights
  samp_weights <- g_sampweight_ind(dataObj)
  if (!is.null(samp_weights)) {
    samp_weights <- slot(dataObj, "rawData")[[samp_weights]]
  } else {
    samp_weights <- rep(1, length(numVal))
  }

  # calculate contributing indices
  nr_cells <- g_nrVars(pI)
  indices <- lapply(1:nr_cells, function(x) {
    c_contributing_indices(object, input = list(strIDs[x]))
  })

  # values and totals of contributing units
  inp <- lapply(1:nr_cells, function(x) {
    ii <- indices[[x]]
    .comp_weighted_inputs(
      vals = numVal[ii],
      w = samp_weights[ii],
      n = n
    )
  })

  # suppStatus: TRUE: unsafe, FALSE: safe
  supp_state <- sapply(1:nr_cells, function(x) {
    inp <- inp[[x]]
    # cells with value 0!
    if (is.null(inp)) {
      FALSE
    } else {
      fun(params = append(params, inp))
    }
  })

  supp_index <- which(supp_state == TRUE)
  if (length(supp_index) > 0) {
    s_sdcStatus(pI) <- list(
      index = supp_index,
      vals = rep("u", length(supp_index))
    )
  }

  if (isFALSE(params$allowZeros)) {
    if (type %in% c("p", "pq")) {
      ind_zero <- which(g_freq(pI) == 0)
    }
    if (type == "nk") {
      cell_totals <- unlist(lapply(inp, function(x) x$celltot))
      ind_zero <- which(cell_totals == 0 & g_freq(pI) >= 0)
    }
    if (length(ind_zero) > 0) {
      s_sdcStatus(pI) <- list(
        index = ind_zero,
        vals = rep("u", length(ind_zero))
      )
    }
  }
  s_problemInstance(object) <- pI
  validObject(object)
  return(object)
}

