prepare_long_data <- function(data, namvar, occ, ID, covariates = character(),
                              min_categories = 2L, consecutive = FALSE) {
  if (!is.data.frame(data) || !nrow(data) || !is.character(namvar) || !length(namvar)) {
    stop("Supply a nonempty data frame and the names of the item columns.", call. = FALSE)
  }
  key <- c(ID, occ, namvar)
  if (length(ID) != 1L || length(occ) != 1L || anyNA(key) ||
      anyDuplicated(key) || !all(c(key, covariates) %in% names(data))) {
    stop("Supply distinct names for the ID, occasion and item columns, and existing covariates.", call. = FALSE)
  }
  if (anyNA(data[[ID]]) || anyNA(data[[occ]]) ||
      (is.numeric(data[[occ]]) && any(!is.finite(data[[occ]])))) {
    stop("IDs and occasions must be observed and finite.", call. = FALSE)
  }
  if (anyDuplicated(data[c(ID, occ)])) {
    stop("Use one row per subject and occasion.", call. = FALSE)
  }

  dat <- as.data.frame(data[, unique(c(key, covariates)), drop = FALSE])
  times <- sort(unique(dat[[occ]]))
  subjects <- sort(unique(dat[[ID]]))
  row_subject <- match(dat[[ID]], subjects)
  row_time <- match(dat[[occ]], times)
  index <- order(row_subject, row_time)
  dat <- dat[index, , drop = FALSE]
  row_subject <- row_subject[index]
  row_time <- row_time[index]
  dat$..t <- row_time
  nsub <- length(subjects)
  nvar <- length(namvar)

  nmes <- matrix(0L, nsub, nvar)
  labels <- vector("list", nvar)
  y <- matrix(NA_integer_, nrow(dat), nvar)
  for (k in seq_len(nvar)) {
    item <- namvar[k]
    x <- dat[[item]]
    obs <- !is.na(x)
    if (!any(obs)) stop("Item has no observed responses: ", item, call. = FALSE)
    if (is.ordered(x)) {
      labels[[k]] <- levels(x)
      y[, k] <- as.integer(x)
    } else {
      v <- x[obs]
      if (!is.numeric(v) || any(!is.finite(v)) || any(v != floor(v)) ||
          !min(v) %in% c(0L, 1L) || max(v) >= .Machine$integer.max) {
        stop("Use an ordered factor or integer categories starting at 0 or 1 for item ", item, ".",
             call. = FALSE)
      }
      labels[[k]] <- seq.int(min(v), max(v))
      if (consecutive && length(unique(v)) != length(labels[[k]])) {
        stop("Use consecutive numeric categories or an ordered factor for item ", item, ".", call. = FALSE)
      }
      y[obs, k] <- as.integer(v + (min(v) == 0L))
    }
    if (length(labels[[k]]) < min_categories) {
      stop("Item ", item, " requires at least ", min_categories, " categories.", call. = FALSE)
    }
    nmes[, k] <- tabulate(row_subject[obs], nbins = nsub)
  }

  at <- which(!is.na(y), arr.ind = TRUE)
  row <- at[, 1L]
  packed <- order(row_subject[row], at[, 2L], row_time[row])
  list(data = dat, times = times, subjects = subjects,
       nmes = nmes, y = y[at[packed, , drop = FALSE]], time = row_time[row[packed]],
       categories = setNames(labels, namvar), c = lengths(labels),
       informative = sum(rowSums(nmes) > 0L))
}

loading_mask <- function(namvar, Q, simple, loadmask) {
  nvar <- length(namvar)
  if (is.null(Q)) {
    Q <- if (!is.null(loadmask)) ncol(as.matrix(loadmask)) else if (!is.null(simple)) max(simple) else 1L
  }
  Q <- check_count(Q, "Q")
  if (is.null(loadmask)) {
    simple <- if (is.null(simple)) rep(1L, nvar) else simple
    if (length(simple) != nvar || anyNA(simple) || any(!simple %in% seq_len(Q))) {
      stop("simple must give one factor index per item.", call. = FALSE)
    }
    loadmask <- matrix(0L, nvar, Q)
    loadmask[cbind(seq_len(nvar), simple)] <- 1L
  } else {
    loadmask <- as.matrix(loadmask)
    if (!identical(dim(loadmask), c(nvar, Q)) || anyNA(loadmask) ||
        any(!loadmask %in% c(0L, 1L))) {
      stop("loadmask must be a 0/1 matrix with one row per item and one column per factor.",
           call. = FALSE)
    }
    storage.mode(loadmask) <- "integer"
  }
  if (any(rowSums(loadmask) == 0L) || any(colSums(loadmask) == 0L)) {
    stop("Each item and factor must have an active loading in loadmask.", call. = FALSE)
  }
  dimnames(loadmask) <- list(namvar, paste0("F", seq_len(Q)))
  loadmask
}
