#' Validate Event Table
#' 
#' Tests whether the object meets the basic requirements of an event table, i.e. a data frame containing at least two numeric, finite columns named 'from' and 'to' ordered such that \code{to} > or = \code{from} on all rows.
#' 
#' @param x an R object.
#' @param verbose logical value indicating whether to print the reason for test failure.
#' @seealso \code{\link{events}}, \code{\link{as_events}}, and \code{\link{read_events}} for creating valid event tables.
#' @export
#' @examples
#' verbose <- TRUE
#' is_events(c(1, 3), verbose)
#' is_events(data.frame(from = 1, t = 3), verbose)
#' is_events(data.frame(from = 1, from = 1, to = 3), verbose)
#' is_events(data.frame(from = 1, to = TRUE), verbose)
#' is_events(data.frame(from = 1, to = NA), verbose)
#' is_events(data.frame(from = 3, to = 1), verbose)
#' is_events(data.frame(from = 1, to = 3), verbose)   # TRUE
is_events <- function(x, verbose = FALSE) {
  if (!is.data.frame(x)) {
    if (verbose) cat("x is not a data.frame\n")
    return(FALSE)
  }
  occurrence <- lapply(rgrep_exact(c("from", "to"), names(x)), length)
  if (any(occurrence == 0)) {
    if (verbose) cat("columns from and to are missing\n")
    return(FALSE)
  }
  if (any(occurrence > 1)) {
    if (verbose) cat("columns from and to appear more than once\n")
    return(FALSE)
  }
  if (!is.numeric(x$from) || !is.numeric(x$to)) {
    if (verbose) cat("columns from and to are not numeric\n")
    return(FALSE)
  }
  if (any(!is.finite(c(x$from, x$to)))) {
    if (verbose) cat("columns from and to contain non-finite values (i.e. NA, NaN, Inf)\n")
    return(FALSE)
  }
  if (any(x$from > x$to)) {
    if (verbose) cat("from > to on one or more rows\n")
    return(FALSE)
  }
  return(TRUE)
}

#' Validate Bin Event Table
#' 
#' Tests whether an object is an event table, plus whether the event table contains only groups of non-overlapping line events. FIXME: Deprecated, since only one function still actually requires this (\code{\link{seq_bin_events}}).
#' 
#' @param x an R object.
#' @param group.col name or numeric index of the column defining the groupping of bins. If \code{NULL}, the bins are treated as one group. 
#' @param verbose logical value indicating whether to print the reason for test failure.
#' @seealso \code{\link{seq_bin_events}} and \code{\link{group_nonoverlapping_events}} for creating valid bin event tables.
#' @keywords internal
is_bin_events <- function(x, group.col = NULL, verbose = FALSE) {
  if (!is_events(x, verbose = verbose)) {
    if (verbose) cat("x is not an event table\n")
    return(FALSE)
  }
  if (any(x$from == x$to)) {
    if (verbose) cat("x contains point events\n")
    return(FALSE)
  }
  if (!is.null(group.col)) {
    if (is.character(group.col))
      group.col <- match(group.col, names(x))
    overlaps <- unlist(lapply(split(x, x[group.col]), has_overlapping_events))
  } else {
    overlaps <- has_overlapping_events(x)
  }
  if (any(overlaps)) {
    if (verbose) cat("x contains overlapping events\n")
    return(FALSE)
  }
  return(TRUE)
}

#' Sorted Events
#' 
#' \code{sort_events} sorts events by ascending \code{from}, then ascending \code{to}. \code{is_unsorted_events} tests whether the events are not sorted, without the cost of sorting them.
#' 
#' @param e an event table
#' @export
#' @rdname sorted_events
#' @examples
#' e <- events(c(1, 1, 3, 2), c(2, 1, 4, 3))
#' is_unsorted_events(e)
#' sort_events(e)
sort_events <- function(e) {  
  ind <- order(e$from, e$to)
  return(e[ind, , drop = FALSE])
}
#' @export
#' @rdname sorted_events
is_unsorted_events <- function(e) {
  n <- nrow(e)
  if (n < 2)
    return(FALSE)
  else
    return(any(is.unsorted(e$from), any(diff(e$from) == 0 & diff(e$to) < 0)))
}

#' Overlapping Events
#' 
#' \code{group_nonoverlapping_events} assigns each event to a group such that each group contains no overlaps. \code{has_overlapping_events} checks whether an event table has overlapping events.
#' 
#' By convention in \code{linbin}, events are considered overlapping if they are line events that share more than an endpoint, or point events that have equal endpoints. Point events on line event endpoints are not considered overlaps.
#' 
#' @param e an event table.
#' @export
#' @keywords internal
#' @rdname overlapping_events
#' @examples
#' e <- events(c(0, 2, 3), c(3, 4, 5))
#' cbind(group = group_nonoverlapping_events(e), e)  # adjacent lines do not overlap
#' e <- events(c(0, 0, 0, 1, 1), c(0, 0, 1, 1, 2))    
#' cbind(group = group_nonoverlapping_events(e), e)  # equal points do overlap
#' has_overlapping_events(events(c(0, 2), c(2, 4)))  # adjacent lines
#' has_overlapping_events(events(c(0, 2), c(3, 4)))  # has overlapping lines
#' has_overlapping_events(events(c(0, 5, 5, 10)))    # points adjcent to lines
#' has_overlapping_events(events(c(0, 5, 5, 5, 10))) # has overlapping points
group_nonoverlapping_events = function(e) {
  # Sort bins as needed
  if (is_unsorted_events(e)) {
    ids <- order(e$from, e$to)
    e <- e[ids, c("from", "to")]
    reorder <- TRUE
  } else {
    reorder <- FALSE
  }
  # Loop through bins, assigning each to a group
  N <- nrow(e)
  groups <- numeric(N)
  s <- 1
  i <- 1
  n <- 0
  repeat {
    # assign bin to set
    groups[i] <- s
    n <- n + 1
    if (n == N)
      # all bins assigned
      break
    k <- i + 1
    # move forward to nearest, unassigned, non-overlapping bin
    while ((e$from[k] < e$to[i] || groups[k] > 0 || (e$from[k] == e$from[i] && e$to[k] == e$to[i])) && k <= N)
      k <- k + 1
    if (k > N) {
      # start new set
      s <- s + 1 
      # back to top
      i <- match(0, groups)
    } else {
      i <- k
    }
  }
  # Reorder as needed
  if (reorder) {
    return(groups[order(ids)])
  } else {
    return(groups)
  }
}
#' @export
#' @rdname overlapping_events
has_overlapping_events <- function(e) {
  ne <- nrow(e)
  if (ne < 2)
    return(FALSE)
  if (is_unsorted_events(e))
    events <- sort_events(e)
  overlines <- e$from[-1] < e$to[-ne]
  overpoints <- e$from[-1] == e$from[-ne] & e$to[-1] == e$to[-ne]
  return(any(overlines | overpoints))
}