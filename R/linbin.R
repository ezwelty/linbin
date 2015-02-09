# === EVENT CREATION ===================================================================================================

#' Event Tables
#' 
#' This function creates an event table, a custom \code{data.frame} used throughout the \code{linbin} package to store and manipulate linearly referenced data. Each row includes an event's endpoints \code{from} and \code{to} (which can be equal, to describe a point, or non-equal, to describe a line) and the values of any variables measured on that interval.
#' 
#' Event endpoints (and any additional arguments) are coerced to a data frame by calling \code{\link{data.frame}}, then coerced to an event table by calling \code{\link{as_events}}. A valid event table has two columns named "from" and "to" containing only finite numeric values (i.e., no \code{NA}, \code{NaN}, or \code{Inf}) and ordered such that \code{to} ≥ \code{from}. \code{\link{is_events}} tests for these requirements. The other columns in the event table can be of any type supported by the \code{\link{data.frame}} class. The core routines of \code{linbin} — particularly \code{\link{sample_events}} and \code{\link{plot_events}} — rely on "bins", event tables consisting of only non-overlapping, line events (\code{\link{is_bin_events}}) to summarize and plot event data.
#' 
#' @param from,to event endpoints, in any format coercible to single data frame columns. \code{from} and \code{to} are swapped as needed so that \code{to} ≥ \code{from} on all rows. If only \code{from} is passed, \code{\link{as_events}} is dispatched for object coercion.
#' @param ... additional arguments, which can be either of the form \code{value} or \code{tag = value}, to be  passed directly to \code{\link{data.frame}} following \code{from} and \code{to}. Component names are created based on the tag (if present) or the deparsed argument itself.
#' @return An event table, the \code{data.frame} object used by \code{linbin} to describe interval data.
#' @seealso \code{\link{data.frame}}.
#' @seealso \code{\link{as_events}} and \code{\link{read_events}} for coercing objects and files to event tables, \code{\link{is_events}} to validate event tables, and \code{\link{seq_bin_events}} and \code{\link{is_bin_events}} to generate and validate bin event tables.
#' @family event creation
#' @export
#' @examples
#' events(1, 10)
#' events(1:10)
#' events(c(0, 15, 25), c(10, 30, 35), x = 1, f = c('a', 'b', 'c'))
events = function(from = integer(), to = integer(), ...) {
  if (!length(to) && !length(list(...)) && length(from))
    return(as_events(from))
  x = data.frame(from, to, ...)
  return(as_events(x))
}

#' Coerce to an Event Table
#' 
#' Attempts to coerce an object to an event table.
#' 
#' @param x object to be coerced to an event table.
#' @param from.col,to.col names or indices of the columns in \code{x} containing the event endpoints. Values are swapped as needed to ensure that \code{to ≥ from} on all rows.
#' @return An event table, the \code{data.frame} object used by \code{linbin} to describe interval data.
#' @seealso \code{\link{events}}, \code{\link{read_events}} for reading files as event tables.
#' @family event creation
#' @export
#' @examples
#' as_events(1)
#' as_events(1:10)
#' as_events(cbind(1:10, 1:10), from.col = 1, to.col = 2)
#' as_events(data.frame(from = 1:10, to = 1:10))
as_events = function(x, ...) {
  if (is.null(x)) 
    return(data.frame(from = integer(), to = integer()))
  UseMethod("as_events")
}
#' @describeIn as_events Expands a numeric vector into two columns of event endpoints: \code{(x1, x2), (x2, x3), (x3, x4)}.
as_events.numeric = function(x) {
  len = length(x)
  if (len > 1) {
    e = data.frame(from = x[-len], to = x[-1])
    need.flip = e$from > e$to
    if (any(need.flip))
      e[need.flip, c("from", "to")] = e[need.flip, c("to", "from")]
    return(e)
  } else {
    return(data.frame(from = x, to = x))
  }
}
#' @describeIn as_events Converts the matrix to a data frame, then calls the \code{data.frame} method.
as_events.matrix = function(x, from.col = "from", to.col = "to") {
  as_events(as.data.frame(x, check.names = TRUE), from.col, to.col)
}
#' @describeIn as_events Renames \code{from.col} and \code{to.col} to "from" and "to" as needed. Since these column names must be unique, other columns cannot also be called "from" or "to".
as_events.data.frame = function(x, from.col = "from", to.col = "to") {
  # Check presence and uniqueness of from and to columns
  if (is.character(from.col)) from.col = which(names(x) %in% from.col)
  if (is.character(to.col)) to.col = which(names(x) %in% to.col)
  names(x)[c(from.col[1], to.col[1])] = c("from", "to")
  occurrence = lapply(rgrep.exact(c("from", "to"), names(x)), length)
  if (any(occurrence < 1))
    stop("One or both of the required columns (from.col and to.col) are missing")
  if (any(occurrence > 1))
    stop('One or both of the reserved column names ("from" and "to") appear more than once')
  # Coerce to numeric
  x$from = as.numeric(x$from)
  x$to = as.numeric(x$to)
  # Check if finite
  if (!all(is.finite(x$from)) || !all(is.finite(x$to)))
    stop("from.col and to.col cannot contain non-finite values (NA, NaN, and Inf)")
  # Order events from -> to
  need.flip = x$from > x$to
  if (any(need.flip))
    x[need.flip, c("from", "to")] = x[need.flip, c("to", "from")]
  return(x)
}

#' Read File as Event Table
#' 
#' Reads a file in table format and attempts to coerce it to an event table.
#' 
#' The file is read into R by calling \code{\link{read.table}}. Any of its arguments can be set by passing additional \code{tag = value} pairs to \code{read_events}. \code{from.col} and \code{to.col} are renamed to "from" and "to" as needed. Since these column names must be unique, other columns cannot also be called "from" or "to".
#' 
#' @param file name, \code{\link{connection}}, or \code{\link{url}} of the file to be read as an event table.
#' @param from.col,to.col names or indices of the columns containing event endpoints. Values are swapped as needed to ensure that \code{to ≥ from} on all rows.
#' @param header logical value indicating whether the file contains column names as its first line. If \code{FALSE}, columns will be named "V" followed by the column number, unless \code{col.names} (a vector of optional column names) is provided as an additional argument.
#' @param sep character seperating values on each line of the file. If \code{sep = ""} (the default), the separator is 'white space' (that is, any combination of one or more spaces, tabs, newlines and carriage returns).
#' @param ... additional arguments, of the form \code{tag = value}, to be passed directly to \code{\link{read.table}} to control how the file is read.
#' @return An event table, the \code{data.frame} object used by \code{linbin} to describe interval data.
#' @seealso \code{\link{read.table}}.
#' @seealso \code{\link{events}}, \code{\link{as_events}} for coercing existing objects.
#' @family event creation
#' @export
read_events = function(file, from.col = "from", to.col = "to", sep = "", header = TRUE, ...) {
  x = read.table(file, sep = sep, header = header, ...)
  return(as_events(x, from.col, to.col))
}

#' Tests for Event Table Types
#' 
#' \code{is_events} tests whether the object meets the requirements of an event table. \code{is_bin_events} tests all these, plus whether the event table contains only groups of non-overlapping line events.
#' 
#' @param x any R object.
#' @param groups a vector defining the groupping of bins. If \code{NULL}, the bins are treated as one group. 
#' @param verbose logical value indicating whether to print the reason for test failure.
#' @return \code{TRUE} if the argument meets the criteria of the test, and \code{FALSE} otherwise.
#' @seealso \code{\link{events}}, \code{\link{as_events}}, \code{\link{group_nonoverlapping_events}}.
#' @family event creation
#' @export is_events
#' @rdname is_events
is_events = function(x, verbose = FALSE) {
  if (!is.data.frame(x)) {
    if (verbose) warning(paste0("x is not a data.frame"))
    return(FALSE)
  }
  occurrence = lapply(rgrep.exact(c("from", "to"), names(x)), length)
  if (any(occurrence == 0)) {
    if (verbose) warning("columns from and to are missing")
    return(FALSE)
  }
  if (any(occurrence > 1)) {
    if (verbose) warning("columns from and to appear more than once")
    return(FALSE)
  }
  if (!is.numeric(x$from) || !is.numeric(x$to)) {
    if (verbose) warning("columns from and to are not numeric")
    return(FALSE)
  }
  if (any(!is.finite(c(x$from, x$to)))) {
    if (verbose) warning("columns from and to contain non-finite values (i.e. NA, NaN, Inf)")
    return(FALSE)
  }
  if (any(x$from > x$to)) {
    if (verbose) warning(paste("FROM > TO on one or more rows."))
    return(FALSE)
  }
  return(TRUE)
}
#' @family bin event creation
#' @export is_bin_events
#' @rdname is_events
is_bin_events = function(x, groups = NULL, verbose = FALSE) {
  if (!is_events(x, verbose = verbose)) {
    if (verbose) warning("x is not an event.table")
    return(FALSE)
  }
  if (any(x$from == x$to)) {
    if (verbose) warning(paste("x contains point events"))
    return(FALSE)
  }
  if (!is.null(group))
    overlaps = unlist(lapply(split(x, groups), has.overlaps))
  else
    overlaps = has.overlaps(x)
  if (any(overlaps)) {
    if (verbose) warning("x contains overlapping events")
    return(FALSE)
  }
  return(TRUE)
}

# === EVENT METRICS ===================================================================================================

#' Event Range
#' 
#' Returns the minimum and maximum endpoints of all the events in an event table.
#' 
#' @param e an event table.
#' @return An endpoint-only event table with a single row.
#' @seealso \code{\link{event_coverage}}
#' @family event metrics
#' @export
event_range = function(e) {
  return(as_events(range(e$from, e$to)))
}

#' Event Coverage
#' 
#' Returns the intervals over which the number of events is always one or greater.
#' 
#' FIXME: Currently drops any isolated points on the boundaries of the event range, and drops all isolated points when \code{closed = FALSE}. What should be the behavior for points?
#'
#' @param e an event table.
#' @param closed logical value indicating whether events should be interpreted as closed intervals. If \code{TRUE}, coverage is continuous at breaks between two adjacent events (the default).
#' @return An endpoint-only event table.
#' @seealso \code{\link{event_gaps}} for gaps (the inverse of coverage), \code{\link{event_range}} for range (coverage with gaps ignored).
#' @family event metrics
#' @export
event_coverage = function(e, closed = TRUE) {
  e.gaps = event_gaps(e, closed = closed)
  return(event_gaps(e.gaps, bounds = event_range(e), closed = closed))
}

#' Event Gaps
#' 
#' Returns the intervals over which there are no events.
#'
#' @param e an event table.
#' @param closed logical value indicating whether events should be interpreted as closed intervals. If \code{TRUE}, no gaps are returned at breaks between two adjacent events (the default).
#' @param range an event table specifying, by its \code{\link{event_range}}, the interval within which to check for gaps.
#' @return An endpoint-only event table.
#' @seealso \code{\link{event_coverage}} for coverage (the inverse of gaps), \code{\link{event_range}} for range (coverage with gaps ignored).
#' @family event metrics
#' @export
event_gaps = function(e, closed = TRUE, range = NULL) {
  # Crop to bounds if not event range
  if (!is.null(range)) {
    range = event_range(range)
    if (!nrow(e))
      # no data in bounds
      return(bounds)
    e.range = event_range(e)
    if (!(range$from <= e.range$from && range$to >= e.range$to))
      # bounds intersect event range
      e = crop_events(e, range)
    if (range$from < e.range$from)
      # bounds extend past event range (from)
      e = rbind(rep(e.range$from, 2), e[c("from", "to")])
    if (range$to > e.range$to)
      # bounds extend past event range (to)
      e = rbind(e[c("from", "to")], rep(range$to, 2))
  }
  # Track overlaps by extending events by cumulative max
  if (is_unsorted_events(e))
    e = sort_events(e)
  e$to = cummax(e$to)
  # Gaps occur when from[i+1] > to[i] if closed intervals
  # (and when == if open intervals)
  if (closed) { 
    isgap = which(e$from[-1] > e$to[-nrow(x)])
  } else { 
    isgap = which(e$from[-1] >= e$to[-nrow(x)])
  }
  # Build gaps event table
  return(events(e$to[isgap], e$from[isgap + 1]))
}

# === EVENT OPERATIONS ===================================================================================================

#' Fill Event Gaps
#' 
#' Fills gaps below a maximum length with events with empty (NA) variables.
#' 
#' @param e an event table.
#' @param max.length the maximum length of gaps to be filled.
#' @return An event table with new rows appended to the bottom.
#' @seealso \code{\link{event_gaps}}
#' @family event operations
#' @export
fill_event_gaps = function(e, max.length = Inf) {
  e.gaps = event_gaps(e)
  gap.lengths = e.gaps$to - e.gaps$from
  fill = gap.lengths < max.length
  ngaps = sum(fill)
  nevents = nrow(e)
  if (ngaps > 0) {
    e = e[c(seq_len(nrow(e)), rep(NA, ngaps)), , drop = FALSE]
    e[nevents + (1:ngaps), c("from", "to")] = gaps[fill, ]
  }
  return(e)
}

#' Sort Events
#' 
#' Sorts events by ascending \code{from}, then ascending \code{to}. The function \code{is_unsorted_events} checks whether the events are not sorted  without the cost of sorting them.
#' 
#' @param e an event table.
#' @return An event table with sorted rows.
#' @family event operations
#' @export
#' @rdname sort_events
sort_events = function(e) {  
  ind = order(e$from, e$to)
  return(e[ind, , drop = FALSE])
}
#' @return \code{TRUE} if the event table is not sorted, and \code{FALSE} otherwise.
#' @export
#' @rdname sort_events
is_unsorted_events = function(e) {
  n = nrow(e)
  if (n < 2)
    return(FALSE)
  else
    return(any(is.unsorted(e$from), any(diff(e$from) == 0 & diff(e$to) < 0)))
}

#' Transform Events
#' 
#' Transforms events by scaling, then translating their endpoint positions. That is, the transformed \code{[from, to] = scale * [from, to] + translate}.
#' 
#' @param e an event table.
#' @param scale number by which event endpoints should be scaled.
#' @param translate number by which event endpoints should be translated.
#' @return An event table with scaled and translated \code{from} and \code{to} columns.
#' @family event operations
#' @export
transform_events = function(e, scale = 1, translate = 0) {
  e[c("from","to")] = scale[1] * e[c("from","to")] + translate[1]
  return(e)
}

#' Crop Events
#' 
#' Crops events to the specified intervals. Events are cut at interval endpoints and any whole or partial events lying outside the intervals are removed.
#'
#' @param e an event table.
#' @param coverage an event table specifying, by its \code{\link{event_coverage}}, the intervals for cropping.
#' @param closed logical indicating whether the events should be cut at breaks between adjacent cropping intervals.
#' @param scaled.cols names or indices of the event table columns to be rescaled after cutting (see \code{\link{cut_events}}). Names are interpreted as regular expressions (\code{\link{regex}}) matching full column names.
#' @return An event table with rows added by cutting and removed by cropping.
#' @seealso \code{\link{cut_events}}.
#' @family event operations
#' @export
crop_events = function(e, coverage, closed = FALSE, scaled.cols = NULL) {  
  coverage = event_coverage(coverage, closed = closed)
  e.cut = cut_events(e, cuts = coverage, scaled.cols = scaled.cols)
  inx = find_intersecting_events(coverage, e.cut)
  keep = .rowSums(inx, nrow(inx), ncol(inx)) > 0
  return(e.cut[keep, , drop = FALSE])
}

#' Cut Events
#' 
#' Cuts events at the specified locations.
#'
#' Line events are cut into multiple events Columns \code{scaled.cols} are scaled by the fraction of the original event length in each resulting event (which assumes that these variables were uniformly distributed over the original interval). To record the parents of cut events, append an unique identification field to the event table before calling this function.
#'
#' @param e an event table.
#' @param cuts the cut locations. May be either a numeric vector or an event table.
#' @param scaled.cols names or indices of the event table columns to be scaled to their new length after cutting. Names are interpreted as regular expressions (\code{\link{regex}}) matching full column names.
#' @return An event table with rows added by cutting.
#' @seealso \code{\link{crop_events}}.
#' @family event operations
#' @export
cut_events = function(e, cuts, scaled.cols = NULL) {
  # Initialize inputs
  if (is.numeric(cuts))
    cuts = sort(unique(cuts))
  else
    cuts = sort(c(cuts$from, cuts$to))
  if (is.character(scaled.cols))
    scaled.cols = unique(unlist(rgrep_exact(scaled.cols, names(e))))
  if (!length(scaled.cols))
    scaled.cols = NULL
  # Find events straddling the cuts
  incuts = find_intersecting_events(events(cuts, cuts), e)
  ncuts = rowSums(incuts)
  # Expand events to accomodate new event segments
  ei = seq_len(nrow(e))
  e = e[rep(ei, ncuts + 1), ]
  # Indices of cut segments
  ind = ncuts > 0
  ncuts = ncuts[ind]
  i = ei[ind]
  i0 = i + c(0, cumsum(ncuts))[seq_along(i)]
  i1 = i0 + ncuts
  # For each cut event
  for (n in seq_along(i)) {
    # Update measures
    pts = c(e$from[i0[n]], cuts[incuts[i[n],]], e$to[i0[n]])
    e$from[i0[n]:i1[n]] = pts[-length(pts)]
    e$to[i0[n]:i1[n]] = pts[-1]
    # Rescale columns
    if (!is.null(scaled.cols)) {
      l = diff(pts)
      e[i0[n]:i1[n], scaled.cols] = e[i0[n]:i1[n], scaled.cols] * (l / sum(l))
    }
  }
  return(e)
}

# === BIN CREATION ===================================================================================================

#' Generate Sequential Bins
#' 
#' Generates groups of regularly sequenced events fitted to the specified intervals. Intended for use as bins with \code{\link{sample_events}}.
#'
#' @param bin.coverage a bin event table specifying the intervals for sequencing. Gaps in coverage do not count towards bin length.
#' @param bin.count number of bins, a numeric vector, from which the bin length is calculated.
#' @param bin.length length of bins, a numeric vector. Ignored if bin.count is defined. When bin.length does not evenly divide the coverage, a shorter bin is appended to the end of the sequence.
#' @param adaptive logical value indicating whether the lengths of bins should be adjusted locally so that a whole number of bins fit within each coverage interval, preserving breaks and gaps.
#' @return An endpoint-only bin event table. An additional group field is appended if the length of \code{bin.count} or \code{bin.length} is > 1.
#' @seealso \code{\link{event_range}}, \code{\link{event_coverage}}, and  \code{\link{crop_events}} for interval operations useful for building bins tailored to a particular situation.
#' @family bin event creation
#' @export
seq_bin_events = function(bin.coverage, bin.count = NULL, bin.length = NULL, adaptive = FALSE) {
  
  # Check inputs
  if (!is_bin_events(bin.coverage))
    stop('bin.coverage is not a bin event table')
  # Flatten event table, calculate total coverage
  total.length = sum(bin.coverage$to - bin.coverage$from)
  if (!is.null(bin.count))
    bin.length = total.length / round(as.numeric(bin.count))
  
  # Generate bins for each bin length
  seq.bins = lapply(bin.length, function(bin.length) {
    if (!adaptive) {
      # Build initial from and to values for the bins
      from = min(bin.coverage$from)
      to = from + total.length
      binseq = seq(from, to, bin.length)
      if (max(binseq) < to)
        # Add smaller bin to reach end of coverage
        binseq[length(binseq) + 1] = to
      # Reinject gaps into bins
      gaps = event_gaps(bin.coverage)
      if (nrow(gaps)) {
        # Locate start of gaps in bin sequence
        pos = findInterval(gaps$from, binseq, rightmost.closed = TRUE)
        temp = numeric(length(binseq))
        temp[pos + 1] = gaps$to - gaps$from
        binseq = binseq + cumsum(temp)
      }
      return(cbind(binseq[-length(binseq)], binseq[-1]))
    } else {
      # Fit bins to intervals of coverage
      seg.length = bin.coverage$to - bin.coverage$from
      # Find evenly dividing length closest to nominal length
      r = seg.length / bin.length
      # Require at least one bin
      r[r < 1] = 1
      l1 = seg.length / ceiling(r)
      l2 = seg.length / floor(r)
      d1 = abs(l1 - bin.length)
      d2 = abs(l2 - bin.length)
      smaller.d2 = d1 < d2
      l1[smaller.d2] = l2[smaller.d2]
      binmat = do.call(rbind, lapply(seq_len(nrow(bin.coverage)), 
             function(i) {
               binseq = seq(bin.coverage$from[i], bin.coverage$to[i], l1[i])
               cbind(binseq[-length(binseq)], binseq[-1])
             }))
      return(binmat)
    }
  })
  # Format to event.table with group id
  if (length(seq.bins) > 1) {
    seq.bins = Map(cbind, seq_along(seq.bins), seq.bins)
    seq.bins = as_events(do.call(rbind, seq.bins), 2, 3)
    names(seq.bins)[1] = "group"
  } else {
    seq.bins = as_events(do.call(rbind, seq.bins), 1, 2)
  }
  return(seq.bins)
}

#' Overlapping Events
#' 
#' \code{group_nonoverlapping_events} assigns each event to a group such that each group contains no overlaps. \code{has_overlapping_events} checks whether an event table has overlaps, without the cost of computing non-overlapping groups
#' 
#' Proceeding from the first event (after sorting), the next still-unassigned event that doesn't overlap is added to the group. When the last event is reached, a new group is created starting from the top-most unassigned event until all events have been assigned to a group. The routine is used internally by \code{\link{sample_events}} to iterate over groups of bins for efficiency.
#' 
#' FIXME: Intervals are assumed open, so point events with equal endpoints are not considered overlaps. Single point events adjacent to line events may be desirable in some cases, however, so maybe point events should only be closed with respect to each other? What is the appropriate behavior?
#' 
#' @param e an event table.
#' @return A numeric vector indicating the group assignment for each event.
#' @seealso \code{\link{seq_bin_events}}, \code{\link{sample_events}}.
#' @family bin event creation
#' @export
#' @keywords internal
#' @rdname overlapping_events
group_nonoverlapping_events = function(e) {
  # Sort bins as needed
  if (is_unsorted_events(e)) {
    ids = order(e$from, e$to)
    e = e[ids, c("from", "to")]
    reorder = TRUE
  } else {
    reorder = FALSE
  }
  # Loop through bins, assigning each to a group
  N = nrow(e)
  groups = numeric(N)
  s = 1
  i = 1
  n = 0
  repeat {
    # assign bin to set
    groups[i] = s
    n = n + 1
    if (n == N)
      # all bins assigned
      break
    k = i + 1
    # move forward to nearest, unassigned, non-overlapping bin
    while ((e$from[k] < e$to[i] || groups[k] > 0) && k <= N)
      k = k + 1
    if (k > N) {
      # start new set
      s = s + 1 
      # back to top
      i = match(0, groups)
    } else {
      i = k
    }
  }
  # Reorder as needed
  if (reorder) {
    return(groups[ids[ids]])
  } else {
    return(groups)
  }
}
#' @return \code{TRUE} if the event table contains overlapping events, and \code{FALSE} otherwise.
#' @export
#' @rdname overlapping_events
has_overlapping_events = function(e) {
  n = nrow(e)
  if (n < 2)
    return(FALSE)
  if (is_unsorted.events(e))
    events = sort_events(e)
  return(any(e$from[-1] < e$to[-n]))
}

# === SAMPLING ===================================================================================================

#' Find Intersecting Events
#' 
#' Returns a logical matrix indicating whether or not each pair of events intersect.
#' 
#' @param ex,ey event tables.
#' @param closed.from,closed.to logical values indicating whether to treat events \code{ex} as closed intervals at their \code{from} and \code{to} endpoints, respectively. If \code{FALSE}, events \code{ey} only sharing an endpoint with events in \code{ex} are not reported as intersecting (the default).
#' @return A logical matrix with \code{ey} events as rows and \code{ex} events as columns.
#' @seealso \code{\link{sample_events}}.
#' @family event sampling
#' @export
find_intersecting_events = function(ex, ey, closed.from = FALSE, closed.to = FALSE) {
  inbins = apply(ex[c("from", "to")], 1, function(ex.each) {
    if (closed.from) {
      inleft = ex.each[1] <= ey$to
    } else {
      inleft = ex.each[1] < ey$to
    }
    if (closed.to) {
      inright = ex.each[2] >= ey$from
    } else {
      inright = ex.each[2] > ey$from
    }
    inbin = inleft & inright
  })
  return(inbins)
  # Converting from logical to a row column 
  # (for single matches): apply(inbins, 1, function(x) match(TRUE, x))
  # (for multiple matches): which(inbins, arr.ind = TRUE)
}

#' Sample Events
#' 
#' Sample event table variables over the specified intervals.
#' 
#' Events are cut at bin endpoints, and any \code{scaled.cols} columns are rescaled to the length of the resulting event segments. The event segments falling into each bin are passed to the sampling functions to compute variables for each bin. 
#' 
#' Sampling functions are specified in lists with the format \code{list(FUN, data.cols, by = group.cols, ...)}. The first element in the list is the function to use. It must accept one or multiple vectors of the same length and output a single vector of the same length. The following unnamed element is a vector specifying the event column names or indices to recursively pass as the first argument of the function. Names are interpreted as regular expressions (\code{\link{regex}}) matching full column names. Additional unnamed elements are vectors specifying additional event columns to pass as the second, third, ... argument of the function. The first "by" element is a vector of event column names or indices used as grouping variables. Any additional named arguments are passed directly to the function. For example:
#' 
#' list(sum, 1:2, na.rm = TRUE) => sum(events[1], na.rm = TRUE), sum(events[2], na.rm = TRUE)
#' list(sum, 1, 3:4, 5) => sum(events[1], events[3], events[4], events[5]), ...
#' list(sum, c('x', 'y'), by = 3:4) => list(sum, 'x'), list(sum, 'y') grouped into all combinations of 3 and 4
#' 
#' Using the last example above, column names are taken from the first argument (e.g. \code{x, y}), and all grouping variables are appended (e.g. \code{x.a.a, x.a.b, x.b.a, y.b.b}), where \code{a} and \code{b} are the levels of columns 3 and 4. \code{NA} is also treated as a factor level. Columns are added left to right in order of the sampling function arguments. Finally, names are made unique by appending sequence numbers to duplicates (using \code{\link{make.unique}}).
#' 
#' FIXME: Sampling at points should be supported. By default, should sample from overlapping line events. Should it also sample from co-located point events and adjacent line events?
#' 
#' @param e an event table.
#' @param bins a bin event table.
#' @param ... lists specifying the sampling functions to be used (see the \code{Details}).
#' @param scaled.cols names or indices of the event columns to be rescaled after cutting (see \code{\link{cut_events}}). Names are interpreted as regular expressions (\code{\link{regex}}) matching full column names.
#' @param col.names a character vector of the column names overriding the automatic naming of data columns output by the sampling functions (see the \code{Details}).
#' @param drop.empty logical value indicating whether or not to drop empty bins. If \code{FALSE}, empty bins are retained with \code{NA} columns (the default).
#' @return A bin event table with dava columns appended.
#' @seealso \code{\link{seq_bin_events}} to generate sequential bins.
#' @family event sampling
#' @export
sample_events = function(e, bins, ..., scaled.cols = NULL, col.names = NULL, drop.empty = FALSE) {
  
  ### Check inputs
  if (!is_bin_events(bins))
    stop('bins os not a bin event table')
  
  ### Initialize
  # Sampling functions
  functions = sampling.functions(names(e), ...)
  # Calculate non-overlapping bin groups
  bin.groups = nonoverlapping.groups(bins)
  # Append original index to bins
  bins = cbind(seq_len(nrow(bins)), bins)
  
  ### For each bin set
  L = lapply(split(bins, bin.groups), function(bins) { 
    
    # Cut events at bin endpoints
    e.cut = cut.events(e, bins)
    # Get index of segments intersecting with each bin
    inbins = find.intersecting(bins, e.cut)
    # Grab first match for each event to convert to bin indices column
    # (each event in only one non-overlapping bin)
    # No match = NA
    bid = apply(inbins, 1, function(x) match(TRUE, x))
    keep = !is.na(bid)
    # Append assignments as last column, removing unassigned events
    e.assigned = cbind(e.cut[keep, ], bid[keep])
    # Apply the functions
    d = do.call(cbind, lapply(functions, function(f) f(e.assigned)))
    # Reinsert empty bins as NA
    kept = sort(unique(bid[keep]))
    if (drop.empty) {
      bins = bins[kept, ]
    } else {
      nb = nrow(bins)
      ind = numeric(nb)
      ind[kept] = seq_along(kept)
      ind[ind == 0] = NA
      d = d[ind, , drop = FALSE]
    }
    # Append bin info
    d = cbind(bins, d)
  })
  
  ### Prepare final output
  # Row bind list of outputs
  # REQUIRES: data.table
  # FIXME: Use conversion to numeric instead?
  # http://stackoverflow.com/questions/5980240/performance-of-rbind-data-frame
  D = as.data.frame(rbindlist(L))
  # Reorder result by original bin index and drop index
  # (has the effect of also making unique names)
  D = D[order(D[[1]]), -1]
  # Apply user names
  if (!is.null(col.names)) {
    names(D)[(length(bins)):length(D)] = col.names
  }
  return(D)
}

# === PLOTTING ===================================================================================================

#' Plot Bins as Barplot
#' 
FIXME: Filling gaps means bar appears when plotting FROM or TO field. Fix by using rect to plot?
#' Send barplot to the barplot grid.
#' @param bins a bin event table.
#' @param i 
#' @param xlim,ylim
#' @param xticks,yticks
#' @param sigfigs
#' @param col
#' @family event plotting
#' @keywords internal
plot_events_single = function(bins, i, 
                       xlim = event_range(bins), ylim = range(c(0, rowSums(bins[i], na.rm = TRUE))),
                       xticks = xlim, yticks = ylim,
                       sigfigs = c(3, 3), col = "grey", ...) {
  
  # REQUIRES sorted, no overlaps
  # Can contain points
  
  ### Fill gaps (with NA)
  # FIXME: as seperate so that width is known but height can be NA even for FROM and TO columns
  bins = sort_events(fill_event_gaps(bins))
  
  ### Barplot
  barplot2(
    t(bins[i]), width = bins$to - bins$from, space = 0, beside = FALSE, 
    xpd = FALSE, axes = FALSE, xlim = as.numeric(xlim - event_range(bins)[[1]]), ylim = as.numeric(ylim), col = col, names.arg = NULL, ...
  )
  
  ### Calculate ticks and labels
  old = par("usr")
  par(usr = c(old[1:2] + event_range(bins)[[1]], old[3:4]))
  
  if (is.null(xticks)) xticks = axTicks(1)
  if (is.null(yticks)) yticks = axTicks(2)
  if (is.null(sigfigs)) {
    xlabels = xticks
    ylabels = yticks
  } else {
    if (length(sigfigs) == 1) sigfigs = rep(sigfigs, 2)
    xlabels = prettyNum(signif(xticks, sigfigs[1]), big.mark = ",")
    ylabels = prettyNum(signif(yticks, sigfigs[2]), big.mark = ",")
  }
  
  ### Draw ticks and labels
  axis(side = 1, at = xticks, labels = xlabels, line = 0.25, tcl = -0.40, col = par("bg"), col.tick = par("fg"), mgp = c(3, 0.65, 0))
  axis(side = 2, at = yticks, labels = ylabels, las = 1)
}

# Barplot grid
# FIXME: Add support for recursive limits, ticks, labels, titles, grid order
# Produces a series of barplots arranged in a grid. The grid is populated from top-left to lower-right.
#' @family event plotting
plot_events = function(events, groups, ...) {
  
  # REQUIRES sorted, no overlaps
  
  # Check input
  nr = length(events)
  nc = length(unique(groups))
  
  # Initialize graphic window and grid layout
  options(warn = -1)
  library(gplots, warn.conflicts = FALSE)
  # mfrow: filled by rows. first field, then group
  # mfcol: filled by cols. first group, then field
  par(mfcol = c(nr, nc), mai = c(0.4,0.4,0.3,0.3), oma = c(0.8,1.6,1.2,0.8), cex.axis = 1.1)
  
  # For each event group
  for (group in unique(groups)) {
    ind = group == groups
    # For each set of variables
    for (i in 1:length(events)) {
      plot.events.grid(events[ind,], i, ...)
    }
  }
  detach(package:gplots)
}