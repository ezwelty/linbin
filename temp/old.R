find.intersecting = function(bins, events, closed.from = FALSE, closed.to = FALSE) {
  
  # EventTables?
  #validate.EventTable(events)
  #validate.EventTable(bins)
  
  ### Initialize bins
  # convert to data matrix
  #mb = data.matrix(bins[c("FROM", "TO")])
  # sort by decreasing TO
  #mb = mb[order(mb[, 1], mb[, 2])]
  
  ### Initialize events
  #me = data.matrix(events[c("FROM", "TO")])
  # sort
  #me = me[order(me[, 1], me[, 2])]
  
  #limits = findInterval(mb, events[, 1])
  #ne = nrow(events)
  
  # Sort events
  # (store original ID in first column)
  #events = cbind(seq_len(ne), events)
  #if (is.unsorted.EventTable(events))
  #  events = .sort.EventTable(events)
  # Sort bins
  #bins = cbind(seq_len(ne), events)
  #if (is.unsorted.EventTable(events))
  #  events = .sort.EventTable(events)
  
  # Pre-allocate start of search at 1st event
  #start.from = 1
  
  ### 
  inbin = apply(bins[c("FROM", "TO")], 1, function(bin) {
    
    # Pre-filter
    # REQUIRES sorted events
    #limits = findInterval(bin, mevents[seq.int(start.from, ), 1])
    #     if (limits[2] == 0) {
    #       # No events!
    #       return(inbin)
    #     }
    #     if (limits[1] == 0) {
    #       limits[1] == 1
    #     }
    #     subset = seq.int(limits[1], limits[2])
    
    # Find events intersecting the bin
    if (closed.from) {
      inleft = bin[1] <= events$TO
    } else {
      inleft = bin[1] < events$TO
    }
    if (closed.to) {
      inright = bin[2] >= events$FROM
    } else {
      inright = bin[2] > events$FROM
    }
    inbin = inleft & inright
    
  })
  
  #   # Return indices (for original, unsorted events)
  #   if (indices) {
  #     # row indices, vector, matrix, or list
  #     #return(apply(inbin, 2, which))
  #     # always list:
  #     return(lapply(split(inbin, col(inbin)), which))
  #   } else {
  #     # logical, matrix
  #     return(inbin)
  #   }
  return(inbin)
}

#' Event Coverage Level (INCOMPLETE)
#' 
#' \code{ncoverage.EventTable} returns the flattened coverage, or union, of all the events in an \code{EventTable}. In other words, the result consists of the intervals, of maximal length, over which the number of events is always 1 or greater. Conversely, \code{gaps.EventTable} returns the gaps in event coverage; the intervals, of maximal length, over which there are 0 events.
#'
#' @param events object of class \code{EventTable}.
#' @param verbose logical value indicating whether to print verbose messages.
#' @return \code{coverage.EventTable} and \code{gaps.EventTable} both return an \code{EventTable} containing only event endpoints (\code{FROM} and \code{TO}). If no gaps are present, \code{gaps.EventTable} returns \code{NA}.
#' @seealso \code{\link{EventTable}}, \code{\link{range.EventTable}} for event range (coverage with gaps ignored).
#' @export ncoverage.EventTable
#' @rdname ncoverage.EventTable
ncoverage.EventTable = function(events, plot = FALSE, verbose = FALSE) {
  
  # Verify EventTable
  if (!is.EventTable(events)) stop("Argument is not an EventTable.\n")
  
  # Cut events at event endpoints
  temp = cut.EventTable(events, events)
  temp = as.EventTable(aggregate(list(events = rep(1, nrow(temp))), temp[c("FROM","TO")], length))
  ptrows = which(find.points(temp))
  for (i in ptrows) {
    temp$events[i] = sum(find.intersecting(temp$FROM[i], temp, closed.lines = T))
  }
  
  # plot (crude)
  if (plot) {
    barplot(temp$events, width = temp$TO - temp$FROM, col = 'black', )  
  }
  
  # Return result
  return(as.EventTable(temp))
}

#' Resample Events at Bins
#' 
#' \code{sample.EventTable} samples events at the specified bins.
#' 
#' Totals (lengths, species counts, etc) are summed, maximums and minimums are recalculated, and factors (eg. unit type) are split apart (eg. UNIT_TYPE.RI, UNIT_TYPE.RA, ...) and summarized by either total length or # of contributing events for each possible value (including NA). Means (eg. % substrate, average depth) are calculated by weighing in each contributing unit by geographic length. All calculations inherit NA handling from na.rm.
#'
#' @param events object of class \code{EventTable}.
#' @param bins object of class \code{EventTable} containing only non-overlapping line events.
#' @param totals names of columns that represent totals.
#' @param means names of columns that represent means.
#' @param maximums names of columns that represent maximums.
#' @param minimums names of columns that represent minimums.
#' @param factors names of columns that represent factors.
#' @param factor.type 'weights','length','coverage'
#' @param mean.type 'weights',length','coverage'
#' @param mean.weights name of column to use for weighing event contributions when computing bin means.
#' @param na.rm a logical value indicating whether NA values should be stripped before the computation of bin totals, means, maximums, and minimums.
#' @param verbose logical value indicating whether to print verbose messages.
#' @return An \code{EventTable} with the same number of rows as \code{bins} and the resampled events data columns. Columns other than FROM and TO in bins are also included if these do not conflict with the data column names.
#' @seealso \code{\link{crop.EventTable}}.
#' @export
resample.EventTable = function(events, bins, totals = NULL, means = NULL, maximums = NULL, minimums = NULL, factors = NULL, weightCol = rep(1,nrow(events)), sliding.mean = FALSE, na.rm = TRUE, verbose = FALSE) {
  
  # Verify events and bins
  if (!is.EventTable(events)) stop("events is not a valid EventTable.\n")
  if (!validate.binEventTable(bins)) stop("bins is not a valid bin EventTable.\n")
  
  # Verify data columns
  validCols = .validate.dataCols(events, totals = totals, means = means, maximums = maximums, minimums = minimums, factors = factors, weightCol = weightCol, verbose = verbose)
  if (!validCols && verbose) stop("Column names are not valid.")
  if (!validCols) stop("Column names are not valid. Use verbose = TRUE for more details.")
  
  ## Crop events to bins
  events = crop.EventTable(events, bins, totals = totals)
  
  ## Add new columns to bins
  # data columns
  dataCols = c(totals, means, maximums, minimums)
  
  # factor columns
  # convert factor columns to factors first, just in case supplied as other type
  events[factors] = lapply(events[factors], factor, exclude = NULL)
  flevels = lapply(events[factors], levels)
  factorCols = paste(rep(factors,lapply(flevels,length)), unlist(as.matrix(flevels)), sep=".")
  
  # add columns (as NA)
  # (remove any existing columns that conflict)
  newCols = c(dataCols, factorCols)
  oldCols = !(names(bins) %in% newCols)
  bins = bins[names(oldCols)]
  bins[newCols] = NA
  
  ## Find points lying on shared boundaries of bins
  # (divide their totals by half)
  dbledges = duplicated(c(bins$FROM, bins$TO))
  pointpos = find.points(events)
  sharedpts = pointpos & events$FROM %in% dbledges
  events[sharedpts,totals] = 0.5 * events[sharepts,totals]
  
  ## For each bin, recompute data fields
  N = nrow(bins)
  for (i in 1:N) {
    
    # Find intersecting events
    # write internal find.intersecting for cut events...
    # or just code inline
    # (including edge points, hence dividing their totals by half first)
    xpos = find.intersecting(bins[i,], events, closed.lines = FALSE, closed.points = TRUE)
    
    # Totals
    bins[i,totals] = colSums(events[xpos,totals], na.rm = na.rm)
    
    # Minimums
    bins[i,minimums] = lapply(events[xpos,minimums], min, na.rm = na.rm)
    
    # Maximums
    bins[i,maximums] = lapply(events[xpos,maximums], max, na.rm = na.rm)
    
    # Means
    # type: length
    # point events ignored
    xlpos = xpos & !pointpos
    weights = events$TO[xlpos] - events$FROM[xlpos]
    weightsums = colSums(as.matrix(!is.na(events[xlpos,means])) * weights)
    bins[i,means] = colSums(weights * events[xlpos,means], na.rm = na.rm) / weightsums
    
    # Factors
    
    # count
    f = 1
    by(events, events[factors[f]], nrow)
    paste(factors[f],names(by(events,events$FACTOR,nrow)),sep=".")
    
    # length
    
    by(events$TO-events$FROM, events[factors[f]], sum)
    
    # coverage
    
    lapply(by(events, events[factors[f]], coverage.EventTable), length.EventTable)
    
    
  }
  
  
  # For sliding mean, dissolve events 
  if (sliding.mean && length(means) > 0) {
    events = cut.EventTable(events, events, totals = totals)
  }
  
}
