
#===================================================================================================

#' Finds all overlaps between an EventTable and binEventTable. Each bin is uniquely identified by a binID, corresponding to either the index or, if present, the ID field of the corresponding row in the binEventTable. For each overlapping event, likewise identified by an eventID, is listed the overlapping bin-event range (overFROM, overTO) and the ratio of event overlap. For lines, ratio = length of overlap / length of event; for points, typically ratio = 1, but points lying at shared endpoints of adjacent bins are included in both bins with ratio = 0.5. In addition to ratio, the gap to the start of the next event (or bin endpoint) is listed. This is needed for later calculating bin COVERAGE, the total flattened length of events in each bin. Gaps between the bin startpoint and first event are listed on a row with eventID = NA that precedes the event listings. Empty bins are retained as a single row with eventID = NA and a gap equal to the full length of the bin.
#'
#' events - (EventTable) Event data to be binned
#' bins - (binEventTable) Bins to be used for binning
find.binEventOverlaps = function(events, bins) {
  
  # Verify events and bin EventTables
  if (!is.EventTable(events)) stop("Events are not a valid EventTable.\n")
  if (!is.binEventTable(bins)) stop("Bins are not a valid binEventTable.\n")
  
  # If ID not found, hardcode user-supplied row order into ID
  if (!"ID" %in% names(events)) events = cbind(ID = 1:nrow(events), events)
  if (!"ID" %in% names(bins)) bins = cbind(ID = 1:nrow(bins), bins)
  
  # Sort events and bin EventTables
  events = .sort.EventTable(events)
  bins = .sort.EventTable(bins)
  
  # Create an empty bin reference table
  binID = FROM = TO = eventID = overFROM = overTO = ratio = gap = integer(0)
  overlaps = data.frame(binID, FROM, TO, eventID, overFROM, overTO, ratio, gap)
  
  # For each bin...
  for (i in 1:nrow(bins)) {
    
    # Select all the point events within range, including at bin endpoints. 
    # (points lying at shared endpoints of adjacent bins are included twice, these are "halved" later with r = 0.5)
    ispoint = find.points(events)
    incpoints = ispoint & bins$TO[i] >= events$FROM & bins$FROM[i] <= events$TO
    
    # Select all the line events within range, discarding those that touch but do not overlap bin.
    isline = find.lines(events)
    inclines = isline & bins$TO[i] > events$FROM & bins$FROM[i] < events$TO
    
    # Build table of included point and line events
    binEvents = events[incpoints | inclines,]
    
    # If the bin contains events...
    if (nrow(binEvents) > 0) {
      
      # Initialize maxTO at beginning of bin (for gap detection)
      maxTO = bins$FROM[i]
      
      # For each unit...  
      for (j in 1:nrow(binEvents)) {
        
        # Determine the ratio of inclusion based on relative size and position of the unit and bin.
        # if point, ratio = 1 unless lying on shared endpoint of adjacent bins, then ratio = 0.5
        if (binEvents$TO[j] == binEvents$FROM[j]) {
          sharedFROM = i > 1 && binEvents$FROM[j] == bins$FROM[i] && binEvents$FROM[j] == bins$TO[i-1]
          sharedTO = i < nrow(bins) && binEvents$TO[j] == bins$TO[i] && binEvents$TO[j] == bins$FROM[i+1]
          if (sharedFROM | sharedTO) ratio = 0.5
          else ratio = 1
          overlap = c(binEvents$FROM[j], binEvents$TO[j])
          
          # if line, intersect range with bin, and divide result by total length of unit
        } else {
          binRange = bins[i,c("FROM","TO")]
          unitRange = binEvents[j,c("FROM","TO")]
          overlap = intersect.ranges(binRange, unitRange)
          ratio = diff(overlap) / diff(as.numeric(unitRange))
        }
        
        # Initialize gap default
        gap = 0
        
        # Measure any gaps between the current maxTO and the next event or bin endpoint
        # (later, these gaps are added up and inversed to provide COVERAGE, the flattened total length of events in each bin)
        # if unit reaches beyond current maxTO, increase maxTO correspondingly
        if (maxTO < binEvents$TO[j]) maxTO = binEvents$TO[j]
        # gap between beginning of bin and first unit is written to results with eventID = overFROM = overTO = NA, and ratio = 0.
        if (j == 1 && (startGap = binEvents$FROM[j] - bins$FROM[i]) > 0) overlaps = rbind(overlaps, data.frame(binID = bins$ID[i], FROM = bins$FROM[i], TO = bins$TO[i], eventID = NA, overFROM = NA, overTO = NA, ratio = 0, gap = startGap))
        # gap between current maxTO and beginning of next unit
        if (j < nrow(binEvents) && (unitGap = binEvents$FROM[j+1] - maxTO) > 0) gap = unitGap
        # gap between end of bin and last unit
        if (j == nrow(binEvents) && (endGap = bins$TO[i] - maxTO) > 0) gap = endGap
        
        # Append the result to the bin reference table
        overlaps = rbind(overlaps, data.frame(binID = bins$ID[i], FROM = bins$FROM[i], TO = bins$TO[i], eventID = binEvents$ID[j], overFROM = overlap[1], overTO = overlap[2], ratio = ratio, gap = gap))
      }
      
      # If the bin is empty...	
    } else {
      
      # Append the null result to the bin reference table
      overlaps = rbind(overlaps, data.frame(binID = bins$ID[i], FROM = bins$FROM[i], TO = bins$TO[i], eventID = NA, overFROM = NA, overTO = NA, ratio = 0, gap = bins$TO[i] - bins$FROM[i]))
    }
  }
  
  # Return result
  return(overlaps)
}

#' Bins the events in the provided EventTable according to the provided binEventOverlaps. Totals (lengths, species counts, etc) are summed, maximums and minimums are recalculated, and factors (eg. unit type) are split apart (eg. UNIT_TYPE.RI, UNIT_TYPE.RA, ...) and summarized by either total length or # of contributing events for each possible value (including NA). Means (eg. % substrate, average depth) are calculated by weighing in each contributing unit by geographic length. For any given mean field, if a null value is encountered, the contribution of that unit in a bin is ignored by normalizing the contributions of the other units by a reduced bin LENGTH. Bin SIZE is simply TO - FROM. Missing values are ignored in calculating bin counts and maxima.
#'
#' events - (data.frame) EventTable
#' overlaps - (data.frame) binEventOverlaps
#' counts - (character vector) names of fields that represent counts (includes lengths, totals). events are assumed to be uniformly distributed (contribution to bin scaled by ratio of inclusion).
#' means - (character vector) names of fields that represent averages (includes proportions, percentages)
#' maximums - (character vector) names of fields that represent maximums
#' minimums - (character vector) names of fields that represent minimums
#' factors - (character vector) names of fields that represent factors
#' factor.weight - (boolean) if TRUE, factor levels are summarized by weights, scaled by ratio of inclusion (overlapping length being the default) of contributing events; if FALSE, factor levels are summarized by the total number of contributing events.
#' equal.weight - (boolean) if TRUE, means are calculated without weights; if FALSE, means are calculated by weighing each contribution by fraction unit length / bin length (default) or unit weight / bin weight (specified by weightField)
#' na.rm - (boolean) a logical value indicating whether NA values should be stripped before the computation of bin counts, means, maximums, and minimums.
#' weightField - (character) name of field to use for weighting events in the computation of bin means instead of the standard overlapping event length (overTO - overFROM)
apply.binEventOverlaps = function(events, overlaps, counts = NULL, means = NULL, maximums = NULL, minimums = NULL, factors = NULL, factor.weight = TRUE, equal.weight = FALSE, na.rm = TRUE, weightField = NULL) {
  
  # Verify EventTable and binEventOverlaps
  if (!is.EventTable(events)) stop("Events are not a valid EventTable.\n")
  if (!is.binEventTable(overlaps)) stop("Overlaps are not supplied as valid binEventOverlaps, as would be returned by find.binEventOverlaps().\n")
  
  # Verify field names are character vectors
  fieldArguments = c("counts","means","maximums","minimums","factors","weightField")
  characterFields = as.logical(lapply(list(counts,means,maximums,minimums,factors,weightField),FUN=is.character))
  nullFields = as.logical(lapply(list(counts,means,maximums,minimums,factors,weightField),FUN=is.null))
  validFields = characterFields | nullFields
  if (any(!validFields)) stop("One or more field name arguments (",paste(fieldArguments[!validFields & !nullFields], collapse=", "),") are not character vectors.\n")
  
  # Verify field names exist
  userFields = c(counts, means, maximums, minimums, factors, weightField) # NULL values are dropped
  supplied = userFields %in% names(events)
  if (any(!supplied)) stop("One or more requested fields (",paste(userFields[!supplied], collapse=", "),") were not found in the provided EventTable.\n")
  
  # Verify fields names are not assigned to multiple types
  uniqueFields = c(counts, means, maximums, minimums, factors)
  isduplicated = duplicated(uniqueFields)
  if (any(isduplicated)) stop("One or more fields (",paste(uniqueFields[isduplicated], collapse=", "),") were assigned to multiple field type arguments.\n")
  
  # Verify field types (numeric for all but factors)
  numericFields = c(counts, means, maximums, minimums, weightField)
  isnumeric = as.logical(lapply(events[c(numericFields)],FUN=is.numeric))
  if (any(!isnumeric)) stop("All but factor fields must be numeric. However, one or more fields (",paste(numericFields[!isnumeric], collapse=", "),") were found to contain non-numeric values.\n")
  
  # Verify logicals
  fieldArguments = c("factor.weight","equal.weight","na.rm")
  logicalFields = as.logical(lapply(list(factor.weight,equal.weight,na.rm),FUN=is.logical))
  if (any(!logicalFields)) stop("One or more arguments which must be either TRUE or FALSE (",paste(fieldArguments[!logicalFields], collapse=", "),") are not type logical.\n")
  
  # If ID not found, hardcode user-supplied row order into ID. Expected to be same order as when supplied to find.binEventOverlaps().
  if (!"ID" %in% names(events)) events = cbind(ID = 1:nrow(events), events)
  
  # # If weightField not specified, use LENGTH = TO - FROM
  # # (what if LENGTH already exists?)
  # if (is.null(weightField)) {
  # events = transform(events, LENGTH = TO - FROM)
  # weightField = "LENGTH"
  # }
  
  # Order and expand EventTable to match binEventOverlaps eventID, but keep seperate.
  # (avoids duplicate or renamed fields if the two tables are merged)
  data = events[match(overlaps$eventID, events$ID),]
  isevent = !is.na(overlaps$eventID)
  
  # Only fields requested by the user are kept, all others are dropped for processing.
  data = data[c(userFields)]
  
  ## Counts
  # Multiply counts by ratio of inclusion
  # (assumes uniform distribution)
  data[c(counts)] = overlaps$ratio * data[c(counts)]
  
  ## Factors
  # For each factor...
  for (i in factors) {
    
    # convert field to factor
    # (just in case supplied as other type)
    data[,i] = factor(data[[i]])
    
    # find all unique values
    fields = c(levels(data[[i]]))
    
    # if any NA entries exist, add "NA" to fields
    if (sum(is.na(data[i])) > 0) fields = c(fields, NA)
    
    # generate unique field names of form UNIT_TYPE.RI, UNIT_TYPE.RA, ...
    fieldNames = paste(i, ".", fields, sep = "")
    
    # create empty data frame with new fields
    dt = data.frame(t(rep(0,length(fieldNames))))
    names(dt) = fieldNames
    
    # For each unit...
    for (j in 1:nrow(data)) {
      
      # look up value
      value = as.character(binData[i][j,])
      
      # If value is valid, add to corresponding new field...
      if (value %in% fields) {
        value = paste(i, ".", value, sep = "")
        
        # ...either the weight of the current unit
        if (factor.length) {
          binData[value][j,] = binData$LENGTH[j]
          
          # ...or the fraction r of the current unit
        } else {
          binData[value][j,] = binData$r[j]
        }
      }
    }
    
    # Delete factor field
    binData[i] = NULL
  }
  
  ## means
  # Create vector of bin length by summing together the LENGTH of rows with the same BID
  # Multiply means by the ratio of unit length to bin length
  
  for (i in means) {
    binTemp = binData
    if (equal.weight) binTemp$LENGTH = 1
    binTemp$LENGTH[is.na(binTemp[i])] = 0
    binLength = rowsum(binTemp, binTemp$BID, na.rm = TRUE)$LENGTH
    lengthDF = data.frame(BID = unique(binTemp$BID), LENGTH = binLength)
    lengths = merge(binTemp["BID"], lengthDF, by="BID", sort=FALSE)$LENGTH
    binData[i] = (binTemp$LENGTH / lengths) * binTemp[i]
    binData[,i][is.nan(binData[[i]])] = NA
  }
  
  # Add field to count number of all contributing units in a bin
  binData = cbind(binData, r.all = 1)
  binData$r.all[binData$r == 0] = 0
  
  # Sum together rows of equal BID across all fields. 
  binTotals = rowsum(binData, binData$BID, na.rm = TRUE)
  
  # Fix BID field
  binTotals$BID = unique(binData$BID)
  
  # Change false zeroes in counts and means to NA
  for (i in binTotals$BID) {
    sub = subset(binData, BID == i)
    for (j in c(means,counts)) {
      if (sum(is.na(sub[j])) == nrow(sub)) binTotals[binTotals$BID == i,j] = NA
    }
  }
  
  ## maximums
  # For each maximum...
  for (i in maximums) {
    
    # For each bin...
    for (j in binTotals$BID) {
      
      # find bin maximum among the units it contains, ignoring NA values unless only NA are present
      binMax = max(subset(binData, BID == j)[i], na.rm = TRUE)
      if (!is.finite(binMax)) binMax = NA 
      binTotals[j, i] = binMax
    }
  }
  
  ## minimums:
  # For each minimum...
  for (i in minimums) {
    
    # For each bin...
    for (j in binTotals$BID) {
      
      # find bin minimum among the units it contains, ignoring NA values unless only NA are present
      binMin = min(subset(binData, BID == j)[i], na.rm = TRUE)
      if (!is.finite(binMin)) binMin = NA 
      binTotals[j, i] = binMin
    }
  }
  
  key = key[!duplicated(key$BID),]
  
  # Merge the BID, bin FROM, TO and bin totals calculated above. Replace TOTALGAP with its opposite, COVERAGE. Add SIZE (TO - FROM).
  binTotals = transform(binTotals, BID = key$BID, FROM = key$FROM, TO = key$TO, COVERAGE = abs(key$TO - key$FROM) - binTotals$TOTALGAP)
  binTotals$TOTALGAP = NULL
  binTotals = cbind(binTotals, SIZE = binTotals$TO - binTotals$FROM)
  
  # Return result.
  order = unique(c("BID","FROM","TO","SIZE","COVERAGE","LENGTH","r","r.all",counts,means,maximums,minimums,names(binTotals)))
  binTotals = binTotals[,order]
}




#' Find Intersecting Events
#' 
#' \code{find.intersectingEvents} returns a logical index vector of the events in \code{x} that intersect any of the events in \code{y}.
#' 
#' @param x either an \code{EventTable} or a numeric vector. Numeric vectors are interpreted as the positions of point events.
#' @param y either an \code{EventTable} or a numeric vector. Numeric vectors are interpreted as the positions of point events.
#' @param closed logical value indicating whether to treat events as closed intervals. If \code{FALSE} (the default), a shared endpoint is not sufficient for an intersection. Point-point matches are only returned if \code{closed = TRUE}.
#' @param closed.from,closed.to logical values indicating whether to treat events in \code{x} as closed at their \code{FROM} and \code{TO} endpoints, respectively, and events in \code{y} as closed at their \code{TO} and \code{FROM} endpoints, respectively.
#' @param verbose logical value indicating whether to print verbose messages.
#' @return A logical vector of the same length as the number of rows in \code{x}.
#' @seealso \code{\link{EventTable}}.
#' @export
# find.intersecting = function(x, y, closed = FALSE, closed.from = closed, closed.to = closed, verbose = FALSE) {
#   
#   # Verify x, coercing to EventTable if necessary
#   if (is.numeric(x)) y = EventTable(x,x)
#   if (!is.EventTable(x)) stop("X is not an EventTable.\n")
#   
#   # Verify y, coercing to EventTable if necessary
#   if (is.numeric(y)) y = EventTable(y,y)
#   if (!is.EventTable(y)) stop("Y is not an EventTable.\n")
#   
#   # Verify logicals
#   fieldArguments = c("closed","closed.from","closed.to")
#   logicalFields = as.logical(lapply(list(closed,closed.from,closed.to),FUN=is.logical))
#   if (any(!logicalFields)) stop("One or more arguments which must be either TRUE or FALSE (",paste(fieldArguments[!logicalFields], collapse=", "),") are not type logical.\n")
#   
#   # Flatten y to coverage
#   if (nrow(y) > 1) ranges = coverage.EventTable(y)
#   
#   # For each event in x, find whether it intersects one or more ranges of y
#   N = nrow(x)
#   xpos = rep(FALSE, N)
#   for (n in 1:N) {
#     
#     # (simple & inefficient, consider sorting first)
#     if (closed.from) {
#       leftmatch = ranges$TO >= x$FROM[n]
#     } else {
#       leftmatch = ranges$TO > x$FROM[n]
#     }
#     if (closed.to) {
#       rightmatch = ranges$FROM <= x$TO[n]
#     } else {
#       rightmatch = ranges$FROM < x$TO[n]
#     }
#     
#     # Return logical vector of intersecting events
#     xpos[n] = xpos[n] || any(leftmatch & rightmatch)
#   }
#   
#   # Return logical of events in x intersecting y
#   return(xpos)
# }