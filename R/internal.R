# Recursive grep returning only exact matches to elements of character vector x
#' @keyword internal
rgrep.exact = function(patterns, x, ..., simplify = TRUE, USE.NAMES = FALSE) {
  sapply(patterns, function(pattern) grep(paste0("^", pattern, "$"), x, ...), simplify = simplify, USE.NAMES = USE.NAMES)
}

# Tests if vector is empty or non-numeric, or contains any non-integer numbers
#' @keyword internal
is.not.integer = function(x) {
  !is.numeric(x) || length(x) == 0 || any(x %% 1 != 0)
}

# Helper for build.function.call
#' @keyword internal
ident <- function(x) {
  y <- as.integer(as.factor(x))
  z <- gsub(" ", "0", format(y, scientific = FALSE))
  return(z)
}

# Create sampling functions for sample_events in nested environments enclosing all their fixed arguments
#' @keyword internal
build.function.call = function(fun, bin.col, data.cols, group.cols = NULL, arg.cols = NULL, arglist = NULL) {
  
  ### Isolate function arguments from executing environment
  e = new.env()
  e$fun = fun
  e$bin.col = bin.col
  e$data.cols = data.cols
  e$group.cols = group.cols
  e$arg.cols = arg.cols
  e$arglist = arglist
  
  ### Create function call where only input is the data
  e$call = function(events) {
    
    # Build table of unique group combinations
    y = events[c(bin.col, group.cols)]
    ynames = names(y)
    if (ncol(y)) {
      grp <- rank(do.call(paste, c(lapply(rev(y), ident), list(sep = "."))), ties.method = "min")
    } else { 
      grp <- integer(NROW(events))
    }
    y <- y[match(sort(unique(grp)), grp, 0L), , drop = FALSE]
    
    # Pre-split additional argument columns
    var.args = lapply(unname(events[arg.cols]), split, f = grp)
    
    # Apply function over each data column, then each subset of data and arguments
    x = events[data.cols]
    z <- lapply(x, function(e) {
      args = c(
        FUN = fun, 
        list((split(e, grp))),
        var.args,
        MoreArgs = list(arglist)
      )
      ans <- do.call(mapply, args)
      ans
    })
    
    # Assemble back into named dataframe
    len <- length(y)
    for (i in seq_along(z)) y[[len + i]] <- z[[i]]
    names(y) <- c(ynames, names(x))
    #row.names(y) <- NULL
    
    # Reshape long result to wide format
    # (if group cols present) 
    if (length(group.cols)) {
      len = length(y)
      # Merge group levels
      groups = seq_along(group.cols) + 1
      y[len + 1] = do.call(paste, c(y[groups], sep = "."))
      # Reshape
      y = reshape(y[-groups], idvar = names(y[1]), timevar = names(y[len + 1]), direction = 'wide')
    }
    
    # Drop bin (id) column
    y = y[-1]
    
    # Sort data columns x.a, x.b, y.a, y.b
    # (rather than x.a, y.a, x.b, y.b)
    return(y[order(names(y))])
  }
  
  # Return function call
  return(e$call)
}

# Parse the function call parameters (from sample.EventTable) into self-enclosed functions
# unnamed arguments are interpreted as column names or indices
# first unnamed argument can be a vector of rows, all others are flattened and passed individually to function
# named arguments are taken literally and passed named and unchanged to the function
# first "by" argument is expected to be vector of column name or indices
# list(sum, 1:2, na.rm = TRUE) => sum(events[1], na.rm = TRUE), sum(events[2], na.rm = TRUE)
# list(sum, 1, 3:4, 5, na.rm = TRUE) => sum(events[1], events[3], events[4], events[5], na.rm = TRUE), ...
# NOTE: Assumes bin column will be appended to the end of the events
#' @keyword internal
sampling.functions = function(col.names, ...) {
  
  ### Initialize inputs
  calls = list(...)
  functions = vector("list", length(calls))
  
  ### Parse each function call
  # FIXME: Replace with an lapply?
  for (i in seq_along(calls)) {
    call = calls[[i]]
    
    ## Pull out group.cols, if present
    group.ind = match("by", names(call))
    if (!is.na(group.ind)) {
      group.cols = rapply(call[group.ind], rgrep.exact, classes = c("character"), how = "replace", x = col.names)
      call[group.ind] = NULL
    } else {
      group.cols = NULL
    }
    
    ## Identify fixed arguments
    args = call[2:length(call)]
    if (is.null(names(args))) 
      names(args) = rep("", length(args))
    fixed = !is.null(names(args)) & names(args) != ""
    
    ## Identify data and argument columns
    # Flatten to single level list (just in case)
    args[!fixed] = lapply(args[!fixed], unlist)
    # Replace field names with corresponding numeric indices
    args[!fixed] = rapply(args[!fixed], rgrep.exact, classes = c("character"), how = "replace", x = col.names)
    # Require integer numeric indices
    if (any(unlist(lapply(args[!fixed], is.not.integer))))
      stop('Required: non-empty integer vectors (after column name matching)')
    
    ## Build function call
    functions[[i]] = build.function.call(
      fun = call[[1]], 
      bin.col = length(col.names) + 1,
      data.cols = args[!fixed][[1]],
      group.cols = unlist(group.cols),
      arg.cols = unlist(args[!fixed][-1]),
      arglist = args[fixed]
      )
  }
  
  ### Return list of all functions
  return(functions)
}