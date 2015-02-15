#' Plot Events as Bar Plots
#' 
#' Plots an event table as a grid of bar plots.
#' 
#' This function is intended for plotting the output of \code{\link{sample_events}}. Given a groupping variable for the rows of the event table (e.g., groups of bins of different sizes), and groups of columns to plot, bar plots are drawn for each combination onto a grid event groups high by column groups wide. In each plot, the specified event table columns are plotted together as stacked bars. Negative and positive values are stacked seperately from the \code{y = 0} baseline. Events with \code{NA} are not shown, differentiating them from zero-valued events which are drawn as thin black lines. Point events are drawn as thin vertical lines. Overlapping events are drawn as overlapping bars, so it is best to use \code{\link{sample_events}} first to flatten the data before plotting.
#'
#' @param e an event table.
#' @param group.col name or index of event table column defining the event groupping for plotting. If \code{NULL}, the events are treated as one group. Group \code{NA} is not plotted.
#' @param groups vector of values from \code{group.col} specifying which groups to plot. If \code{NULL}, all groups are plotted by order of first appearance in \code{group.col}.
#' @param data.cols names or indices of columns to plot, given as a list of character or numeric vectors. If multiple columns are specified, their bars are stacked together in one plot. Names are interpreted as regular expressions (\code{\link{regex}}) matching full column names. If \code{NULL}, all columns not named \code{from}, \code{to}, or \code{group.col} are each plotted seperately.
#' @param dim the row and column dimensions of the grid. If \code{NULL}, the grid is column groups (rows) by event groups (columns).
#' @param byrow if \code{TRUE}, plots will be added by rows, rather than columns, to the grid.
#' @param main titles for each plot (recycled as necessary). If \code{NULL}, plots are titled by the column names. Set \code{main = ""} to hide.
#' @param xlim,ylim limits for the x and y axes for all plots. If \code{NULL}, limits are set to the range of the data and the y limits extended as needed to include 0.
#' @param xticks,yticks the values to label on the x and y axes of all plots. If \code{NULL}, only the min and max x and y are labeled (and 0 as needed for y). If \code{\link{axTicks}}, the function will be used to generate R default tick marks.
#' @param xlabels,ylabels the character strings to display at the xticks and yticks, coerced to character vectors and recyled as necessary. If \code{NULL}, the coordinate positions of the ticks (user or auto) are used as labels, modified by \code{sigfigs}.
#' @param plot.grid if \code{TRUE}, a lined horizontal grid is plotted at the yticks.
#' @param sigfigs the maximum significant figures of x and y axis labels.
#' @param col color(s) for the bars in each plot. If \code{NULL}, bars are transparent. By default, a grey palette is used.
#' @param border color(s) for bar borders in each plot. Use border = NA to omit borders.
#' @param lty line type(s) for bar borders in each plot.
#' @param lwd line width(s) for bar borders in each plot.
#' @param xpd logical value or \code{NA}. If \code{FALSE}, all plotting is clipped to the plot region, if \code{TRUE}, all plotting is clipped to the figure region, and if \code{NA}, all plotting is clipped to the device region.
#' @param mar numerical vector of the form c(bottom, left, top, right) giving the size of the inner margins for each plot in lines of text.
#' @param oma numeric vector of the form c(bottom, left, top, right) giving the size of the outer margins in lines of text.
#' @param ... additional arguments passed to \code{\link{plot}}.
#' @seealso \code{\link{seq_bin_events}} for generating groups of sequential bins, \code{\link{sample_events}} to populate bins with event data.
#' @export
#' @examples
#' e = events(from = c(0, 10, 15, 25), to = c(10, 20, 25, 40), length = c(10, 10, 10, 15), x = c(1, 2, 1, 1), f = c('a', 'b', 'a', 'a'))
#' b = seq_bin_events(event_coverage(e), c(8, 4, 2, 1))
#' bins = sample_events(e, b, list(sum, 'length'), scaled.cols = 'length')
#' plot_events(bins, group.col = 'group', data.cols = 'length', dim = c(2,2))
plot_events = function(e, group.col = NULL, groups = NULL, data.cols = NULL, dim = NULL, byrow = FALSE, main = NULL, xlim = NULL, ylim = NULL, xticks = NULL, yticks = NULL, xlabels = NULL, ylabels = NULL, plot.grid = FALSE, sigfigs = c(3, 3), col = NULL, border = par("fg"), lty = par("lty"), lwd = par("lwd"), xpd = FALSE, mar = c(2.1, 2.75, 1.5, 0.5),  oma = rep(2, 4), ...) {
  
  ### Initialize
  if (is.null(group.col)) {
    # all events in one group
    group.seq = numeric(nrow(e))
    groups = 0
  } else {
    group.seq = e[[group.col]]
    if (is.null(groups)) {
      # all event groups
      groups = unique(group.seq)
      groups = groups[!is.na(groups)]
    } else {
      # user-specified event groups
      groups = unique(groups)
    }
  }
  if (is.null(data.cols)) {
    group.col = if_else(is.character(group.col), group.col, names(e)[group.col])
    data.cols = c(list(), setdiff(names(e), c("from", "to", group.col)))
  } else {
    data.cols = rapply(as.list(data.cols), rgrep_exact, x = names(e), classes = "character", how = "replace")
    ind2name = function(x) names(e)[ind]
    data.cols = rapply(as.list(data.cols), ind2name, classes = "numeric", how = "replace")
    data.cols = lapply(data.cols, unlist)
  }
  if (is.null(col)) {
    col = grey.colors(length(data.cols))
  }
  
  # grid dimensions
  nr = length(data.cols)
  nc = length(groups)
  n = nr * nc
  if (!byrow) {
    par(mfcol = if_else(is.null(dim), c(nr, nc), dim))
  } else {
    par(mfrow = if_else(is.null(dim), c(nr, nc), dim))
  }
  par(mar = mar, oma = oma)
  #xlab = if_else(is.null(xlab), xlab, rep_len(xlab, n))
  #ylab = if_else(is.null(ylab), ylab, rep_len(ylab, n))
  main = rep_len(if_else(is.null(main), sapply(data.cols, paste, collapse = "+"), main), n)
  
  ### Plot each combination
  j = 1
  for (group in groups) {
    ind = group == group.seq
    # For each set of variables
    for (i in seq_along(data.cols)) {
      plot_events_single(e[ind, , drop = FALSE], data.cols[[i]], main = main[j], xlim = xlim, ylim = ylim, xticks = xticks, yticks = yticks, sigfigs = sigfigs, col = col, border = border, lty = lty, lwd = lwd, xpd = xpd, plot.grid = plot.grid, xlabels = xlabels, ylabels = ylabels, xlab = "", ylab = "", ...)
      j = j + 1
    }
  }
}

#' Plot Events as Bars
#' 
#' Plots event table columns as vertical bars.
#' 
#' The specified event table columns are plotted together as stacked bars. Negative and positive values are stacked seperately from the \code{y = 0} baseline. Events with \code{NA} are not shown, differentiating them from zero-valued events which are drawn as thin black lines. Point events are drawn as thin vertical lines. Overlapping events are drawn as overlapping bars, so it is best to use \code{\link{sample_events}} first to flatten the data before plotting.
#' 
#' @param e an event table.
#' @param cols names or indices of the event table columns to plot together as stacked bars.
#' @param xlim,ylim limits for the x and y axes. If \code{NULL}, limits are set to the range of the data and the y limits extended as needed to include 0.
#' @param xticks,yticks the values to label on the x and y axes. If \code{NULL}, only the min and max x and y are labeled (and 0 as needed for y). If \code{\link{axTicks}}, the function will be used to generate R default tick marks.
#' @param plot.grid if \code{TRUE}, a lined horizontal grid is plotted at the yticks.
#' @param sigfigs the maximum significant figures of x and y axis labels.
#' @param main an overall title for the plot.
#' @param xlab,ylab titles for the x and y axes.
#' @param col color(s) for the bars. If \code{NULL}, bars are transparent. By default, a grey palette is used.
#' @param border color(s) for bar borders. Use border = NA to omit borders.
#' @param lty line type(s) for bar borders.
#' @param lwd line width(s) for bar borders.
#' @param xpd logical value or \code{NA}. If \code{FALSE}, all plotting is clipped to the plot region, if \code{TRUE}, all plotting is clipped to the figure region, and if \code{NA}, all plotting is clipped to the device region.
#' @param ... additional arguments passed to \code{\link{plot}}.
#' @seealso \code{\link{plot_events}}.
#' @keywords internal
#' @examples
#' e = events(from = c(0, 10, 20, 30), to = c(10, 20, 30, 40), x = c(2, -1, -1, NA) * 100, y = c(1, 1, -1, 0))
#' plot_events_single(e, c("x", "y"), plot.grid = TRUE, yticks = axTicks, xticks = axTicks, main = "Hello")
plot_events_single = function(e, cols, xlim = NULL, ylim = NULL, xticks = NULL, yticks = NULL, xlabels = NULL, ylabels = NULL, main = NULL, xlab = NULL, ylab = NULL, plot.grid = FALSE, sigfigs = c(3, 3), col = grey.colors(length(cols)), border = par("fg"), lty = par("lty"), lwd = par("lwd"), xpd = FALSE, ...) {
  
  # Compute plot limits
  if (is.null(xlim)) {
    xlim = c(min(e$from), max(e$to))
  }
  if (is.null(ylim)) {
    sum.neg = function(x, ...) sum(x[x < 0], ...)
    sum.pos = function(x, ...) sum(x[x > 0], ...)
    y.neg = apply(e[cols], 1, sum.neg, na.rm = TRUE)
    y.pos = apply(e[cols], 1, sum.pos, na.rm = TRUE)
    ylim = range(c(0, y.neg, y.pos))
  }
  
  # Initialize plot
  plot(xlim, ylim, type = 'n', axes = FALSE, xlim = xlim, ylim = ylim, main = main, ylab = ylab, xlab = xlab, ...)
  
  # Axis ticks
  if (is.null(xticks)) {
    xticks = xlim
  } else if (identical(xticks, axTicks)) {
    xticks = axTicks(1)
  }
  if (is.null(yticks)) {
    yticks = unique(c(0, ylim))
  } else if (identical(yticks, axTicks)) {
    yticks = unique(c(0, axTicks(2)))
  }
  if (length(sigfigs) == 1) {
    sigfigs = rep(sigfigs, 2)
  }
  if (is.null(xlabels)) {
    if (is.null(sigfigs)) {
      xlabels = xticks
    } else {
      xlabels = prettyNum(signif(xticks, sigfigs[1]), big.mark = ",")
    }
  } else {
    xlabels = rep_len(xlabels, length(xticks))
  }
  if (is.null(ylabels)) {
    if (is.null(sigfigs)) {
      ylabels = yticks
    } else {
      ylabels = prettyNum(signif(yticks, sigfigs[2]), big.mark = ",")
    }
  } else {
    ylabels = rep_len(ylabels, length(yticks))
  }
  
  # Draw axes, ticks, and labels
  if (!identical(NA, xticks)) {
    axis(side = 1, at = xticks, labels = xlabels, col = NA, col.ticks = par("fg"))
  }
  if (!identical(NA, yticks)) {
    if (plot.grid) {
      axis(side = 2, at = yticks, labels = FALSE, tck = 1, lty = 3)
    }
    axis(side = 2, at = yticks, labels = ylabels, las = 1)
  }
  
  # Plot stacked barplots
  # (for each set of bars, stack positive on positive and negative on negative)
  nc = length(cols)
  ne = nrow(e)
  col = if_else(is.null(col), col, rep_len(col, nc))
  border = if_else(is.null(border), border, rep_len(border, nc))
  lwd = if_else(is.null(lwd), lwd, rep_len(lwd, nc))
  lty = if_else(is.null(lty), lty, rep_len(lty, nc))
  base.pos = rep(0, ne)
  base.neg = rep(0, ne)
  for (i in seq_len(nc)) {
    h = e[[cols[i]]]
    y0 = ifelse(h > 0, base.pos, base.neg)
    y1 = y0 + h
    rect(e$from, y0, e$to, y1, col = col[i], border = border[i], lwd = lwd[i], lty = lty[i], xpd = xpd)
    base.pos = ifelse(!is.na(h) & h > 0, base.pos + h, 0)
    base.neg = ifelse(!is.na(h) & h < 0, base.neg + h, 0)
  }
}

plot_event_endpoints = function(e, ...) {
  e = cbind(e, 1)
  plot_events_single(e, length(e), col = "grey", border = par("bg"), lwd = 3, xticks = axTicks, xlab = '', ylab = '', xlim = c(0, 100), ...)
}