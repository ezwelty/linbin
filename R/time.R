#' Convert event endpoints to date-times
#' 
#' @param e Event table or atomic vector.
#' @param tz Time zone (see \code{\link[base]{timezones}}).
#' @param origin Date-time object (see \code{\link[base]{as.POSIXct}}).
#' @export
#' @examples
#' t <- as.POSIXct("1970-01-01", tz = "UTC") + 0:4
#' e <- events(t)
#' to_datetime(e)
#' to_datetime(e$from)
to_datetime <- function(e, tz = "UTC", origin = as.POSIXct("1970-01-01", tz = "UTC")) {
  if (is.atomic(e)) {
    as.POSIXct(e, tz = tz, origin = origin)
  } else {
    e$from <- as.POSIXct(e$from, tz = tz, origin = origin)
    e$to <- as.POSIXct(e$to, tz = tz, origin = origin)
    e
  }
}

#' Convert event endpoints to dates
#' 
#' @param e Event table or atomic vector.
#' @param origin Date object (see \code{\link[base]{as.Date}}).
#' @export
#' @examples
#' t <- as.Date("1970-01-01") + 0:4
#' e <- events(t)
#' to_date(e)
#' to_date(e$from)
to_date <- function(e, origin = as.Date("1970-01-01")) {
  if (is.atomic(e)) {
    as.Date(e, origin = origin)
  } else {
    e$from <- as.Date(e$from, origin = origin)
    e$to <- as.Date(e$to, origin = origin)
    e 
  }
}