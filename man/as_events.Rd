% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/event_creation.R
\name{as_events}
\alias{as_events}
\alias{as_events.data.frame}
\alias{as_events.matrix}
\alias{as_events.numeric}
\title{Coerce to an Event Table}
\usage{
as_events(x, ...)

\method{as_events}{numeric}(x, ...)

\method{as_events}{matrix}(x, from.col = 1, to.col = 2, ...)

\method{as_events}{data.frame}(x, from.col = 1, to.col = 2, ...)
}
\arguments{
\item{x}{object to be coerced to an event table.}

\item{...}{Additional arguments passed to or used by methods.}

\item{from.col,to.col}{names or indices of the columns in \code{x} containing the event endpoints. Values are swapped as needed to ensure that \code{to > or = from} on all rows.}
}
\description{
Attempts to coerce an object to an event table.
}
\section{Methods (by class)}{
\itemize{
\item \code{numeric}: Expands a numeric vector into two columns of event endpoints.

\item \code{matrix}: Converts the matrix to a data frame, then calls the \code{data.frame} method.

\item \code{data.frame}: Renames \code{from.col} and \code{to.col} to "from" and "to" as needed. Since these column names must be unique, other columns cannot also be called "from" or "to".
}}
\examples{
as_events(1)
as_events(1:5)
as_events(cbind(1:5, 1:5), 1, 2)
as_events(data.frame(x = 1, start = 1:5, stop = 1:5), "start", "stop")
}
\seealso{
\code{\link{events}} for creating event tables and \code{\link{read_events}} for reading files as event tables.
}
