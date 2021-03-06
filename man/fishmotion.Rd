% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets.R
\name{fishmotion}
\alias{fishmotion}
\title{Fish Movements}
\format{
Two data frames \code{motion} and \code{origin} with 1,140 rows and 149 rows of 4 variables, respectively.
}
\source{
Armstrong, J. B., D. E. Schindler, C. P. Ruff, G. T. Brooks, K. E. Bentley, and C. E. Torgersen. 2013. Diel horizontal migration in streams: juvenile fish exploit spatial heterogeneity in thermal and trophic resources. Ecology 94:2066-2075. DOI: \doi{10.1890/12-1200.1}
}
\description{
A pair of event tables (in a list) documenting the movements of tagged Coho Salmon (\emph{Oncorhynchus kisutch}) in Bear Creek (Southwest Alaska, USA) for 29 July - 19 August 2008. Table \code{motion} lists individual fish residence time intervals in each of three stream regions, while table \code{origin} lists the study-wide residence time of each fish and the stream region in which the fish was first tagged.
}
\details{
\itemize{
  \item from, to - start and end times as seconds since 1970-01-01 UTC (POSIXct)
  \item fish.id - unique identifier for each fish
  \item region - stream region (1 = 0 - 930 m, a cold downstream region with abundant and spawning sockeye salmon; 2 = 930 - 1360 m, a cold middle region with few if any sockeye salmon; 3 = > 1360 m, a warm upstream region where sockeye salmon were absent)
}
}
