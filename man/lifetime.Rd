% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/track.R
\name{lifetime}
\alias{lifetime}
\title{Lifetime of a track object.}
\usage{
lifetime(tracks, t_var = "Metadata_timePoint")
}
\arguments{
\item{tracks}{data frame with track objects}

\item{t_var}{variable name / columne name used for time coordinates}
}
\value{
Calculate life time of each track object
}
\description{
Lifetime of a track object.
}
\examples{
 data <- tibble::tibble(
   Metadata_timePoint = c(1:5),
   Location_Center_X = c(1, 2, 3, 4, 5),
   Location_Center_Y = c(1, 1, 1, 1, 1),
   TrackObjects_Label = c(rep(1, 5))
 )
 tracks <- migrationminer::displace(data,"TrackObjects_Label")
 lifetime <-  migrationminer::lifetime(tracks)

}
