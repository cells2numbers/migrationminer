% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/track.R
\name{valid_observation_time}
\alias{valid_observation_time}
\title{calculate valid observation time as sum of the length of all
valid tracks divided by the sum of the length of all tracks}
\usage{
valid_observation_time(tracks, min_path_length = 19)
}
\arguments{
\item{tracks}{data frame with track objects}

\item{min_path_length}{minimum length of a valid track}
}
\value{
valid_observation_time
}
\description{
calculate valid observation time as sum of the length of all
valid tracks divided by the sum of the length of all tracks
}
\examples{
 data <- tibble::tibble(
   Metadata_timePoint = c(1:5),
   Location_Center_X = c(1, 2, 3, 4, 5),
   Location_Center_Y = c(1, 1, 1, 1, 1),
   TrackObjects_Label = c(rep(1, 5))
 )
 data <- dplyr::group_by_(data, "TrackObjects_Label")
 tracks <- migrationminer::track(data, "TrackObjects_Label")
 min_path_length <- 5
 vot <-   migrationminer::valid_observation_time(tracks, min_path_length)
}
