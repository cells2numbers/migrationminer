#' Compute track statistics
#'
#' @param population, single cell data
#' @param strata, column name storing the track label
#' @return track
#' @importFrom magrittr %>%
#' @examples
#'  data <- tibble::data_frame(
#'    Metadata_timePoint = c(1:5),
#'    Location_Center_X = c(1, 2, 3, 4, 5),
#'    Location_Center_Y = c(1, 1, 1, 1, 1),
#'    TrackObjects_Label = c(rep(1, 5))
#'  )
#'  data <- dplyr::group_by_(data,'TrackObjects_Label')
#'  tracks <- neutrominer::track(data,'TrackObjects_Label')
#' @export
track <- function(population, strata) {
  # process `population`, which is the data you get from CellProfiler
  tracks <- displace(population, strata)

  features <- list(
    angle(tracks),
    chemotaxis_index(tracks),
    directionality(tracks),
    distance(tracks),
    directional_persistence(tracks),
    forward_migration_index(tracks),
    lifetime(tracks),
    mean_squared_displacement(tracks, tau = 2),
    sector_analysis(tracks),
    speed(tracks))

  return(Reduce(function(...) merge(..., all = TRUE, by = strata), features))
}

#' Add spatial displacement per frame for each track object
#'
#' @param population, data frame storing single cell data
#' @param strata, column name storing the track label
#' @return displacement
#' @examples
#'  data <- tibble::data_frame(
#'    Metadata_timePoint = c(1:5),
#'    Location_Center_X = c(1, 2, 3, 4, 5),
#'    Location_Center_Y = c(1, 1, 1, 1, 1),
#'    TrackObjects_Label = c(rep(1, 5))
#'  )
#'  tracks <- neutrominer::displace(data,'TrackObjects_Label')
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @export
#'


displace <- function(population, strata) {
  dplyr::right_join(
    population %>%
      dplyr::mutate(Helper_order =  order(Metadata_timePoint)) %>%
      dplyr::select_(.dots = c(strata, 'Location_Center_X', 'Location_Center_Y','Helper_order')) %>%
      #dplyr::mutate(Metadata_timePoint2 =  c(0,Metadata_timePoint[1:(length(Metadata_timePoint) - 1)]) ) %>%
      #dplyr::mutate(Metadata_timePoint2 =  (Metadata_timePoint - 1) ) %>%
      dplyr::mutate(Helper_order =  Helper_order - 1),
      #dplyr::filter(Metadata_timePoint2 != -1) %>%
      #dplyr::select(-Metadata_timePoint),
    population %>% dplyr::mutate(Helper_order =  order(Metadata_timePoint)),
    by = (.dots = c(strata, "Helper_order"))
    ) %>%
  dplyr::mutate(Track_dX = Location_Center_X.x - Location_Center_X.y) %>%
  dplyr::mutate(Track_dY = Location_Center_Y.x - Location_Center_Y.y) %>%
  dplyr::select(-Location_Center_X.x, -Location_Center_Y.x) %>%
  dplyr::rename(Location_Center_X = Location_Center_X.y) %>%
  dplyr::rename(Location_Center_Y = Location_Center_Y.y) %>%
  #dplyr::rename(Metadata_timePoint = Metadata_timePoint2) %>%
  dplyr::select(-Helper_order) %>%
  dplyr::mutate(TrackObjects_Distance_Traveled = sqrt(Track_dX^2 + Track_dY^2))
}

#' Add spatial displacement per frame for each track object
#'
#' @param tracks data frame with single cell data
#' @return displacement
#' @examples
#'  data <- tibble::data_frame(
#'    Metadata_timePoint = c(1:5),
#'    Location_Center_X = c(1, 2, 3, 4, 5),
#'    Location_Center_Y = c(1, 1, 1, 1, 1),
#'    TrackObjects_Label = c(rep(1, 5))
#'  )
#'  tracks <- neutrominer::displace(data,'TrackObjects_Label')
#'  speed <- neutrominer::speed(tracks)
#' @importFrom magrittr %>%
#' @export
speed <- function(tracks) {
  tracks %>%
    dplyr::summarize(Track_Length = n(),
      Track_Speed = sum(TrackObjects_Distance_Traveled, na.rm = TRUE) / (n() - 1),
      Track_Speed_max = max(TrackObjects_Distance_Traveled, na.rm = TRUE),
      Track_Speed_X = sum(Track_dX, na.rm = TRUE) / (n() - 1),
      Track_Speed_Y = sum(Track_dY, na.rm = TRUE) / (n() - 1)) %>%
    dplyr::select(-Track_Length)

}

#' Compute the forward migration index of a tracked object
#'
#' @param tracks data frame with track objects
#' @return forward migration index
#' @examples
#'  data <- tibble::data_frame(
#'    Metadata_timePoint = c(1:5),
#'    Location_Center_X = c(1, 2, 3, 4, 5),
#'    Location_Center_Y = c(1, 1, 1, 1, 1),
#'    TrackObjects_Label = c(rep(1, 5))
#'  )
#'  tracks <- neutrominer::displace(data,'TrackObjects_Label')
#'  forward_migration_index <- neutrominer::forward_migration_index(tracks)
#'
#' @importFrom magrittr %>%
#' @importFrom utils tail
#' @export
forward_migration_index <- function(tracks) {
  s <- tracks %>%
    dplyr::summarize(Track_Integrated_Distance_Traveled = sum(TrackObjects_Distance_Traveled, na.rm = TRUE),
      Track_Displacement_X = tail(Location_Center_X, n = 1) - Location_Center_X[1],
      Track_Displacement_Y = tail(Location_Center_Y, n = 1) - Location_Center_Y[1]
    ) %>%
    dplyr::mutate(Track_xFMI = Track_Displacement_X / Track_Integrated_Distance_Traveled,
      Track_yFMI = Track_Displacement_Y / Track_Integrated_Distance_Traveled) %>%
    dplyr::select(-Track_Integrated_Distance_Traveled, -Track_Displacement_X, -Track_Displacement_Y )
}

#' Calculate lifetime of a track object.
#'
#' @param tracks data frame with track objects
#' @return Calculate life time of each track object
#' @examples
#'  data <- tibble::data_frame(
#'    Metadata_timePoint = c(1:5),
#'    Location_Center_X = c(1, 2, 3, 4, 5),
#'    Location_Center_Y = c(1, 1, 1, 1, 1),
#'    TrackObjects_Label = c(rep(1, 5))
#'  )
#'  tracks <- neutrominer::displace(data,'TrackObjects_Label')
#'  lifetime <-  neutrominer::lifetime(tracks)
#'
#' @importFrom magrittr %>%
#' @export
lifetime  <- function(tracks) {
  tracks %>%
    dplyr::summarize(
      Track_Length = n(),
      Track_Life_Time = length(unique(Metadata_timePoint)),
      Track_One_Cell = length(unique(Metadata_timePoint)) == length(Metadata_timePoint) )
}

#' Calculate angle of a track object.
#'
#' @param tracks data frame with track objects
#' @return The angle of each track
#' @examples
#'  data <- tibble::data_frame(
#'    Metadata_timePoint = c(1:5),
#'    Location_Center_X = c(1, 2, 3, 4, 5),
#'    Location_Center_Y = c(1, 1, 1, 1, 1),
#'    TrackObjects_Label = c(rep(1, 5))
#'  )
#'  tracks <- neutrominer::displace(data,'TrackObjects_Label')
#'  angle <-  neutrominer::angle(tracks)
#'
#' @importFrom magrittr %>%
#' @importFrom utils tail
#' @export
angle <- function(tracks) {
  tracks %>%
    dplyr::summarize(Track_Angle = atan2(tail(Location_Center_Y, n = 1) - Location_Center_Y[1],
      tail(Location_Center_X, n = 1) - Location_Center_X[1]))
}

#' Calculate distance traveled and the integrated distance traveled of a track object.
#'
#' @param tracks data frame with track objects
#' @return distance traveled
#' @examples
#'  data <- tibble::data_frame(
#'    Metadata_timePoint = c(1:5),
#'    Location_Center_X = c(1, 2, 3, 4, 5),
#'    Location_Center_Y = c(1, 1, 1, 1, 1),
#'    TrackObjects_Label = c(rep(1, 5))
#'  )
#'  tracks <- neutrominer::displace(data,'TrackObjects_Label')
#'  distance <-  neutrominer::distance(tracks)
#'
#' @importFrom magrittr %>%
#' @importFrom utils tail
#' @export
distance <- function(tracks) {
  tracks %>%
    dplyr::summarize(Track_Integrated_Distance_Traveled = sum(TrackObjects_Distance_Traveled, na.rm = TRUE),
      Track_Distance_Traveled = sqrt( (tail(Location_Center_Y, n = 1) - Location_Center_Y[1] )^2 +
          (tail(Location_Center_X, n = 1) - Location_Center_X[1] )^2 ))
}

#' Calculate the directionality of a track object.
#'
#' @param tracks data frame with track objects
#'
#' @return directionality
#' @examples
#'  data <- tibble::data_frame(
#'    Metadata_timePoint = c(1:5),
#'    Location_Center_X = c(1, 2, 3, 4, 5),
#'    Location_Center_Y = c(1, 1, 1, 1, 1),
#'    TrackObjects_Label = c(rep(1, 5))
#'  )
#'  tracks <- neutrominer::displace(data,'TrackObjects_Label')
#'  directionality <-  neutrominer::directionality(tracks)
#'
#' @importFrom magrittr %>%
#' @export
directionality <- function(tracks) {
  tracks %>%
    distance() %>%
    dplyr::mutate( Track_Directionality = Track_Distance_Traveled / Track_Integrated_Distance_Traveled) %>%
    dplyr::select(-Track_Distance_Traveled, -Track_Integrated_Distance_Traveled)
}

#' Calculate the mean squared displacement of a track object.
#'
#' @param tracks data frame with track objects
#' @param tau delta t
#' @return mean_squared_displacement
#' @examples
#'  data <- tibble::data_frame(
#'    Metadata_timePoint = c(1:5),
#'    Location_Center_X = c(1, 2, 3, 4, 5),
#'    Location_Center_Y = c(1, 1, 1, 1, 1),
#'    TrackObjects_Label = c(rep(1, 5))
#'  )
#'  tau <- 2
#'  tracks <- neutrominer::displace(data,'TrackObjects_Label')
#'  mean_squared_displacement <-  neutrominer::mean_squared_displacement(tracks,tau)
#'
#' @importFrom magrittr %>%
#' @export
mean_squared_displacement <- function(tracks,tau) {
  tracks %>%
    dplyr::summarize(Track_MSD = (Location_Center_X[tau] - Location_Center_X[1])^2 +
        (Location_Center_Y[tau] - Location_Center_Y[1])^2)
}

#' Calculate the mean directional_persistence of a track object.
#'
#' @param tracks data frame with track objects
#' @return directional persistence
#' @examples
#'  data <- tibble::data_frame(
#'    Metadata_timePoint = c(1:5),
#'    Location_Center_X = c(1, 2, 3, 4, 5),
#'    Location_Center_Y = c(1, 1, 1, 1, 1),
#'    TrackObjects_Label = c(rep(1, 5))
#'  )
#'  tracks <- neutrominer::displace(data,'TrackObjects_Label')
#'  directional_persistence <-  neutrominer::directional_persistence(tracks)
#'
#' @importFrom magrittr %>%
#' @export
directional_persistence <- function(tracks) {
  directional_persistence <- tracks %>%
    directionality %>%
    dplyr::mutate(Track_DP = ceiling(3 * Track_Directionality)) %>%
    dplyr::select(-Track_Directionality)
}

#' Calculate the mean chemotaxis index of a track object.
#'
#' @param tracks data frame with track objects
#' @return chemotaxis_index
#' @examples
#' data <- tibble::data_frame(
#'   Metadata_timePoint = c(1:5),
#'   Location_Center_X = c(1, 2, 3, 4, 5),
#'   Location_Center_Y = c(1, 1, 1, 1, 1),
#'   TrackObjects_Label = c(rep(1, 5))
#' )
#'  tracks <- neutrominer::displace(data,'TrackObjects_Label')
#'  chemotaxis_index <-  neutrominer::chemotaxis_index(tracks)
#' @importFrom magrittr %>%
#' @export
chemotaxis_index <- function(tracks) {
  chemotaxis_index <- tracks %>%
    angle() %>%
    dplyr::mutate(Track_CI = -cos(Track_Angle) ) %>%
    dplyr::select(-Track_Angle)
}

#
#   \  3  /
# 1   \ /   2
#     / \
#   /  4  \
#' perform sector analysis and label each track according to its direction of movement.
#
#' @param tracks data frame with track objects
#' @return sector
#' @examples
#'  data <- tibble::data_frame(
#'    Metadata_timePoint = c(1:5),
#'    Location_Center_X = c(1, 2, 3, 4, 5),
#'    Location_Center_Y = c(1, 1, 1, 1, 1),
#'    TrackObjects_Label = c(rep(1, 5))
#'  )
#'  tracks <- neutrominer::displace(data,'TrackObjects_Label')
#'  sector_analysis <-  neutrominer::sector_analysis(tracks)
#'
#' @importFrom magrittr %>%
#' @export
sector_analysis <- function(tracks) {
  sector_analysis <- tracks %>%
    angle() %>%
    dplyr::mutate(Track_Positive_Sector     = as.numeric( abs(Track_Angle) > (3 * pi / 4)),
      Track_Negative_Sector     = as.numeric( abs(Track_Angle) < pi / 4),
      Track_Neutral_Sector_Up   = as.numeric( (Track_Angle >= pi / 4) & (Track_Angle < 3 * pi / 4 )),
      Track_Neutral_Sector_Down = as.numeric( (Track_Angle <= -pi / 4) & (Track_Angle >= -3 * pi / 4 ))
    ) %>%
    dplyr::mutate(Track_Sector = Track_Positive_Sector + 2 * Track_Negative_Sector + 3 * Track_Neutral_Sector_Up + 4 * Track_Neutral_Sector_Down) %>%
    dplyr::select(-Track_Angle)
}


#' calculate valid observation time as sum of the length of all
#' valid tracks divided by the sum of the length of all tracks
#'
#' @param tracks data frame with track objects
#' @param min_path_length minimum length of a valid track
#' @return valid_observation_time
#' @examples
#'  data <- tibble::data_frame(
#'    Metadata_timePoint = c(1:5),
#'    Location_Center_X = c(1, 2, 3, 4, 5),
#'    Location_Center_Y = c(1, 1, 1, 1, 1),
#'    TrackObjects_Label = c(rep(1, 5))
#'  )
#'  data <- dplyr::group_by_(data,'TrackObjects_Label')
#'  tracks <- neutrominer::track(data,'TrackObjects_Label')
#'  min_path_length <- 5
#'  vot <-   neutrominer::valid_observation_time(tracks, min_path_length)
#' @importFrom magrittr %>%
#' @export
valid_observation_time <- function(tracks, min_path_length) {
  valid_observation_time <- merge(tracks %>%
      dplyr::filter(Track_Length > min_path_length) %>%
      dplyr::summarise(sum_track_valid = sum(Track_Length)) ,
    tracks %>%
      dplyr::summarise(sum_track = sum(Track_Length))) %>%
    dplyr::mutate(VOT = (sum_track_valid / sum_track)) %>%
    dplyr::select(-sum_track_valid)
}


#' Identify valid tracks. Valid tracks are defined as tracks with a life time longer then a predefined value.
#'
#' @param tracks data frame with track objects
#' @param min_path_length minimum length of a valid track
#' @return valid_observation_time
#' @examples
#'  data <- tibble::data_frame(
#'    Metadata_timePoint = c(1:5),
#'    Location_Center_X = c(1, 2, 3, 4, 5),
#'    Location_Center_Y = c(1, 1, 1, 1, 1),
#'    TrackObjects_Label = c(rep(1, 5))
#'  )
#'  data <- dplyr::group_by_(data,'TrackObjects_Label')
#'  tracks <- neutrominer::track(data,'TrackObjects_Label')
#'  min_path_length <- 5
#'  validate_tracks <-   neutrominer::validate_tracks(tracks, min_path_length)
#' @importFrom magrittr %>%
#' @export
validate_tracks <- function(tracks, min_path_length){
  tracks %>%
    dplyr::mutate(Track_Valid = as.numeric(Track_Length > min_path_length)) %>%
    dplyr::summarize(
      Exp_Tracks = n(),
      Exp_Valid_Tracks = sum(Track_Valid),
      Exp_Valid_Track_Fraction = sum(Track_Valid) / n(),
      Exp_Mean_Track_Length = mean(Track_Length),
      Exp_Mean_Track_Life_Time = mean(Track_Life_Time)
      )
}


#' Assess track quality.
#' @param tracks data frame with track objects
#' @param min_path_length minimum length of a valid track
#' @param strata column name of track index column
#' @return valid_observation_time
#' @examples
#'  data <- tibble::data_frame(
#'    Metadata_timePoint = c(1:5),
#'    Location_Center_X = c(1, 2, 3, 4, 5),
#'    Location_Center_Y = c(1, 1, 1, 1, 1),
#'    TrackObjects_Label = c(rep(1, 5))
#'  )
#'  strata <- 'TrackObjects_Label'
#'  data <- dplyr::group_by_(data,strata)
#'  tracks <- neutrominer::track(data,strata)
#'  min_path_length <- 5
#'  trackQuality <- neutrominer::assess(tracks,min_path_length,strata)
#' @importFrom magrittr %>%
#' @export
assess <- function(tracks, min_path_length, strata) {
  track_info <- list(valid_observation_time(tracks, min_path_length),
    validate_tracks(tracks,min_path_length))
  return(Reduce(function(...) merge(..., all = TRUE, by_ = strata), track_info))
}
