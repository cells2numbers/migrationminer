# https://github.com/tidyverse/magrittr/issues/29
#
#
if (getRversion() >= "2.15.1")  utils::globalVariables(
  c(":=", "n", "Location_Center_X", "Location_Center_Y", "Metadata_timePoint",
    "TrackObjects_Distance_Traveled", "Track_Angle", "Track_Directionality",
    "Track_Displacement_X", "Track_Displacement_Y", "Track_Distance_Traveled",
    "Track_Negative_Sector", "Track_Neutral_Sector_Down",
    "Track_Neutral_Sector_Up", "Track_Positive_Sector",  "Track_Valid",
    "Track_dX", "Track_dY", "sum_track", "sum_track_valid",
    "Track_Sector", "n_per_sector", "sector_down", "sector_left", "sector_right",
    "sector_up", "helper_order", ".")
  )

#' Compute track statistics
#'
#' @param population, single cell data
#' @param x_var variable name / columne name used for x-coordinates
#' @param y_var variable name / columne name used for y-coordinates
#' @param t_var variable name / columne name used for time coordinates
#' @param strata, column name storing the track label
#' @return track
#' @importFrom magrittr %>%
#' @examples
#'  data <- tibble::tibble(
#'    Metadata_timePoint = c(1:5),
#'    Location_Center_X = c(1, 2, 3, 4, 5),
#'    Location_Center_Y = c(1, 1, 1, 1, 1),
#'    TrackObjects_Label = c(rep(1, 5))
#'  )
#'  data <- dplyr::group_by_(data,"TrackObjects_Label")
#'  tracks <- migrationminer::track(data,"TrackObjects_Label")
#' @export
track <- function(population, strata,
                  x_var = "Location_Center_X",
                  y_var = "Location_Center_Y",
                  t_var = "Metadata_timePoint") {

  # process `population`, which is the data you get from CellProfiler
  tracks <- displace(population, strata,
    t_var = t_var,
    x_var = x_var,
    y_var = y_var)

  tracks %<>% dplyr::group_by_(.dots = strata)

  features <- list(
    angle(tracks, x_var = x_var, y_var = y_var),
    chemotaxis_index(tracks),
    directionality(tracks, x_var = x_var, y_var = y_var),
    distance(tracks, x_var = x_var, y_var = y_var),
    directional_persistence(tracks),
    forward_migration_index(
      tracks = tracks,
      x_var = x_var,
      y_var = y_var),
    lifetime(tracks, t_var = t_var),
    mean_squared_displacement(tracks, tau = 2),
    sector_analysis(tracks),
    speed(tracks),
    mean_position(tracks, x_var = x_var, y_var = y_var))

  return(Reduce(function(...) merge(..., all = TRUE, by = strata), features))
}


#' Add spatial displacement per frame for each track object
#'
#' @param population, data frame storing single cell data
#' @param x_var variable name / columne name used for x-coordinates
#' @param y_var variable name / columne name used for y-coordinates
#' @param t_var variable name / columne name used for time coordinates
#' @param strata, column name storing the track label
#' @return displacement
#' @examples
#'  data <- tibble::tibble(
#'    Metadata_timePoint = c(1:5),
#'    Location_Center_X = c(1, 2, 3, 4, 5),
#'    Location_Center_Y = c(1, 1, 1, 1, 1),
#'    TrackObjects_Label = c(rep(1, 5))
#'  )
#'  tracks <- migrationminer::displace(data,"TrackObjects_Label")
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @export
#'
displace <- function(population, strata,
                      x_var = "Location_Center_X",
                      y_var = "Location_Center_Y",
                      t_var = "Metadata_timePoint") {
  t_var <- as.name(t_var)

  displacement <- dplyr::right_join(
    population %>%
      dplyr::mutate(helper_order =  order(!! t_var)) %>%
      dplyr::select_( .dots = c(strata,  x_var, y_var, "helper_order")) %>%
      dplyr::mutate(helper_order =  helper_order - 1),
    population %>% dplyr::mutate(helper_order =  order(!! t_var)),
    by = c(strata, "helper_order")
    ) %>%
  dplyr::mutate(
    Track_dX =
      (!!(as.name(stringr::str_c(x_var, ".x")))) -
      #paste0(x_var,".x") -
      (!!(as.name(stringr::str_c(x_var, ".y"))))
  ) %>%
  dplyr::mutate(
    Track_dY =
      (!!(as.name(stringr::str_c(y_var, ".x")))) -
      (!!(as.name(stringr::str_c(y_var, ".y"))))
  ) %>%
  dplyr::select(
    - (!!(as.name(stringr::str_c(x_var, ".x")))),
    - (!!(as.name(stringr::str_c(y_var, ".x"))))
  ) %>%
    #
    # when using rename with expression "=" needs to be replaced with ":=", see
    # https://github.com/tidyverse/dplyr/issues/1600
    # also see vignette("programming")
  dplyr::rename(
    !!x_var := !!stringr::str_c(x_var, ".y")
  ) %>%
  dplyr::rename(
    !!y_var := !!stringr::str_c(y_var, ".y")
  ) %>%
  # old code without tidyeval
  #dplyr::mutate(Track_dX = Location_Center_X.x - Location_Center_X.y) %>%
  #dplyr::mutate(Track_dY = Location_Center_Y.x - Location_Center_Y.y) %>%
  #dplyr::select(-Location_Center_X.x, -Location_Center_Y.x) %>%
  #dplyr::rename(Location_Center_X = Location_Center_X.y) %>%
  #dplyr::rename(Location_Center_Y = Location_Center_Y.y) %>%
  #dplyr::rename(Metadata_timePoint = Metadata_timePoint2) %>%
  dplyr::select(-helper_order) %>%
  dplyr::mutate(
    TrackObjects_Distance_Traveled = sqrt( (Track_dX ^ 2) + (Track_dY ^ 2))
    )
}

#' Calculate speed for each track object in each frame
#'
#' @param tracks data frame with single cell data
#' @return displacement
#' @examples
#'  data <- tibble::tibble(
#'    Metadata_timePoint = c(1:5),
#'    Location_Center_X = c(1, 2, 3, 4, 5),
#'    Location_Center_Y = c(1, 1, 1, 1, 1),
#'    TrackObjects_Label = c(rep(1, 5))
#'  )
#'  tracks <- migrationminer::displace(data,"TrackObjects_Label")
#'  speed <- migrationminer::speed(tracks)
#' @importFrom magrittr %>%
#' @export
speed <- function(tracks) {
  tracks %>%
    dplyr::summarize(Track_Length = dplyr::n(),
      Track_Speed = sum(TrackObjects_Distance_Traveled, na.rm = TRUE) /
        (dplyr::n() - 1),
      Track_Speed_max = max(TrackObjects_Distance_Traveled, na.rm = TRUE),
      Track_Speed_X = sum(Track_dX, na.rm = TRUE) / (dplyr::n() - 1),
      Track_Speed_Y = sum(Track_dY, na.rm = TRUE) / (dplyr::n() - 1)) %>%
    dplyr::select(-Track_Length)
}


#' Compute the forward migration index of a tracked object
#'
#' @param tracks data frame with track objects
#' @param x_var variable name / columne name used for x-coordinates
#' @param y_var variable name / columne name used for y-coordinates
#' @return forward migration index
#' @examples
#'  data <- tibble::tibble(
#'    Metadata_timePoint = c(1:5),
#'    Location_Center_X = c(1, 2, 3, 4, 5),
#'    Location_Center_Y = c(1, 1, 1, 1, 1),
#'    TrackObjects_Label = c(rep(1, 5))
#'  )
#'  tracks <- migrationminer::displace(data,"TrackObjects_Label")
#'  forward_migration_index <- migrationminer::forward_migration_index(tracks)
#'
#' @importFrom magrittr %>%
#' @importFrom utils tail
#' @export
forward_migration_index <- function(tracks,
                                    x_var = "Location_Center_X",
                                    y_var = "Location_Center_Y") {

  x_var <- as.name(x_var)
  y_var <- as.name(y_var)

  tracks %>%
    dplyr::summarize(
      Track_Integrated_Distance_Traveled =
        sum(TrackObjects_Distance_Traveled, na.rm = TRUE),
      Track_Displacement_X =
        tail(!!x_var, n = 1) - (!!x_var)[1],
      Track_Displacement_Y =
        tail(!!y_var, n = 1) - (!!y_var)[1]
    ) %>%
    dplyr::mutate(
      Track_xFMI = Track_Displacement_X / Track_Integrated_Distance_Traveled,
      Track_yFMI = Track_Displacement_Y / Track_Integrated_Distance_Traveled) %>%
    dplyr::select(
      -Track_Integrated_Distance_Traveled,
      -Track_Displacement_X,
      -Track_Displacement_Y
      )

}


#' Lifetime of a track object.
#'
#' @param tracks data frame with track objects
#' @param t_var variable name / columne name used for time coordinates
#' @return Calculate life time of each track object
#' @examples
#'  data <- tibble::tibble(
#'    Metadata_timePoint = c(1:5),
#'    Location_Center_X = c(1, 2, 3, 4, 5),
#'    Location_Center_Y = c(1, 1, 1, 1, 1),
#'    TrackObjects_Label = c(rep(1, 5))
#'  )
#'  tracks <- migrationminer::displace(data,"TrackObjects_Label")
#'  lifetime <-  migrationminer::lifetime(tracks)
#'
#' @importFrom magrittr %>%
#' @export
lifetime  <- function(tracks,  t_var = "Metadata_timePoint") {

  t_var <- as.name(t_var)

  tracks %>%
    dplyr::summarize(
      Track_Length = dplyr::n(),
      Track_Life_Time =
        length(unique(!!t_var)),
      Track_One_Cell =
        length(unique(!!t_var)) == length(!!t_var) )
}

#' Angle of a track object.
#'
#' @param tracks data frame with track objects
#' @param x_var variable name / columne name used for x-coordinates
#' @param y_var variable name / columne name used for y-coordinates
#' @return The angle of each track
#' @examples
#'  data <- tibble::tibble(
#'    Metadata_timePoint = c(1:5),
#'    Location_Center_X = c(1, 2, 3, 4, 5),
#'    Location_Center_Y = c(1, 1, 1, 1, 1),
#'    TrackObjects_Label = c(rep(1, 5))
#'  )
#'  tracks <- migrationminer::displace(data,"TrackObjects_Label")
#'  angle <-  migrationminer::angle(tracks)
#'
#' @importFrom magrittr %>%
#' @importFrom utils tail
#' @export
angle <- function(tracks,
                  x_var = "Location_Center_X",
                  y_var = "Location_Center_Y") {

  x_var <- as.name(x_var)
  y_var <- as.name(y_var)

  tracks %>%
    dplyr::summarize(
      Track_Angle =
        atan2(tail(!!y_var, n = 1) - (!!y_var)[1],
              tail(!!x_var, n = 1) - (!!x_var)[1])
      )
}
#   tail(!!x_var, n = 1) - (!!x_var)[1],

#' Distance traveled and integrated distance traveled of a track object.
#'
#' @param tracks data frame with track objects
#' @param x_var variable name / columne name used for x-coordinates
#' @param y_var variable name / columne name used for y-coordinates
#' @return distance traveled
#' @examples
#'  data <- tibble::tibble(
#'    Metadata_timePoint = c(1:5),
#'    Location_Center_X = c(1, 2, 3, 4, 5),
#'    Location_Center_Y = c(1, 1, 1, 1, 1),
#'    TrackObjects_Label = c(rep(1, 5))
#'  )
#'  tracks <- migrationminer::displace(data,"TrackObjects_Label")
#'  distance <-  migrationminer::distance(tracks)
#'
#' @importFrom magrittr %>%
#' @importFrom utils tail
#' @export
distance <- function(tracks,
                      x_var = "Location_Center_X",
                      y_var = "Location_Center_Y") {

  x_var <- as.name(x_var)
  y_var <- as.name(y_var)

  tracks %>%
    dplyr::summarize(
      "Track_Integrated_Distance_Traveled" =
        sum(TrackObjects_Distance_Traveled, na.rm = TRUE),
      "Track_Distance_Traveled" =
        sqrt( (tail(!!y_var, n = 1) - (!!y_var)[1] ) ^ 2 +
              (tail(!!x_var, n = 1) - (!!x_var)[1] ) ^ 2 ))
}

#' Directionality of a track object.
#'
#' @param tracks data frame with track objects
#' @param x_var variable name / columne name used for x-coordinates
#' @param y_var variable name / columne name used for y-coordinates
#' @return directionality
#' @examples
#'  data <- tibble::tibble(
#'    Metadata_timePoint = c(1:5),
#'    Location_Center_X = c(1, 2, 3, 4, 5),
#'    Location_Center_Y = c(1, 1, 1, 1, 1),
#'    TrackObjects_Label = c(rep(1, 5))
#'  )
#'  tracks <- migrationminer::displace(data,"TrackObjects_Label")
#'  directionality <-  migrationminer::directionality(tracks)
#'
#' @importFrom magrittr %>%
#' @export
directionality <- function(tracks,
                            x_var = "Location_Center_X",
                            y_var = "Location_Center_Y") {
  tracks %>%
    distance(., x_var = x_var, y_var = y_var) %>%
    dplyr::mutate(
      Track_Directionality =
        Track_Distance_Traveled / Track_Integrated_Distance_Traveled) %>%
    dplyr::select(
      -Track_Distance_Traveled,
      -Track_Integrated_Distance_Traveled
      )
}

#' Mean squared displacement of a track object.
#'
#' @param tracks data frame with track objects
#' @param tau delta t
#' @return mean_squared_displacement
#' @examples
#'  data <- tibble::tibble(
#'    Metadata_timePoint = c(1:5),
#'    Location_Center_X = c(1, 2, 3, 4, 5),
#'    Location_Center_Y = c(1, 1, 1, 1, 1),
#'    TrackObjects_Label = c(rep(1, 5))
#'  )
#'  tau <- 2
#'  tracks <- migrationminer::displace(data,"TrackObjects_Label")
#'  mean_squared_displacement <-  migrationminer::mean_squared_displacement(tracks,tau)
#'
#' @importFrom magrittr %>%
#' @export
mean_squared_displacement <- function(tracks, tau = 10) {
  tracks %>%
    dplyr::summarize(Track_MSD =
        (Location_Center_X[tau] - Location_Center_X[1]) ^ 2 +
        (Location_Center_Y[tau] - Location_Center_Y[1]) ^ 2)
}

#' Calculate the mean directional_persistence of a track object.
#'
#' @param tracks data frame with track objects
#' @return directional persistence
#' @examples
#'  data <- tibble::tibble(
#'    Metadata_timePoint = c(1:5),
#'    Location_Center_X = c(1, 2, 3, 4, 5),
#'    Location_Center_Y = c(1, 1, 1, 1, 1),
#'    TrackObjects_Label = c(rep(1, 5))
#'  )
#'  tracks <- migrationminer::displace(data,"TrackObjects_Label")
#'  directional_persistence <-  migrationminer::directional_persistence(tracks)
#'
#' @importFrom magrittr %>%
#' @export
directional_persistence <- function(tracks) {
  tracks %>%
    directionality %>%
    dplyr::mutate(Track_DP = ceiling(3 * Track_Directionality)) %>%
    dplyr::select(-Track_Directionality)
}

#' Calculate the mean chemotaxis index of a track object.
#'
#' @param tracks data frame with track objects
#' @param x_var variable name / columne name used for x-coordinates
#' @param y_var variable name / columne name used for y-coordinates
#' @return chemotaxis_index
#' @examples
#' data <- tibble::tibble(
#'   Metadata_timePoint = c(1:5),
#'   Location_Center_X = c(1, 2, 3, 4, 5),
#'   Location_Center_Y = c(1, 1, 1, 1, 1),
#'   TrackObjects_Label = c(rep(1, 5))
#' )
#'  tracks <- migrationminer::displace(data,"TrackObjects_Label")
#'  chemotaxis_index <-  migrationminer::chemotaxis_index(tracks)
#' @importFrom magrittr %>%
#' @export
chemotaxis_index <- function(tracks,
                              x_var = "Location_Center_X",
                              y_var = "Location_Center_Y") {
  chemotaxis_index <- tracks %>%
    angle(., x_var = x_var, y_var = y_var) %>%
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
#' @param x_var variable name / columne name used for x-coordinates
#' @param y_var variable name / columne name used for y-coordinates
#' @return sector
#' @examples
#'  data <- tibble::tibble(
#'    Metadata_timePoint = c(1:5),
#'    Location_Center_X = c(1, 2, 3, 4, 5),
#'    Location_Center_Y = c(1, 1, 1, 1, 1),
#'    TrackObjects_Label = c(rep(1, 5))
#'  )
#'  tracks <- migrationminer::displace(data,"TrackObjects_Label")
#'  sector_analysis <-  migrationminer::sector_analysis(tracks)
#'
#' @importFrom magrittr %>%
#' @export
sector_analysis <- function(tracks,
                            x_var = "Location_Center_X",
                            y_var = "Location_Center_Y") {
  tracks %>%
    angle(., x_var = x_var, y_var = y_var) %>%
    dplyr::mutate(
      "Track_Positive_Sector" = as.numeric( abs(Track_Angle) > (3 * pi / 4)),
      "Track_Negative_Sector" = as.numeric(
        abs(Track_Angle) < pi / 4
        ),
      "Track_Neutral_Sector_Up" = as.numeric(
        (Track_Angle >= pi / 4) &
        (Track_Angle < 3 * pi / 4 )
        ),
      "Track_Neutral_Sector_Down" = as.numeric(
        (Track_Angle <= -pi / 4) &
        (Track_Angle >= -3 * pi / 4 )
        )
    ) %>%
    dplyr::mutate(Track_Sector =
        1 * Track_Positive_Sector +
        2 * Track_Negative_Sector +
        3 * Track_Neutral_Sector_Up +
        4 * Track_Neutral_Sector_Down
      ) %>%
    dplyr::select(-Track_Angle)
}


#' calculate valid observation time as sum of the length of all
#' valid tracks divided by the sum of the length of all tracks
#'
#' @param tracks data frame with track objects
#' @param min_path_length minimum length of a valid track
#' @return valid_observation_time
#' @examples
#'  data <- tibble::tibble(
#'    Metadata_timePoint = c(1:5),
#'    Location_Center_X = c(1, 2, 3, 4, 5),
#'    Location_Center_Y = c(1, 1, 1, 1, 1),
#'    TrackObjects_Label = c(rep(1, 5))
#'  )
#'  data <- dplyr::group_by_(data, "TrackObjects_Label")
#'  tracks <- migrationminer::track(data, "TrackObjects_Label")
#'  min_path_length <- 5
#'  vot <-   migrationminer::valid_observation_time(tracks, min_path_length)
#' @importFrom magrittr %>%
#' @export
valid_observation_time <- function(tracks, min_path_length = 19) {
  merge(tracks %>%
      dplyr::filter(Track_Length > min_path_length) %>%
      dplyr::summarise(sum_track_valid = sum(Track_Length)),
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
#'  data <- tibble::tibble(
#'    Metadata_timePoint = c(1:5),
#'    Location_Center_X = c(1, 2, 3, 4, 5),
#'    Location_Center_Y = c(1, 1, 1, 1, 1),
#'    TrackObjects_Label = c(rep(1, 5))
#'  )
#'  data <- dplyr::group_by_(data,"TrackObjects_Label")
#'  tracks <- migrationminer::track(data,"TrackObjects_Label")
#'  min_path_length <- 5
#'  validate_tracks <-   migrationminer::validate_tracks(tracks, min_path_length)
#' @importFrom magrittr %>%
#' @export
validate_tracks <- function(tracks, min_path_length = 19){
  tracks %>%
    dplyr::mutate(Track_Valid = as.numeric(Track_Length > min_path_length)) %>%
    dplyr::summarize(
      "Exp_Tracks" = dplyr::n(),
      "Exp_Valid_Tracks" = sum(Track_Valid),
      "Exp_Valid_Track_Fraction" = sum(Track_Valid) / dplyr::n(),
      "Exp_Mean_Track_Length" = mean(Track_Length),
      "Exp_Mean_Track_Life_Time" = mean(Track_Life_Time)
      )
}


#' Assess track quality.
#' @param tracks data frame with track objects
#' @param min_path_length minimum length of a valid track
#' @param strata column name of track index column
#' @return valid_observation_time
#' @examples
#'  data <- tibble::tibble(
#'    Metadata_timePoint = c(1:5),
#'    Location_Center_X = c(1, 2, 3, 4, 5),
#'    Location_Center_Y = c(1, 1, 1, 1, 1),
#'    TrackObjects_Label = c(rep(1, 5))
#'  )
#'  strata <- "TrackObjects_Label"
#'  data <- dplyr::group_by_(data, strata)
#'  tracks <- migrationminer::track(data, strata)
#'  min_path_length <- 5
#'  trackQuality <- migrationminer::assess(tracks,min_path_length,strata)
#' @importFrom magrittr %>%
#' @export
assess <- function(tracks, min_path_length = 19, strata) {
  track_info <- list(
    migrationminer::valid_observation_time(tracks, min_path_length),
    migrationminer::validate_tracks(tracks, min_path_length)
    )

  return(Reduce(function(...) merge(..., all = TRUE, by_ = strata), track_info))
}

#'Get mean position
#' @param tracks data frame with track objects
#' @param strata column name of track index column
#' @param x_var variable name / columne name used for x-coordinates
#' @param y_var variable name / columne name used for y-coordinates
#' @return data frame with mean x and y positions
#' @examples
#'  data <- tibble::tibble(
#'    Metadata_timePoint = c(1:5),
#'    Location_Center_X = c(1, 2, 3, 4, 5),
#'    Location_Center_Y = c(1, 1, 1, 1, 1),
#'    TrackObjects_Label = c(rep(1, 5))
#'  )
#'  strata <- 'TrackObjects_Label'
#'  data <- dplyr::group_by_(data,strata)
#'  position <- migrationminer::mean_position(data,strata)
#' @importFrom magrittr %>%
#' @export
mean_position <- function(tracks, strata,
                          x_var = "Location_Center_X",
                          y_var = "Location_Center_Y") {

  x_var <- as.name(x_var)
  y_var <- as.name(y_var)

  tracks %>% dplyr::summarise(
        "Track_Pos_X" = mean(!!x_var),
        "Track_Pos_Y" = mean(!!y_var)
    )
}


#'Summarize sector analysis per experiment
#' @param tracks data frame with track objects
#' @param strata column name of track index column
#' @return sector_summary
#' @examples
#'  data <- tibble::tibble(
#'    Metadata_timePoint = c(rep(0, 5),rep(1,5)),
#'    Location_Center_X = c(0, 0, 0, 0, 0, -1,  1,  0,  0, 0),
#'    Location_Center_Y = c(0, 0, 0, 0, 0,  0,  0,  1, -1, 1),
#'    TrackObjects_Label = c(1:5,1:4,1),
#'    Metadata_experiment = c(rep(0,10))
#'  )
#'
#'  strata <- c('TrackObjects_Label','Metadata_experiment')
#'  tracks <- track(data, strata = strata)
#'  sector <- summarize_sectors(tracks,'Metadata_experiment')
#' @importFrom magrittr %>%
#' @export
summarize_sectors <- function(tracks, strata) {
  results <- tracks %>%
    dplyr::ungroup() %>%
    dplyr::group_by_(.dots = c(strata,'Track_Sector')) %>%
    dplyr::count() %>%
    dplyr::rename(n_per_sector  = n ) %>%
    tidyr::spread(key = Track_Sector, value = n_per_sector)

  # replace na with 0
  results[is.na(results)] <- 0

  # add missing cols
  missing_cols <- setdiff(c("1","2","3","4"),setdiff(colnames(results),strata))
  results[missing_cols] <- 0

  sector_summary <- results %>%
    dplyr::rename('sector_left' = '1') %>%
    dplyr::rename('sector_right' = '2') %>%
    dplyr::rename('sector_up' = '3') %>%
    dplyr::rename('sector_down' = '4') %>%
    dplyr::mutate(sector_left_fraction =
        sector_left  / (sector_left + sector_right + sector_up + sector_down)) %>%
    dplyr::mutate(sector_right_fraction =
        sector_right / (sector_left + sector_right + sector_up + sector_down)) %>%
    dplyr::mutate(sector_up_fraction =
        sector_up    / (sector_left + sector_right + sector_up + sector_down)) %>%
    dplyr::mutate(sector_down_fraction =
        sector_down  / (sector_left + sector_right + sector_up + sector_down))
}


#'Windrose style plot
#'
#' @param tracks data frame with track objects
#' @param spdres speed resolution used for coloring
#' @param dirres angle resolution of the plot
#' @param spdmin minimum speed (used for color mapping of speed)
#' @param spdmax max speed (used for color mapping of speed )
#' @param palette color map (used for color mapping of speed)
#' @param title_name title name of the plot
#' @param scale_name legend name / title for the scale
#' @return windrose_plot
#' @examples
#'  df <- dplyr::tibble(
#'    Track_Speed = runif(1000, max = 5),
#'    Track_Angle =  runif(1000, min = -pi, max = pi)
#'  )
#'  plot_windrose(df)
#' @importFrom magrittr %>% %<>%
#' @export
plot_windrose <- function(tracks,
  spdres = 1,
  dirres = 30,
  spdmin = 0,
  spdmax = 5,
  palette = "YlGnBu",
  scale_name = "speed in pixel/frame",
  title_name = "Distribution of angle and speed"){

  tracks %<>%
    dplyr::mutate(Track_Angle = Track_Angle + pi/2) %>% # rotate all angles by 90 degree or pi/2
    dplyr::mutate(Track_Angle = ifelse(Track_Angle > pi, Track_Angle - 2*pi, Track_Angle) ) %>%
    dplyr::mutate(Track_Angle = ifelse(Track_Angle < 0 , Track_Angle + 2*pi, Track_Angle) ) %>%
    dplyr::filter(!is.na(Track_Speed), !is.na(Track_Angle))


  data <- data.frame(
    spd = tracks$Track_Speed ,
    dir = (180 * (tracks$Track_Angle) / pi)
  )

  spd = "spd"
  dir = "dir"


  # Tidy up input data ----
  n.in <- NROW(data)
  dnu <- (is.na(data[[spd]]) | is.na(data[[dir]]))
  data[[spd]][dnu] <- NA
  data[[dir]][dnu] <- NA

  # figure out the wind speed bins ----
  spdseq <- seq(spdmin,spdmax,spdres)

  # get some information about the number of bins, etc.
  n.spd.seq <- length(spdseq)
  n.colors.in.range <- n.spd.seq - 1

  # create the color map
  spd.colors <- colorRampPalette(RColorBrewer::brewer.pal(min(max(3,
    n.colors.in.range),
    min(9,
      n.colors.in.range)),
    palette))(n.colors.in.range)

  if (max(data[[spd]],na.rm = TRUE) > spdmax) {
    spd.breaks <- c(spdseq,
      max(data[[spd]],na.rm = TRUE))
    spd.labels <- c(paste(c(spdseq[1:n.spd.seq - 1]),
      '-',
      c(spdseq[2:n.spd.seq])),
      paste(spdmax,
        "-",
        max(data[[spd]],na.rm = TRUE)))
    spd.colors <- c(spd.colors, "grey50")
  } else{
    spd.breaks <- spdseq
    spd.labels <- paste(c(spdseq[1:n.spd.seq - 1]),
      '-',
      c(spdseq[2:n.spd.seq]))
  }
  data$spd.binned <- cut(x = data[[spd]],
    breaks = spd.breaks,
    labels = spd.labels,
    ordered_result = TRUE)
  # clean up the data
  data. <- na.omit(data)

  # figure out the wind direction bins
  dir.breaks <- c(-dirres/2,
    seq(dirres/2, 360 - dirres/2, by = dirres),
    360 + dirres/2)
  #
  dir.labels <- c(paste(360 - dirres/2,"-",dirres/2),
    paste(seq(dirres/2, 360 - 3*dirres/2, by = dirres),
      "-",
      seq(3*dirres/2, 360 - dirres/2, by = dirres)),
    paste(360 - dirres/2,"-",dirres/2))

  # assign each wind direction to a bin
  dir.binned <- cut(data[[dir]],
    breaks = dir.breaks,
    ordered_result = TRUE)
  levels(dir.binned) <- dir.labels
  data$dir.binned <- dir.binned


  # deal with change in ordering introduced somewhere around version 2.2
  if (packageVersion("ggplot2") > "2.2") {
    #cat("Hadley broke my code\n")
    data$spd.binned = with(data, factor(spd.binned, levels = rev(levels(spd.binned))))
    spd.colors = rev(spd.colors)
  }

  # create the plot ----
  windrose_plot <- ggplot2::ggplot(data = data,
    ggplot2::aes(x = dir.binned,
      fill = spd.binned)) +
    ggplot2::geom_bar() +
    ggplot2::scale_x_discrete(drop = FALSE,
      labels = ggplot2::waiver()) +
    ggplot2::coord_polar(start = -((dirres/2)/360) * 2*pi) +
    ggplot2::scale_fill_manual(name = scale_name,
      values = spd.colors,
      drop = FALSE) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank()) +
    ggplot2::labs(title = title_name)

  # return the handle to the wind rose
  return(windrose_plot)
}

#'3D trajectory plot
#'
#' @param tracks data frame with track objects
#' @param strata
#' trajectory_plot3d() <- function()
