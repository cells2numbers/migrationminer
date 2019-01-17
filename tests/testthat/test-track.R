context("track")
test_that("`track` collapse single cell data to track objects", {
  # sample data set with one simple tracks
  # todo: add more tracks to cover different scenarios inlcuding
  # *invalid tracks
  data <- tibble::tibble(
    Metadata_timePoint = c(1:5),
    Location_Center_X = c(1, 2, 3, 4, 5),
    Location_Center_Y = c(1, 1, 1, 1, 1),
    TrackObjects_Label = c(rep(1, 5))
  )

  data2 <- tibble::tibble(
    Metadata_timePoint = c(1:5),
    Location_Center_X = c(1, 2, 3, 4, 5),
    Location_Center_Y = c(1, 1, 1, 1, 1),
    TrackObjects_Label = c(rep(1, 5)),
    Metadata_condition = c("a", "a", "a", "a", "a")
  )


  data2 <- dplyr::group_by_(data2,
    .dots = c("TrackObjects_Label", "Metadata_condition")
    )

  f2 <- migrationminer::track(data2,
    c("TrackObjects_Label", "Metadata_condition")
    )

  # define results for test data
  distances <- tibble::tibble(
    TrackObjects_Label = c(1),
    TrackObjects_Distance_Traveled = c(1, 1, 1, 1)
  )

  track_angle <- tibble::tibble(
    TrackObjects_Label = c(1),
    Track_Angle = c(0)
  )


  track_ci  <- tibble::tibble(
    TrackObjects_Label = c(1),
    Track_CI = c(-1)
  )

  track_directionality  <- tibble::tibble(
    TrackObjects_Label = c(1),
    Track_Directionality = c(1)
  )

  track_distance  <- tibble::tibble(
    TrackObjects_Label = c(1),
    Track_Integrated_Distance_Traveled = c(4),
    Track_Distance_Traveled = c(4)
  )

  track_dp  <- tibble::tibble(
    TrackObjects_Label = c(1),
    Track_DP = c(3)
  )

  track_fmi <- tibble::tibble(
    TrackObjects_Label = c(1),
    Track_xFMI = c(1),
    Track_yFMI = c(0)
  )

  track_life_time  <- tibble::tibble(
    TrackObjects_Label = c(1),
    Track_Length = c(as.integer(5)),
    Track_Life_Time = c(as.integer(5)),
    Track_One_Cell = TRUE
  )

  track_msd  <- tibble::tibble(
    TrackObjects_Label = c(1),
    Track_MSD = c(1)
  )

  track_sectors  <- tibble::tibble(
    TrackObjects_Label = c(1),
    Track_Positive_Sector = c(0),
    Track_Negative_Sector = c(1),
    Track_Neutral_Sector_Up = c(0),
    Track_Neutral_Sector_Down = c(0),
    Track_Sector = c(2)
  )

  track_speed  <- tibble::tibble(
    TrackObjects_Label = c(1),
    Track_Speed = c(1),
    Track_Speed_max = c(1),
    Track_Speed_X = c(1),
    Track_Speed_Y = c(0)
  )

  vot <- tibble::tibble(
    sum_track = as.integer(5),
    VOT = 1
  )

  valid_tracks <- tibble::tibble(
    Exp_Tracks = 1,
    Exp_Valid_Tracks = 1,
    Exp_Valid_Track_Fraction = 1,
    Exp_Mean_Track_Length = 5,
    Exp_Mean_Track_Life_Time = 5
  )

  track_quality <- tibble::tibble(
    sum_track = as.integer(5),
    VOT = 1,
    Exp_Tracks = 1,
    Exp_Valid_Tracks = 1,
    Exp_Valid_Track_Fraction = 1,
    Exp_Mean_Track_Length = 5,
    Exp_Mean_Track_Life_Time = 5
  )

  track_position <- tibble::tibble(
    TrackObjects_Label = c(1),
    Track_Pos_X = 3,
    Track_Pos_Y = 1
  )

  # create test data for track command
  feature_list <- list(track_angle,
    track_ci,
    track_directionality,
    track_distance,
    track_dp,
    track_fmi,
    track_life_time,
    track_msd,
    track_sectors,
    track_speed,
    track_position)

  strata <- "TrackObjects_Label"
  features <- Reduce(
    function(...) merge(..., all = TRUE, by = strata),
    feature_list
    )

  track_data <- migrationminer:::displace(data, strata) %>%
    dplyr::group_by_(strata)

  expect_equal(
    track_data %>%
      track(., strata),
    features
  )

  expect_equal(
    track_data %>%
      dplyr::select(TrackObjects_Distance_Traveled) %>%
      na.omit(),
    distances
  )

  expect_equal(
    track_data %>% migrationminer::angle(.),
    track_angle
  )

  expect_equal(
    track_data %>%
      migrationminer::chemotaxis_index(.),
    track_ci
  )

  expect_equal(
    track_data %>%
      migrationminer::directionality(.),
    track_directionality
  )

  expect_equal(
    track_data %>%
      migrationminer::distance(.),
    track_distance
  )

  expect_equal(
    track_data %>%
      migrationminer::directional_persistence(.),
    track_dp
  )

  expect_equal(
    track_data %>%
      migrationminer::forward_migration_index(.),
    track_fmi
  )

  expect_equal(
    track_data %>%
      migrationminer::lifetime(.),
    track_life_time
  )

  expect_equal(
    track_data %>%
      migrationminer::mean_squared_displacement(., 2),
    track_msd
  )

  expect_equal(
    track_data %>%
      migrationminer::sector_analysis(),
    track_sectors
  )

  expect_equal(
    track_data %>%
      migrationminer::speed(),
    track_speed
  )

  expect_equivalent(
    features %>%
      migrationminer::valid_observation_time(., 3),
    vot
  )

  expect_equivalent(
    f2 %>%
      migrationminer::validate_tracks(., 2),
    valid_tracks
  )

  expect_equivalent(
    f2 %>% migrationminer::assess(., 2),
    track_quality
  )

})
