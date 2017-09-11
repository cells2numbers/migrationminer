context("track")
test_that("`track` collapse single cell data to track objects", {
  # sample data set with one simple tracks
  # todo: add more tracks to cover different scenarios inlcuding
  # *invalid tracks
  data <- tibble::data_frame(
    Metadata_timePoint = c(1:5),
    Location_Center_X = c(1, 2, 3, 4, 5),
    Location_Center_Y = c(1, 1, 1, 1, 1),
    TrackObjects_Label = c(rep(1, 5))
  )

  data2 <- tibble::data_frame(
    Metadata_timePoint = c(1:5),
    Location_Center_X = c(1, 2, 3, 4, 5),
    Location_Center_Y = c(1, 1, 1, 1, 1),
    TrackObjects_Label = c(rep(1, 5)),
    Metadata_condition = c('a','a','a','a','a')
  )

  # ##############################################################################################################
  # # debug data to find inf values in directionality
  # debug_data <- tibble::data_frame(
  #   Metadata_timePoint = c(1,2,3,5,6),
  #   Location_Center_X = c(1, 2, 3, 4, 5),
  #   Location_Center_Y = c(1, 1, 1, 1, 1),
  #   TrackObjects_Label = c(2,2,2,2,2),
  #   Metadata_condition = c('a','a','a','a','a')
  # ) %>%
  #   print
  #
  # displace <- function(population, strata) {
  #   dplyr::right_join(
  #     population %>%
  #       dplyr::mutate(Helper_order =  order(Metadata_timePoint)) %>%
  #       dplyr::select_(.dots = c(strata, 'Location_Center_X', 'Location_Center_Y','Helper_order')) %>%
  #       #dplyr::mutate(Metadata_timePoint2 =  c(0,Metadata_timePoint[1:(length(Metadata_timePoint) - 1)]) ) %>%
  #       #dplyr::mutate(Metadata_timePoint2 =  (Metadata_timePoint - 1) ) %>%
  #       dplyr::mutate(Helper_order =  Helper_order - 1),
  #       #dplyr::filter(Metadata_timePoint2 != -1) %>%
  #       #dplyr::select(-Metadata_timePoint),
  #       population %>% dplyr::mutate(Helper_order =  order(Metadata_timePoint)),
  #     by = (.dots = c(strata, "Helper_order"))
  #   ) %>%
  #     dplyr::mutate(Track_dX = Location_Center_X.x - Location_Center_X.y) %>%
  #     dplyr::mutate(Track_dY = Location_Center_Y.x - Location_Center_Y.y) %>%
  #     dplyr::select(-Location_Center_X.x, -Location_Center_Y.x) %>%
  #     dplyr::rename(Location_Center_X = Location_Center_X.y) %>%
  #     dplyr::rename(Location_Center_Y = Location_Center_Y.y) %>%
  #     #dplyr::rename(Metadata_timePoint = Metadata_timePoint2) %>%
  #     dplyr::mutate(TrackObjects_Distance_Traveled = sqrt(Track_dX^2 + Track_dY^2))
  # }
  #
  # strata <- c('TrackObjects_Label', 'Metadata_condition')
  #
  # displace(debug_data, strata)

  # debug_data <- dplyr::group_by_(debug_data,.dots = c('TrackObjects_Label', 'Metadata_condition'))
  # debug_feature <- neutrominer::track(debug_data, c('TrackObjects_Label', 'Metadata_condition'))
  #
  # print(debug_feature)
  #
  # displacement <- neutrominer::displace(debug_data, c('TrackObjects_Label', 'Metadata_condition')) %>%
  # print(displacement)
  #
  # debug_data %>%
  #   dplyr::ungroup() %>%
  #   dplyr::select_(.dots = c('TrackObjects_Label', 'Location_Center_X', 'Location_Center_Y','Metadata_timePoint')) %>%
  #   dplyr::mutate(Metadata_timePoint2 =  (Metadata_timePoint - 1) ) %>%
  #   #dplyr::mutate(Metadata_timePoint3 =  c(-1,Metadata_timePoint[1:(length(Metadata_timePoint) - 1)]) ) %>%
  #   dplyr::mutate(Metadata_timePoint4 =  order(Metadata_timePoint) - 1 ) %>%
  #   dplyr::filter(Metadata_timePoint2 != -1) %>%
  #   #dplyr::select(-Metadata_condition) %>%
  #   print
  #
  #

  # Metadata_timePoint[1:length(Metadata_timePoint)]
  # end dubug code
  ##############################################################################################################

  data2 <- dplyr::group_by_(data2,.dots = c('TrackObjects_Label', 'Metadata_condition'))
  f2 <- neutrominer::track(data2, c('TrackObjects_Label', 'Metadata_condition'))

  # define results for test data
  distances <- tibble::data_frame(
    TrackObjects_Label = c(1),
    TrackObjects_Distance_Traveled = c(1,1,1,1)
  )

  track_angle <- tibble::data_frame(
    TrackObjects_Label = c(1),
    Track_Angle = c(0)
  )


  track_ci  <- tibble::data_frame(
    TrackObjects_Label = c(1),
    Track_CI = c(-1)
  )

  track_directionality  <- tibble::data_frame(
    TrackObjects_Label = c(1),
    Track_Directionality = c(1)
  )

  track_distance  <- tibble::data_frame(
    TrackObjects_Label = c(1),
    Track_Integrated_Distance_Traveled = c(4),
    Track_Distance_Traveled = c(4)
  )

  track_dp  <- tibble::data_frame(
    TrackObjects_Label = c(1),
    Track_DP = c(3)
  )

  track_fmi <- tibble::data_frame(
    TrackObjects_Label = c(1),
    Track_xFMI = c(1),
    Track_yFMI = c(0)
  )

  track_life_time  <- tibble::data_frame(
    TrackObjects_Label = c(1),
    Track_Length = c(as.integer(5)),
    Track_Life_Time = c(as.integer(5)),
    Track_One_Cell = TRUE
  )

  track_msd  <- tibble::data_frame(
    TrackObjects_Label = c(1),
    Track_MSD = c(1)
  )

  track_sectors  <- tibble::data_frame(
    TrackObjects_Label = c(1),
    Track_Positive_Sector = c(0),
    Track_Negative_Sector = c(1),
    Track_Neutral_Sector_Up = c(0),
    Track_Neutral_Sector_Down = c(0),
    Track_Sector = c(2)
  )

  track_speed  <- tibble::data_frame(
    TrackObjects_Label = c(1),
    Track_Speed = c(1),
    Track_Speed_max = c(1),
    Track_Speed_X = c(1),
    Track_Speed_Y = c(0)
  )

  vot <- tibble::data_frame(
    sum_track = as.integer(5),
    VOT = 1
  )

  valid_tracks <- tibble::data_frame(
    Exp_Tracks = 1,
    Exp_Valid_Tracks = 1,
    Exp_Valid_Track_Fraction = 1,
    Exp_Mean_Track_Length = 5,
    Exp_Mean_Track_Life_Time = 5
  )

  track_quality <- tibble::data_frame(
    sum_track = as.integer(5),
    VOT = 1,
    Exp_Tracks = 1,
    Exp_Valid_Tracks = 1,
    Exp_Valid_Track_Fraction = 1,
    Exp_Mean_Track_Length = 5,
    Exp_Mean_Track_Life_Time = 5
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
    track_speed)

  strata <- 'TrackObjects_Label'
  features <- Reduce(function(...) merge(..., all = TRUE, by = strata), feature_list)

  track_data <- neutrominer:::displace(data,strata) %>%
    dplyr::group_by_(strata)

  expect_equal(
    track_data %>%
      track(.,strata),
    features
  )

  expect_equal(
    track_data %>%
      dplyr::select(TrackObjects_Distance_Traveled) %>%
      na.omit(),
    distances
  )

  expect_equal(
    track_data %>% neutrominer::angle(.),
    track_angle
  )

  expect_equal(
    track_data %>%
      neutrominer::chemotaxis_index(.),
    track_ci
  )

  expect_equal(
    track_data %>%
      neutrominer::directionality(.),
    track_directionality
  )

  expect_equal(
    track_data %>%
      neutrominer::distance(.),
    track_distance
  )

  expect_equal(
    track_data %>%
      neutrominer::directional_persistence(.),
    track_dp
  )

  expect_equal(
    track_data %>%
      neutrominer::forward_migration_index(.),
    track_fmi
  )

  expect_equal(
    track_data %>%
      neutrominer::lifetime(.),
    track_life_time
  )

  expect_equal(
    track_data %>%
      neutrominer::mean_squared_displacement(.,2),
    track_msd
  )

  expect_equal(
    track_data %>%
      neutrominer::sector_analysis(),
    track_sectors
  )

  expect_equal(
    track_data %>%
      neutrominer::speed(),
    track_speed
  )

  expect_equivalent(
    features %>%
      neutrominer::valid_observation_time(.,3),
    vot
  )

  expect_equivalent(
    f2 %>%
      neutrominer::validate_tracks(.,2),
    valid_tracks
  )

  expect_equivalent(
    f2 %>% neutrominer::assess(.,2),
    track_quality
  )

})
