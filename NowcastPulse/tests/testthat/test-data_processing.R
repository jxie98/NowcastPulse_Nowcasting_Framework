test_that("parse_date_by_type('weekly') parses full M/D/YYYY dates correctly", {
  x <- c("1/2/2010", "1/9/2010", "1/16/2010", "1/23/2010", "1/30/2010")
  out <- parse_date_by_type(x, "weekly")

  expect_s3_class(out, "Date")
  expect_equal(out, as.Date(c("2010-01-02", "2010-01-09", "2010-01-16",
                               "2010-01-23", "2010-01-30")))
  # regression guard: these must NOT be scattered across decades
  expect_equal(format(out, "%Y"), rep("2010", 5))
})

test_that("parse_date_by_type('weekly') still supports legacy 2-token 'M/YY' format", {
  out <- parse_date_by_type(c("1/24", "12/23"), "weekly")
  expect_equal(out, as.Date(c("2024-01-01", "2023-12-01")))
})

test_that("adjust_seasonal_align falls back to unadjusted series when seasonal::final() returns numeric(0)", {
  n <- 60
  dates <- seq(as.Date("2018-01-01"), by = "month", length.out = n)
  set.seed(1)
  x <- 100 + cumsum(rnorm(n)) + 5 * sin(seq_len(n) * 2 * pi / 12)

  testthat::local_mocked_bindings(
    seas  = function(...) structure(list(), class = "seas"),
    .package = "seasonal"
  )
  testthat::local_mocked_bindings(
    final = function(...) numeric(0),
    .package = "seasonal"
  )

  out <- adjust_seasonal_align(x, dates = dates, freq = 12)

  expect_equal(out, x)
})

test_that("adjust_seasonal_align refits with a standard airline model when SEATS selects (0 0 0)", {
  n <- 60
  dates <- seq(as.Date("2018-01-01"), by = "month", length.out = n)
  set.seed(1)
  x <- 100 + cumsum(rnorm(n)) + 5 * sin(seq_len(n) * 2 * pi / 12)
  sa_values <- x - 1

  calls <- character(0)
  testthat::local_mocked_bindings(
    seas = function(..., arima.model = NULL) {
      calls <<- c(calls, if (is.null(arima.model)) NA_character_ else arima.model)
      structure(list(), class = "seas")
    },
    .package = "seasonal"
  )
  testthat::local_mocked_bindings(
    udg = function(...) "(0 0 0)(0 0 1)",
    .package = "seasonal"
  )
  testthat::local_mocked_bindings(
    final = function(...) sa_values,
    .package = "seasonal"
  )

  out <- adjust_seasonal_align(x, dates = dates, freq = 12)

  expect_equal(out, sa_values)
  expect_equal(length(calls), 2)
  expect_true(is.na(calls[1]))
  expect_equal(calls[2], "(0 1 1)(0 1 1)")
})

test_that("adjust_seasonal_align returns the seasonally adjusted series on success", {
  n <- 60
  dates <- seq(as.Date("2018-01-01"), by = "month", length.out = n)
  set.seed(1)
  x <- 100 + cumsum(rnorm(n)) + 5 * sin(seq_len(n) * 2 * pi / 12)
  sa_values <- x - 1 # arbitrary stand-in "adjusted" series of the right length

  testthat::local_mocked_bindings(
    seas  = function(...) structure(list(), class = "seas"),
    .package = "seasonal"
  )
  testthat::local_mocked_bindings(
    final = function(...) sa_values,
    .package = "seasonal"
  )

  out <- adjust_seasonal_align(x, dates = dates, freq = 12)

  expect_equal(out, sa_values)
})

test_that("week_start() anchors dates to the Sunday on/before each date", {
  d <- as.Date(c("2010-01-02", "2010-01-03", "2010-01-06", "2010-01-09"))
  # Jan 2 2010 = Saturday, Jan 3 = Sunday, Jan 6 = Wednesday, Jan 9 = Saturday
  expect_equal(
    week_start(d),
    as.Date(c("2009-12-27", "2010-01-03", "2010-01-03", "2010-01-03"))
  )
})

test_that("aggregate_to_weekly_internal buckets daily data into Sunday-anchored weeks", {
  df <- data.frame(
    date = c("20100103", "20100104", "20100108", "20100109", "20100110"),
    x    = c(1, 2, 3, 10, 20),
    stringsAsFactors = FALSE
  )
  # Jan 3 (Sun), 4 (Mon), 8 (Fri) -> week of 2010-01-03
  # Jan 9 (Sat) -> week of 2010-01-03 too (same week, Sun 3 - Sat 9)
  # Jan 10 (Sun) -> week of 2010-01-10
  out <- aggregate_to_weekly_internal(df, "daily", agg_map = c(x = "Sum"))

  expect_equal(out$date, c("2010-01-03", "2010-01-10"))
  expect_equal(out$x, c(1 + 2 + 3 + 10, 20))
})

test_that("np_aggregate_weekly re-anchors native weekly data and aggregates daily down to weekly", {
  sa_data <- list(
    daily = data.frame(
      date  = c("20100103", "20100106", "20100109", "20100110"),
      var_d = c(1, 2, 3, 100),
      stringsAsFactors = FALSE
    ),
    weekly = data.frame(
      date  = c("1/9/2010", "1/16/2010"),
      var_w = c(5, 6),
      stringsAsFactors = FALSE
    )
  )

  out <- np_aggregate_weekly(sa_data, agg_map = c(var_d = "Sum"))

  expect_equal(out$daily$date, c("2010-01-03", "2010-01-10"))
  expect_equal(out$daily$var_d, c(1 + 2 + 3, 100))
  # native weekly Jan 9 (Sat) -> week start 2010-01-03; Jan 16 -> 2010-01-10
  expect_equal(out$weekly$date, c("2010-01-03", "2010-01-10"))
  expect_null(out$daily_oil)
})

test_that("np_merge_monthly works generically to merge a weekly-aggregated list", {
  weekly_list <- list(
    daily  = data.frame(date = c("2010-01-03", "2010-01-10"), var_d = c(1, 2), stringsAsFactors = FALSE),
    weekly = data.frame(date = c("2010-01-03", "2010-01-10"), var_w = c(5, 6), stringsAsFactors = FALSE)
  )

  out <- np_merge_monthly(weekly_list, start_date = "2010-01-01")

  expect_equal(names(out), c("date", "var_d", "var_w"))
  expect_equal(out$date, c("2010-01-03", "2010-01-10"))
})

test_that("np_process_data(target_freq = 'weekly') builds a weekly combined dataset end-to-end", {
  tmp_dir <- tempfile("npdata_")
  raw_dir <- file.path(tmp_dir, "Raw")
  dir.create(raw_dir, recursive = TRUE)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  # sa_index != "NSA" so apply_x13() passes through unchanged (no X-13 call needed)
  write.csv(
    data.frame(variable_names = c("x", "y"), agg_index = c("Sum", "Average"),
              sa_index = c("SA", "SA"), trans_index = c("none", "none")),
    file.path(raw_dir, "Test_Index.csv"), row.names = FALSE
  )

  daily_dates <- format(seq(as.Date("2020-01-06"), by = "day", length.out = 14), "%Y%m%d")
  daily_df <- do.call(data.frame, c(
    list(variable_names = "x", descriptions = "Daily var",
         check.names = FALSE, stringsAsFactors = FALSE),
    stats::setNames(as.list(seq_along(daily_dates)), daily_dates)
  ))
  write.csv(daily_df, file.path(raw_dir, "Test_Daily_Data.csv"), row.names = FALSE)

  weekly_dates <- c("1/12/2020", "1/19/2020")
  weekly_df <- do.call(data.frame, c(
    list(variable_names = "y", descriptions = "Weekly var",
         check.names = FALSE, stringsAsFactors = FALSE),
    stats::setNames(as.list(c(100, 200)), weekly_dates)
  ))
  write.csv(weekly_df, file.path(raw_dir, "Test_Weekly_Data.csv"), row.names = FALSE)

  result <- suppressWarnings(np_process_data(
    data_dir = tmp_dir, prefix = "Test",
    target_freq = "weekly", save_outputs = FALSE
  ))

  expect_equal(result$target_freq, "weekly")
  expect_true(all(c("x_SA", "y_SA") %in% names(result$combined)))
  # Sunday-anchored weeks: Jan 6-11 -> 2020-01-05; Jan 12-18 -> 2020-01-12; Jan 19 -> 2020-01-19
  expect_equal(result$combined$date, c("2020-01-05", "2020-01-12", "2020-01-19"))
  expect_equal(result$combined$y_SA, c(NA_real_, 100, 200))
})
