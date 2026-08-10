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
