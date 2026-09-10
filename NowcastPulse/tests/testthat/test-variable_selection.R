test_that("strip_variable_suffixes() strips trailing _lagN then _SA", {
  expect_equal(
    strip_variable_suffixes(c("im_SA_lag12", "cpi_SA", "REER_SA",
                              "CRISIS_DUMMY_SA", "loan_transport_SA_lag2")),
    c("im", "cpi", "REER", "CRISIS_DUMMY", "loan_transport")
  )
})

test_that("apply_trans('LOG') returns the log level, distinct from 'DLOG'", {
  x <- c(1, 2.718281828, 20.0855369)

  expect_equal(apply_trans(x, "LOG"), log(x))
  expect_equal(apply_trans(x, "log()"), log(x))
  # DLOG must still take priority over LOG (both start differently, but
  # guard against the "LOG" pattern accidentally matching "DLOG(...)")
  expect_equal(apply_trans(x, "DLOG(x)"), dlog(x))
})

test_that("apply_trans('SAAR') annualises the quarterly log-difference into a percent rate", {
  x <- c(100, 101, 102.5, 101.8, 103.2)

  d        <- c(NA_real_, diff(log(x)))
  expected <- (exp(4 * d) - 1) * 100

  expect_equal(apply_trans(x, "SAAR"), expected)
  expect_equal(apply_trans(x, "saar()"), expected)

  # Sanity check against a hand-computed value: 1% quarterly growth should
  # annualise to a bit over 4%, not exactly 4% (compounding).
  q <- c(100, 101)
  expect_equal(apply_trans(q, "SAAR")[2], (exp(4 * log(1.01)) - 1) * 100)
  expect_gt(apply_trans(q, "SAAR")[2], 4)
})

test_that("apply_trans('DIFY') returns the year-on-year (seasonal) difference", {
  x <- c(1:24)

  expected <- c(rep(NA_real_, 12), x[13:24] - x[1:12])
  expect_equal(apply_trans(x, "DIFY", n_periods = 12L), expected)
  expect_equal(apply_trans(x, "dify()", n_periods = 12L), expected)

  # DIFY must not be shadowed by the "D(" first-difference pattern.
  expect_false(isTRUE(all.equal(apply_trans(x, "DIFY", n_periods = 12L),
                                 apply_trans(x, "D(x)"))))
})

test_that("np_transform_data applies trans_index = 'LOG' as a plain log level", {
  combined <- data.frame(
    date  = c("2020-01", "2020-02", "2020-03"),
    x_SA  = c(10, 20, 40),
    stringsAsFactors = FALSE
  )

  dta_trans <- np_transform_data(
    combined    = combined,
    trans_map   = c(x = "LOG"),
    dep_var     = "x_SA",
    n_ar_lags   = 1L,
    target_freq = "monthly"
  )

  expect_equal(dta_trans$x_SA, log(c(10, 20, 40)))
})

test_that("np_transform_data supports target_freq = 'weekly'", {
  n <- 120
  # combined weekly data uses a full "YYYY-MM-DD" anchor (Sunday-start weeks)
  dates <- format(seq(as.Date("2020-01-05"), by = "week", length.out = n), "%Y-%m-%d")
  combined <- data.frame(
    date  = dates,
    im_SA = 100 + cumsum(rnorm(n)),
    stringsAsFactors = FALSE
  )

  dta_trans <- np_transform_data(
    combined    = combined,
    trans_map   = c(im = "PCHY"),
    dep_var     = "im_SA",
    n_ar_lags   = 2L,
    target_freq = "weekly"
  )

  expect_s3_class(dta_trans$date, "Date")
  expect_equal(dta_trans$date, as.Date(dates))
  expect_true(all(c("im_SA_lag1", "im_SA_lag2", "im_SA_lag52") %in% names(dta_trans)))
})

test_that("np_transform_data still handles 'YYYY-MM' anchors for non-weekly target_freq", {
  combined <- data.frame(
    date  = c("2020-01", "2020-02", "2020-03"),
    im_SA = c(100, 101, 102),
    stringsAsFactors = FALSE
  )

  dta_trans <- np_transform_data(
    combined    = combined,
    trans_map   = c(im = "PCHY"),
    dep_var     = "im_SA",
    n_ar_lags   = 1L,
    target_freq = "monthly"
  )

  expect_equal(dta_trans$date, as.Date(c("2020-01-01", "2020-02-01", "2020-03-01")))
})

test_that("read_variable_descriptions() reads descriptions from data files, Index.csv takes precedence", {
  tmp_dir <- tempfile("npdata_")
  raw_dir <- file.path(tmp_dir, "Raw")
  dir.create(raw_dir, recursive = TRUE)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  # Index.csv defines a description for "im" that should win over the data file
  write.csv(
    data.frame(
      variable_names = c("im", "cpi"),
      agg_index      = c("Average", "Average"),
      sa_index       = c("NSA", "NSA"),
      trans_index    = c("PCHY", "PCHY"),
      descriptions   = c("Imports (from index)", NA)
    ),
    file.path(raw_dir, "Fiji_Index.csv"), row.names = FALSE
  )

  write.csv(
    data.frame(
      variable_names = c("im", "cpi", "petrol"),
      descriptions   = c("Imports (from data file)", "Consumer Price Index", "Petrol price"),
      `202001` = c(1, 2, 3), `202002` = c(4, 5, 6),
      check.names = FALSE
    ),
    file.path(raw_dir, "Fiji_Monthly_Data.csv"), row.names = FALSE
  )

  desc_map <- read_variable_descriptions(tmp_dir, "Fiji")

  expect_equal(unname(desc_map["im"]),     "Imports (from index)")
  expect_equal(unname(desc_map["cpi"]),    "Consumer Price Index")
  expect_equal(unname(desc_map["petrol"]), "Petrol price")
})

test_that("np_select_variables forces base_controls into the final model and never drops them", {
  set.seed(42)
  n     <- 200
  dates <- seq(as.Date("2000-01-01"), by = "month", length.out = n)

  hf1          <- rnorm(n)             # strongly related to y -> should be selected
  hf2          <- rnorm(n)             # unrelated noise -> should not pass corr_threshold
  junk_control <- rnorm(n)             # unrelated noise, but forced as a base control
  y            <- 2 * hf1 + rnorm(n, sd = 0.1)

  dta_trans <- data.frame(
    date = dates, y = y, hf1 = hf1, hf2 = hf2, junk_control = junk_control
  )

  sel <- np_select_variables(
    dta_trans      = dta_trans,
    dep_var        = "y",
    oos_start      = dates[150],
    oos_end        = dates[n],
    corr_threshold = 0.3,
    base_controls  = "junk_control",
    verbose        = FALSE
  )

  # base_controls is forced into the final model ...
  expect_true("junk_control" %in% names(stats::coef(sel$final_model)))
  expect_true("junk_control" %in% sel$selected_vars)
  # ... even though it is unrelated to y and would normally get pruned
  final_pvals <- sel$spec_tbl$p_value[sel$spec_tbl$term == "junk_control"]
  expect_gt(final_pvals, 0.05)

  # base_controls is excluded from the correlation-ranked candidate pool
  expect_false("junk_control" %in% sel$corr_tbl$variable)

  # the real, correlated HF candidate is still found by the stepwise search
  expect_true("hf1" %in% sel$selected_vars)
})

test_that("np_select_variables errors if a base_controls column is missing", {
  set.seed(1)
  n     <- 60
  dates <- seq(as.Date("2000-01-01"), by = "month", length.out = n)
  dta_trans <- data.frame(date = dates, y = rnorm(n), hf1 = rnorm(n))

  expect_error(
    np_select_variables(
      dta_trans     = dta_trans,
      dep_var       = "y",
      oos_start     = dates[40],
      oos_end       = dates[n],
      base_controls = "does_not_exist",
      verbose       = FALSE
    ),
    "base_controls not found"
  )
})

test_that("np_load_baseline_selection() round-trips what np_baseline_selection() saves", {
  tmp_dir <- tempfile("baseline_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  fake_result <- list(
    raw       = list(trans_map = c(x = "PCHY")),
    processed = list(combined = data.frame(date = "2020-01", x = 1)),
    dta_trans = data.frame(date = as.Date("2020-01-01"), y_SA = 1, x_SA = 2),
    selection = list(selected_vars = c("x_SA"), final_model = NULL)
  )
  saveRDS(fake_result, file.path(tmp_dir, "baseline_selection_y_SA.rds"))

  loaded <- np_load_baseline_selection(tmp_dir, dep_var = "y_SA")

  expect_equal(loaded, fake_result)
})

test_that("np_load_baseline_selection() errors when no saved file exists", {
  tmp_dir <- tempfile("baseline_")
  dir.create(tmp_dir)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  expect_error(
    np_load_baseline_selection(tmp_dir, dep_var = "y_SA"),
    "No saved baseline selection found"
  )
})
