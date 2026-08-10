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
