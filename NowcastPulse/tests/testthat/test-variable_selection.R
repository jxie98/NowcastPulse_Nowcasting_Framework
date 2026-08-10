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
