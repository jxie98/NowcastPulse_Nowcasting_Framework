test_that("summary_oos_metrics computes DirAccuracy as sign(predicted)==sign(actual)", {
  oos_preds <- data.frame(
    date      = as.Date("2020-01-01") + (0:5) * 30,
    actual    = c( 1.5, -0.5,  2.0, -1.0,  0.8, -0.2),
    predicted = c( 1.0, -0.3, -0.5, -0.4,  0.1,  0.2)  # rows 3 and 6 disagree in sign
  )

  m <- summary_oos_metrics(oos_preds)

  # sign matches for rows 1,2,4,5 (4 of 6) -> 4/6 * 100
  expect_equal(m$DirAccuracy, 100 * 4 / 6)
})

test_that("summary_oos_metrics excludes actual == 0 from DirAccuracy", {
  oos_preds <- data.frame(
    date      = as.Date("2020-01-01") + (0:3) * 30,
    actual    = c(1, -1, 0, 2),
    predicted = c(1, -1, 5, -2)   # row 3 (actual==0) excluded; row 4 disagrees
  )

  m <- summary_oos_metrics(oos_preds)

  expect_equal(m$DirAccuracy, 100 * 2 / 3)
})

test_that("summary_oos_metrics returns NA metrics (incl. DirAccuracy) when there are no valid observations", {
  oos_preds <- data.frame(date = as.Date("2020-01-01"), actual = NA_real_, predicted = NA_real_)

  m <- summary_oos_metrics(oos_preds)

  expect_equal(m$N_Obs, 0L)
  expect_true(is.na(m$DirAccuracy))
})

test_that("HitRate and DirAccuracy measure different things", {
  # actual grows every period (always positive), predicted always predicts a
  # bigger positive number each period -> HitRate should be 100 (change
  # direction always matches) but predicted also always agrees in sign with
  # actual here too, so add a case where signs disagree despite matching
  # period-over-period direction.
  oos_preds <- data.frame(
    date      = as.Date("2020-01-01") + (0:3) * 30,
    actual    = c(-3, -2, -1, 0.5),   # increasing throughout, sign flips at the end
    predicted = c(-6, -5, -4, -3)     # also increasing throughout, but always negative
  )

  m <- summary_oos_metrics(oos_preds)

  expect_equal(m$HitRate, 100)                 # every period-over-period change direction matches
  expect_lt(m$DirAccuracy, m$HitRate)           # but sign(predicted) never matches the final positive actual
})

make_oos_preds <- function(seed) {
  set.seed(seed)
  n <- 12
  data.frame(
    date      = seq(as.Date("2023-01-01"), by = "month", length.out = n),
    actual    = rnorm(n),
    predicted = rnorm(n)
  )
}

test_that("np_model_summary skips NULL model entries instead of erroring", {
  model_outputs <- list(
    Bridge = list(oos_preds = make_oos_preds(1)),
    PCA    = NULL,                                   # not run yet
    DFM    = list(oos_preds = make_oos_preds(2))
  )

  expect_message(
    result <- np_model_summary(model_outputs, dep_var = "y", verbose = TRUE),
    "Skipping model.*PCA"
  )

  expect_equal(sort(result$eval_metrics$Model), c("Bridge", "DFM"))
  expect_false("PCA" %in% names(result$oos_wide))
})

test_that("np_model_summary skips entries with invalid/missing oos_preds", {
  model_outputs <- list(
    Bridge = list(oos_preds = make_oos_preds(1)),
    XGBoost = list(oos_preds = data.frame(date = Sys.Date())),  # missing actual/predicted cols
    DFM = list(oos_preds = make_oos_preds(2))
  )

  result <- np_model_summary(model_outputs, dep_var = "y", verbose = FALSE)

  expect_equal(sort(result$eval_metrics$Model), c("Bridge", "DFM"))
})

test_that("np_model_summary errors only when every model is invalid", {
  model_outputs <- list(Bridge = NULL, PCA = NULL)

  expect_error(
    np_model_summary(model_outputs, dep_var = "y", verbose = FALSE),
    "None of the supplied model_outputs"
  )
})
