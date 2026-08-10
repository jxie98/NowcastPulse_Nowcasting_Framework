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
