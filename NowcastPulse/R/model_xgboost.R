# ============================================================
# model_xgboost.R
# XGBoost (gradient-boosted trees) model estimation, nowcasting,
# SHAP decomposition, and output generation. eta/max_depth/
# subsample/colsample_bytree tuned once via grid-search
# cross-validation on the full sample; the booster is refit at
# each OOS step with those fixed hyperparameters.
# ============================================================


# ---- Internal: filter candidate RHS vars for one XGBoost fit --------
# Mirrors filter_glmnet_rhs() in model_elasticnet.R and the reference
# vintage scripts' defensive checks: a candidate must be observed in the
# target row, have enough overlapping history with dep_var, and have
# non-zero variance in that history.
#' @noRd
filter_xgb_rhs <- function(candidates, target_row, train_df, dep_var, min_history = 8L) {
  in_target <- vapply(candidates, function(v) !is.na(target_row[[v]][1]), logical(1))

  with_history <- vapply(candidates, function(v) {
    sum(!is.na(train_df[[dep_var]]) & !is.na(train_df[[v]])) >= min_history
  }, logical(1))

  has_variance <- vapply(candidates, function(v) {
    valid <- train_df[[v]][!is.na(train_df[[dep_var]]) & !is.na(train_df[[v]])]
    if (length(valid) < 2) return(FALSE)
    stats::var(valid, na.rm = TRUE) > 1e-10
  }, logical(1))

  candidates[in_target & with_history & has_variance]
}


#' Estimate an XGBoost model and generate nowcasts
#'
#' Runs the full XGBoost pipeline on a transformed dataset: tunes
#' \code{eta}, \code{max_depth}, \code{subsample}, and
#' \code{colsample_bytree} via grid-search \code{k}-fold cross-validation
#' on the full sample (\code{\link[xgboost]{xgb.cv}}), fits a final
#' \code{\link[xgboost]{xgb.train}} model with those hyperparameters,
#' refits an interpretable OLS on all candidate variables for reporting,
#' applies ARIMA imputation to fill predictor trailing NAs, generates an
#' expanding-window OOS nowcast, optionally computes a TreeSHAP
#' decomposition of the full-sample fit, and saves all outputs.
#'
#' Like \code{\link{np_model_elasticnet}}, the tree hyperparameters are
#' tuned \emph{once} on the full sample and held fixed across the OOS
#' window (as in the reference methodology) — re-running a full
#' cross-validated grid search at every period would be far more
#' expensive, and the tuned hyperparameters are generally more stable than
#' the fitted trees themselves across a modest expanding window. The
#' booster itself is still re-estimated at each OOS step using only data
#' strictly before the target period, so no future information leaks into
#' a forecast. Candidate predictors that are missing at the target period,
#' lack sufficient training history, or have zero variance in that history
#' are dropped from that period's fit (mirroring the reference vintage
#' scripts' defensive filtering).
#'
#' Unlike a linear model, XGBoost has no coefficients to report zero/
#' non-zero selection for, so the interpretable OLS used for the
#' regression-table output is fit on \emph{all} \code{candidate_vars}
#' (rather than a selected subset as in \code{\link{np_model_elasticnet}}).
#' Feature importance (Gain/Cover/Frequency) and, optionally, SHAP values
#' are used instead to characterise which predictors drive the booster's
#' predictions.
#'
#' @param dta_trans Data frame as returned by \code{\link{np_transform_data}}.
#'   Must contain a \code{date} column of class \code{Date} and all columns
#'   named in \code{candidate_vars}.
#' @param dep_var Character. Dependent variable column name, e.g.
#'   \code{"im_SA"}.
#' @param candidate_vars Character vector of predictor column names to feed
#'   into XGBoost (AR lags, dummies, and a wide indicator pool can all be
#'   combined here). These must already be present in \code{dta_trans}.
#'   Names of the form \code{X_lagN} whose base \code{X} exists in
#'   \code{dta_trans} are auto-created. At least two are required.
#' @param param_grid Data frame of hyperparameter combinations to
#'   cross-validate over, typically built with \code{\link{expand.grid}}
#'   with columns \code{eta}, \code{max_depth}, \code{subsample},
#'   \code{colsample_bytree}. Defaults to a 4x3x3x3 grid spanning
#'   \code{eta = c(0.01, 0.05, 0.1, 0.3)}, \code{max_depth = c(3, 5, 7)},
#'   \code{subsample = c(0.7, 0.8, 1.0)}, \code{colsample_bytree =
#'   c(0.7, 0.8, 1.0)}.
#' @param nfolds Integer. Number of cross-validation folds. Defaults to
#'   \code{10L}.
#' @param nrounds_max Integer. Maximum boosting rounds per CV fit and per
#'   final fit (subject to \code{early_stopping_rounds}). Defaults to
#'   \code{500L}.
#' @param early_stopping_rounds Integer. Stop CV early if the held-out RMSE
#'   has not improved for this many rounds. Defaults to \code{20L}.
#' @param min_history Integer. Minimum number of overlapping non-NA
#'   observations (with \code{dep_var}) a candidate needs in the training
#'   window to be eligible for that OOS period's fit. Defaults to \code{8L}.
#' @param oos_start Date. Start of the OOS nowcast window.
#' @param oos_end Date. End of the OOS nowcast window (may extend beyond the
#'   last observed actual, producing true nowcasts).
#' @param out_dir Character. Directory where output files are written.
#'   Created if it does not exist. Defaults to \code{NULL} (no files saved).
#' @param model_label Character. Short label used in chart titles and file
#'   names. Defaults to \code{"XGBoost"}.
#' @param compute_shap Logical. If \code{TRUE} (default) and \code{out_dir}
#'   is supplied, compute a TreeSHAP decomposition of the full-sample fit
#'   (beeswarm summary, bar importance, latest-observation waterfall,
#'   top-5 dependence plots, and a time series of top-5 signed
#'   contributions) using the \pkg{shapviz} package. Silently skipped with
#'   a message if \pkg{shapviz} is not installed.
#' @param seed Integer. Random seed for CV fold assignment reproducibility.
#'   Defaults to \code{123}.
#' @param verbose Logical. Print progress messages. Defaults to \code{TRUE}.
#'
#' @return Invisibly, a named list with elements:
#' \describe{
#'   \item{\code{full_fit_xgb}}{The final \code{xgb.Booster} object fit on
#'     all available observations with the tuned hyperparameters.}
#'   \item{\code{full_fit}}{An \code{lm} object on all \code{candidate_vars},
#'     for reporting (t-stats, R²).}
#'   \item{\code{optimal_params}}{One-row data frame with the tuned
#'     \code{eta}, \code{max_depth}, \code{subsample},
#'     \code{colsample_bytree}, \code{best_iteration}, and \code{cv_rmse}.}
#'   \item{\code{importance_df}}{Data frame of \code{Feature}, \code{Gain},
#'     \code{Cover}, \code{Frequency} from \code{\link[xgboost]{xgb.importance}},
#'     sorted by \code{Gain} descending.}
#'   \item{\code{shap_importance}}{Data frame of mean absolute SHAP value
#'     per feature, sorted descending. \code{NULL} if \code{compute_shap}
#'     was \code{FALSE} or \pkg{shapviz} was unavailable.}
#'   \item{\code{nowcast_tbl}}{Data frame combining in-sample fitted values
#'     and OOS expanding-window nowcasts. Columns: \code{date},
#'     \code{actual}, \code{nowcast}, \code{type} (\code{"In-sample"} or
#'     \code{"OOS"}).}
#'   \item{\code{oos_preds}}{Data frame of OOS predictions only, with columns
#'     \code{date}, \code{actual}, \code{predicted}.}
#'   \item{\code{dta_trans}}{The \code{dta_trans} data frame after ARIMA
#'     imputation of predictor trailing NAs.}
#' }
#'
#' @examples
#' \dontrun{
#' result <- np_baseline_selection(
#'   data_dir  = "Data/", prefix = "Fiji", dep_var = "im_SA",
#'   oos_start = as.Date("2020-01-01"), oos_end = as.Date("2025-12-01")
#' )
#'
#' dta_trans <- result$dta_trans
#'
#' xgb_out <- np_model_xgboost(
#'   dta_trans      = dta_trans,
#'   dep_var        = "im_SA",
#'   candidate_vars = c("im_SA_lag1", "im_SA_lag12", "crisis_dummy",
#'                      "cpi_SA", "iip_SA", "petrol_SA", "money_supply_SA"),
#'   oos_start      = as.Date("2024-01-01"),
#'   oos_end        = as.Date("2026-04-01"),
#'   out_dir        = "Outputs/XGBoost",
#'   model_label    = "Baseline"
#' )
#' xgb_out$importance_df
#' }
#'
#' @export
np_model_xgboost <- function(dta_trans,
                             dep_var,
                             candidate_vars,
                             param_grid = expand.grid(
                               eta               = c(0.01, 0.05, 0.1, 0.3),
                               max_depth         = c(3, 5, 7),
                               subsample         = c(0.7, 0.8, 1.0),
                               colsample_bytree  = c(0.7, 0.8, 1.0)
                             ),
                             nfolds                = 10L,
                             nrounds_max           = 500L,
                             early_stopping_rounds = 20L,
                             min_history           = 8L,
                             oos_start,
                             oos_end,
                             out_dir      = NULL,
                             model_label  = "XGBoost",
                             compute_shap = TRUE,
                             seed         = 123,
                             verbose      = TRUE) {

  oos_start <- as.Date(oos_start)
  oos_end   <- as.Date(oos_end)

  if (!dep_var %in% names(dta_trans))
    stop("'", dep_var, "' not found in dta_trans.")

  lag_res        <- auto_create_lags(candidate_vars, dta_trans, verbose = verbose)
  dta_trans      <- lag_res$dta_trans
  candidate_vars <- lag_res$vars

  if (length(candidate_vars) < 2)
    stop("candidate_vars must contain at least 2 variables.")

  if (!is.null(out_dir) && !dir.exists(out_dir))
    dir.create(out_dir, recursive = TRUE)

  safe_label <- gsub("[^A-Za-z0-9_]", "_", model_label)

  # ---- Step 1: Full-sample hyperparameter grid search + fit ----------

  if (verbose) message("Step 1: Full-sample XGBoost tuning (",
                        nrow(param_grid), " combination(s) x ", nfolds, "-fold CV) ...")

  full_train_df <- dta_trans[, c("date", dep_var, candidate_vars), drop = FALSE]
  full_train_df <- stats::na.omit(full_train_df)

  if (nrow(full_train_df) <= length(candidate_vars) + 1)
    stop("Insufficient observations for XGBoost after listwise deletion.")

  X_train <- as.matrix(full_train_df[, candidate_vars, drop = FALSE])
  y_train <- full_train_df[[dep_var]]
  dtrain  <- xgboost::xgb.DMatrix(data = X_train, label = y_train)

  set.seed(seed)

  cv_results <- do.call(rbind, lapply(seq_len(nrow(param_grid)), function(i) {
    params <- list(
      objective         = "reg:squarederror",
      eta               = param_grid$eta[i],
      max_depth         = param_grid$max_depth[i],
      subsample         = param_grid$subsample[i],
      colsample_bytree  = param_grid$colsample_bytree[i]
    )

    cv_fit <- tryCatch(
      xgboost::xgb.cv(
        params = params, data = dtrain, nrounds = nrounds_max, nfold = nfolds,
        early_stopping_rounds = early_stopping_rounds, verbose = 0, metrics = "rmse"
      ),
      error = function(e) NULL
    )
    if (is.null(cv_fit)) return(NULL)

    best_iter <- if (!is.null(cv_fit$best_iteration) && length(cv_fit$best_iteration) > 0) {
      cv_fit$best_iteration
    } else {
      which.min(cv_fit$evaluation_log$test_rmse_mean)
    }

    data.frame(
      eta               = param_grid$eta[i],
      max_depth         = param_grid$max_depth[i],
      subsample         = param_grid$subsample[i],
      colsample_bytree  = param_grid$colsample_bytree[i],
      best_iteration    = best_iter,
      cv_rmse           = cv_fit$evaluation_log$test_rmse_mean[best_iter]
    )
  }))

  if (is.null(cv_results) || nrow(cv_results) == 0)
    stop("XGBoost cross-validation failed for every combination in param_grid.")

  optimal_params <- cv_results[which.min(cv_results$cv_rmse), ]

  if (verbose) {
    message("  Optimal eta               = ", round(optimal_params$eta, 3))
    message("  Optimal max_depth         = ", optimal_params$max_depth)
    message("  Optimal subsample         = ", round(optimal_params$subsample, 2))
    message("  Optimal colsample_bytree  = ", round(optimal_params$colsample_bytree, 2))
    message("  Best iteration            = ", optimal_params$best_iteration)
    message("  CV RMSE                   = ", round(optimal_params$cv_rmse, 6))
  }

  params_final <- list(
    objective         = "reg:squarederror",
    eta               = optimal_params$eta,
    max_depth         = optimal_params$max_depth,
    subsample         = optimal_params$subsample,
    colsample_bytree  = optimal_params$colsample_bytree
  )

  full_fit_xgb <- xgboost::xgb.train(
    params = params_final, data = dtrain,
    nrounds = optimal_params$best_iteration, verbose = 0
  )

  # ---- Step 2: Feature importance -------------------------------------

  importance_mat <- xgboost::xgb.importance(model = full_fit_xgb)
  importance_df <- data.frame(
    Feature   = importance_mat$Feature,
    Gain      = importance_mat$Gain,
    Cover     = importance_mat$Cover,
    Frequency = importance_mat$Frequency,
    stringsAsFactors = FALSE
  )
  importance_df <- importance_df[order(-importance_df$Gain), ]
  rownames(importance_df) <- NULL

  if (verbose) {
    message("  Top feature(s) by Gain:")
    for (i in seq_len(min(10, nrow(importance_df))))
      message("    ", formatC(importance_df$Feature[i], width = 30, flag = "-"),
              " Gain=", round(importance_df$Gain[i], 4))
  }

  if (!is.null(out_dir)) {
    imp_png <- file.path(out_dir, paste0("xgboost_importance_", safe_label, ".png"))
    grDevices::png(imp_png, width = 800, height = 600)
    xgboost::xgb.plot.importance(importance_matrix = importance_mat, top_n = 15)
    grDevices::dev.off()
    if (verbose) message("  Saved: ", basename(imp_png))
  }

  # ---- Step 3: Interpretable OLS (all candidate_vars) + regression table --

  ols_fml  <- stats::as.formula(paste(dep_var, "~", paste(candidate_vars, collapse = " + ")))
  full_fit <- stats::lm(ols_fml, data = full_train_df)

  if (verbose) message("  OLS Adj R² (all candidates, for reference) = ",
                        round(summary(full_fit)$adj.r.squared, 4))

  if (!is.null(out_dir)) {
    sg_html <- file.path(out_dir, paste0("regression_", safe_label, ".html"))
    sg_png  <- file.path(out_dir, paste0("regression_", safe_label, ".png"))

    html_content <- make_reg_html(
      full_fit,
      title     = paste0("XGBoost Model — ", model_label, " (", dep_var,
                         ", eta=", round(optimal_params$eta, 3),
                         ", max_depth=", optimal_params$max_depth, ")"),
      dep_label = dep_var
    )
    writeLines(html_content, sg_html)

    if (requireNamespace("webshot2", quietly = TRUE)) {
      webshot2::webshot(sg_html, sg_png, vwidth = 500, vheight = 600, zoom = 2)
      file.remove(sg_html)
      if (verbose) message("  Saved: ", basename(sg_png))
    } else {
      if (verbose) message("  Saved: ", basename(sg_html),
                           "  (install 'webshot2' to also get PNG)")
    }
  }

  # ---- Step 4: ARIMA imputation for predictor trailing NAs -----------

  if (verbose) message("Step 2: ARIMA imputation for predictor trailing NAs ...")
  n_imputed <- 0L
  for (col in candidate_vars) {
    if (any(is.na(dta_trans[[col]]))) {
      dta_trans[[col]] <- arima_impute_col(dta_trans[[col]])
      n_imputed <- n_imputed + 1L
    }
  }
  if (verbose) message("  Imputed ", n_imputed, " column(s).")

  # ---- Step 5: OOS expanding-window nowcast --------------------------
  # Hyperparameters are held fixed at their full-sample tuned values (see
  # function docs); the booster is re-estimated each period on data
  # strictly before it, so no future information leaks into a forecast.

  if (verbose) message("Step 3: Generating OOS nowcasts (expanding window) ...")

  oos_months <- dta_trans$date[dta_trans$date >= oos_start &
                                 dta_trans$date <= oos_end]

  oos_preds <- do.call(rbind, lapply(oos_months, function(t_month) {

    actual <- dta_trans[dta_trans$date == t_month, dep_var, drop = TRUE]
    actual <- if (length(actual) > 0) actual[1] else NA_real_

    na_row <- function() data.frame(date = t_month, actual = actual,
                                    predicted = NA_real_, stringsAsFactors = FALSE)

    target_row <- dta_trans[dta_trans$date == t_month, candidate_vars, drop = FALSE]
    if (nrow(target_row) == 0) return(na_row())

    train_base <- dta_trans[dta_trans$date < t_month, c(dep_var, candidate_vars), drop = FALSE]

    rhs_vars <- filter_xgb_rhs(candidate_vars, target_row, train_base, dep_var, min_history)

    if (length(rhs_vars) == 0) {
      train_df <- stats::na.omit(train_base[, dep_var, drop = FALSE])
      if (nrow(train_df) < 3) return(na_row())
      return(data.frame(date = t_month, actual = actual,
                        predicted = mean(train_df[[dep_var]], na.rm = TRUE),
                        stringsAsFactors = FALSE))
    }

    train_df <- stats::na.omit(train_base[, c(dep_var, rhs_vars), drop = FALSE])
    min_obs  <- max(5, length(rhs_vars) + 2)
    if (nrow(train_df) < min_obs) return(na_row())

    X_tr <- as.matrix(train_df[, rhs_vars, drop = FALSE])
    y_tr <- train_df[[dep_var]]

    fit <- tryCatch({
      dtrain_tr <- xgboost::xgb.DMatrix(data = X_tr, label = y_tr)
      xgboost::xgb.train(params = params_final, data = dtrain_tr,
                         nrounds = optimal_params$best_iteration, verbose = 0)
    }, error = function(e) NULL)
    if (is.null(fit)) return(na_row())

    X_new <- as.matrix(target_row[, rhs_vars, drop = FALSE])
    if (any(is.na(X_new))) return(na_row())

    data.frame(
      date      = t_month,
      actual    = actual,
      predicted = as.numeric(stats::predict(fit, newdata = xgboost::xgb.DMatrix(data = X_new))),
      stringsAsFactors = FALSE
    )
  }))

  # Restore Date class lost by do.call(rbind) on data frames with Date columns
  oos_preds$date <- as.Date(oos_preds$date, origin = "1970-01-01")

  # ---- Step 6: Combine in-sample + OOS into nowcast_tbl --------------

  full_train_fitted <- as.numeric(stats::predict(full_fit_xgb, newdata = dtrain))

  insample_df <- data.frame(
    date    = full_train_df$date,
    actual  = full_train_df[[dep_var]],
    nowcast = full_train_fitted,
    type    = "In-sample",
    stringsAsFactors = FALSE
  )
  oos_df <- data.frame(
    date    = oos_preds$date,
    actual  = oos_preds$actual,
    nowcast = oos_preds$predicted,
    type    = "OOS",
    stringsAsFactors = FALSE
  )
  nowcast_tbl <- rbind(insample_df, oos_df)
  # Restore Date class — rbind on mixed data frames can strip it
  nowcast_tbl$date <- as.Date(nowcast_tbl$date, origin = "1970-01-01")
  nowcast_tbl <- nowcast_tbl[order(nowcast_tbl$date), ]
  rownames(nowcast_tbl) <- NULL

  if (verbose) {
    oos_complete <- oos_preds[!is.na(oos_preds$predicted), ]
    message("  OOS nowcasts produced: ", nrow(oos_complete), " of ",
            length(oos_months), " period(s)")
    message("  Last OOS nowcast:\n",
            paste(utils::capture.output(
              print(utils::tail(oos_complete, 3), row.names = FALSE)
            ), collapse = "\n"))
  }

  # ---- Step 7: Nowcast chart -------------------------------------------

  if (!is.null(out_dir) && requireNamespace("ggplot2", quietly = TRUE)) {
    chart_data <- nowcast_tbl[!is.na(nowcast_tbl$nowcast), ]
    oos_pts    <- oos_preds[!is.na(oos_preds$predicted), ]

    p <- ggplot2::ggplot(chart_data, ggplot2::aes(x = date)) +
      ggplot2::geom_line(ggplot2::aes(y = actual,  colour = "Actual"),
                         linewidth = 0.8, na.rm = TRUE) +
      ggplot2::geom_line(ggplot2::aes(y = nowcast, colour = "XGBoost",
                                      linetype = type),
                         linewidth = 0.8, na.rm = TRUE) +
      ggplot2::geom_point(data = oos_pts,
                          ggplot2::aes(y = predicted, colour = "XGBoost"),
                          size = 2) +
      ggplot2::geom_vline(xintercept = as.numeric(oos_start),
                          linetype = "dotdash", colour = "grey50",
                          linewidth = 0.6) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dotted",
                          colour = "grey60") +
      ggplot2::scale_colour_manual(
        values = c("Actual" = "#0C4550", "XGBoost" = "#D62828")
      ) +
      ggplot2::scale_linetype_manual(
        values = c("In-sample" = "solid", "OOS" = "dashed"), guide = "none"
      ) +
      ggplot2::scale_x_date(date_breaks = "6 months",
                            date_labels = "%Y-%m") +
      ggplot2::annotate("text", x = oos_start, y = Inf,
                        label = paste0("OOS: ", format(oos_start, "%Y-%m")),
                        hjust = -0.05, vjust = 1.5, size = 3,
                        colour = "grey40") +
      ggplot2::labs(
        title    = paste0(dep_var, " — XGBoost (", model_label,
                          "): Actual vs Nowcast"),
        subtitle = paste0("Solid = in-sample fit  |  Dashed + points = OOS expanding-window\n",
                          "eta=", round(optimal_params$eta, 3),
                          "  max_depth=", optimal_params$max_depth,
                          "  |  ", length(candidate_vars), " candidate(s)"),
        x = NULL, y = paste0("Annual Growth Rate (", dep_var, ")"),
        colour = NULL
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1),
        plot.subtitle   = ggplot2::element_text(size = 7, colour = "grey50")
      )

    chart_path <- file.path(out_dir,
                            paste0("xgboost_nowcast_chart_", safe_label, ".png"))
    ggplot2::ggsave(chart_path, p, width = 12, height = 6, dpi = 150)
    if (verbose) message("  Saved: ", basename(chart_path))

    # ---- OOS-only chart with evaluation metrics ------------------
    oos_eval <- oos_preds[!is.na(oos_preds$predicted) &
                            !is.na(oos_preds$actual), ]

    if (nrow(oos_eval) >= 2) {
      oos_mae  <- mean(abs(oos_eval$predicted - oos_eval$actual))
      oos_rmse <- sqrt(mean((oos_eval$predicted - oos_eval$actual)^2))
      metrics_label <- sprintf(
        "OOS MAE = %.4f  |  RMSE = %.4f  (n = %d evaluated period(s))",
        oos_mae, oos_rmse, nrow(oos_eval)
      )
    } else {
      oos_mae <- oos_rmse <- NA_real_
      metrics_label <- "No evaluated OOS periods yet (actuals not yet available)"
    }

    oos_plot_data      <- oos_preds
    oos_plot_data$date <- as.Date(oos_plot_data$date, origin = "1970-01-01")

    p_oos <- ggplot2::ggplot(oos_plot_data,
                             ggplot2::aes(x = date)) +
      ggplot2::geom_line(ggplot2::aes(y = actual,    colour = "Actual"),
                         linewidth = 0.9, na.rm = TRUE) +
      ggplot2::geom_line(ggplot2::aes(y = predicted, colour = "XGBoost"),
                         linewidth = 0.9, linetype = "dashed", na.rm = TRUE) +
      ggplot2::geom_point(ggplot2::aes(y = predicted, colour = "XGBoost"),
                          size = 2.5, na.rm = TRUE) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dotted",
                          colour = "grey60") +
      ggplot2::scale_colour_manual(
        values = c("Actual" = "#0C4550", "XGBoost" = "#D62828")
      ) +
      ggplot2::scale_x_date(date_breaks = "3 months",
                            date_labels = "%Y-%m") +
      ggplot2::labs(
        title    = paste0(dep_var, " — XGBoost (",
                          model_label, "): OOS Evaluation"),
        subtitle = metrics_label,
        x = NULL, y = paste0("Annual Growth Rate (", dep_var, ")"),
        colour = NULL
      ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          legend.position = "bottom",
          axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1),
          plot.subtitle   = ggplot2::element_text(size = 9, colour = "grey30",
                                                  family = "mono")
        )

    oos_chart_path <- file.path(out_dir,
                                paste0("xgboost_oos_eval_", safe_label, ".png"))
    ggplot2::ggsave(oos_chart_path, p_oos, width = 10, height = 5, dpi = 150)
    if (verbose) message("  Saved: ", basename(oos_chart_path),
                         if (!is.na(oos_mae))
                           paste0("  (MAE=", round(oos_mae, 4),
                                  ", RMSE=", round(oos_rmse, 4), ")")
                         else "  (no actuals yet)")
  }

  # ---- Step 8: SHAP decomposition (optional) --------------------------

  shap_importance <- NULL

  if (compute_shap && !is.null(out_dir)) {
    if (!requireNamespace("shapviz", quietly = TRUE)) {
      if (verbose) message("  SHAP decomposition skipped: install 'shapviz' to enable.")
    } else {
      if (verbose) message("Step 4: TreeSHAP decomposition ...")

      X_train_df <- as.data.frame(X_train)
      sv <- shapviz::shapviz(full_fit_xgb, X_pred = dtrain, X = X_train_df)

      shap_raw <- stats::predict(full_fit_xgb, newdata = dtrain, predcontrib = TRUE)
      shap_mat <- shap_raw[, seq_len(ncol(X_train)), drop = FALSE]
      colnames(shap_mat) <- colnames(X_train)

      shap_importance <- data.frame(
        Feature     = colnames(shap_mat),
        MeanAbsSHAP = colMeans(abs(shap_mat)),
        stringsAsFactors = FALSE
      )
      shap_importance <- shap_importance[order(-shap_importance$MeanAbsSHAP), ]
      rownames(shap_importance) <- NULL

      top5_shap <- shap_importance$Feature[seq_len(min(5, nrow(shap_importance)))]

      # 1) Beeswarm summary
      p_bee <- shapviz::sv_importance(sv, kind = "beeswarm", max_display = 15) +
        ggplot2::labs(title = paste0("XGBoost SHAP Beeswarm — ", model_label)) +
        ggplot2::theme_minimal()
      ggplot2::ggsave(file.path(out_dir, paste0("xgboost_SHAP_beeswarm_", safe_label, ".png")),
                      plot = p_bee, width = 9, height = 6, dpi = 150)

      # 2) Bar chart — mean |SHAP|
      p_bar <- shapviz::sv_importance(sv, kind = "bar", max_display = 15) +
        ggplot2::labs(title = paste0("XGBoost Mean |SHAP| — ", model_label)) +
        ggplot2::theme_minimal()
      ggplot2::ggsave(file.path(out_dir, paste0("xgboost_SHAP_importance_bar_", safe_label, ".png")),
                      plot = p_bar, width = 8, height = 6, dpi = 150)

      # 3) Waterfall for the most recent training observation
      last_idx <- nrow(X_train_df)
      p_wf <- shapviz::sv_waterfall(sv, row_id = last_idx) +
        ggplot2::labs(
          title    = paste0("XGBoost SHAP Waterfall — ", model_label),
          subtitle = sprintf("Most recent training observation: %s",
                             format(full_train_df$date[last_idx], "%Y-%m-%d"))
        ) +
        ggplot2::theme_minimal()
      ggplot2::ggsave(file.path(out_dir, paste0("xgboost_SHAP_waterfall_", safe_label, ".png")),
                      plot = p_wf, width = 8, height = 6, dpi = 150)

      # 4) Dependence plots — top 5 features
      for (feat in top5_shap) {
        p_dep <- shapviz::sv_dependence(sv, v = feat) +
          ggplot2::labs(title = sprintf("XGBoost SHAP Dependence: %s (%s)", feat, model_label)) +
          ggplot2::theme_minimal()
        safe_feat <- gsub("[^A-Za-z0-9_]", "_", feat)
        ggplot2::ggsave(
          file.path(out_dir, sprintf("xgboost_SHAP_dependence_%s_%s.png", safe_feat, safe_label)),
          plot = p_dep, width = 7, height = 5, dpi = 150
        )
      }

      # 5) Time series of signed SHAP contributions — top 5 features
      shap_ts <- as.data.frame(shap_mat)
      shap_ts$date <- full_train_df$date
      shap_ts <- shap_ts[, c("date", top5_shap), drop = FALSE]
      shap_ts <- tidyr::pivot_longer(
        shap_ts, cols = dplyr::all_of(top5_shap),
        names_to = "Feature", values_to = "SHAP"
      )

      p_ts <- ggplot2::ggplot(shap_ts, ggplot2::aes(x = date, y = SHAP, colour = Feature)) +
        ggplot2::geom_line(linewidth = 0.9) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
        ggplot2::labs(
          title = paste0("XGBoost SHAP Contributions Over Time — ", model_label),
          x = NULL, y = "SHAP value", colour = "Feature"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "bottom")
      ggplot2::ggsave(file.path(out_dir, paste0("xgboost_SHAP_timeseries_", safe_label, ".png")),
                      plot = p_ts, width = 10, height = 6, dpi = 150)

      if (verbose) message("  Saved 5 SHAP output file group(s) to ", out_dir)
    }
  }

  # ---- Step 9: Excel export -------------------------------------------

  if (!is.null(out_dir) && requireNamespace("openxlsx", quietly = TRUE)) {
    wb <- openxlsx::createWorkbook()

    openxlsx::addWorksheet(wb, "Nowcast")
    openxlsx::writeDataTable(wb, "Nowcast", nowcast_tbl,
                             tableStyle = "TableStyleMedium9")
    openxlsx::setColWidths(wb, "Nowcast",
                           cols = seq_len(ncol(nowcast_tbl)), widths = "auto")

    latest_oos <- oos_preds[!is.na(oos_preds$predicted) &
                               oos_preds$date == max(oos_preds$date[
                                 !is.na(oos_preds$predicted)]), ]
    openxlsx::addWorksheet(wb, "Latest_Nowcast")
    openxlsx::writeDataTable(wb, "Latest_Nowcast", latest_oos,
                             tableStyle = "TableStyleMedium2")

    openxlsx::addWorksheet(wb, "Feature_Importance")
    openxlsx::writeDataTable(wb, "Feature_Importance", importance_df,
                             tableStyle = "TableStyleMedium2")

    if (!is.null(shap_importance)) {
      openxlsx::addWorksheet(wb, "SHAP_Importance")
      openxlsx::writeDataTable(wb, "SHAP_Importance", shap_importance,
                               tableStyle = "TableStyleMedium2")
    }

    openxlsx::addWorksheet(wb, "CV_Results")
    openxlsx::writeDataTable(wb, "CV_Results", cv_results,
                             tableStyle = "TableStyleMedium2")

    xl_path <- file.path(out_dir,
                         paste0("xgboost_nowcast_", safe_label, ".xlsx"))
    openxlsx::saveWorkbook(wb, xl_path, overwrite = TRUE)
    if (verbose) message("  Saved: ", basename(xl_path))
  }

  # ---- Step 10: Save RDS -----------------------------------------------

  if (!is.null(out_dir)) {
    rds_path <- file.path(out_dir,
                          paste0("xgboost_model_data_", safe_label, ".rds"))
    saveRDS(
      list(
        dta_trans       = dta_trans,
        candidate_vars  = candidate_vars,
        dep_var         = dep_var,
        optimal_params  = optimal_params,
        cv_results      = cv_results,
        importance_df   = importance_df,
        shap_importance = shap_importance,
        full_fit_xgb    = full_fit_xgb,
        full_fit        = full_fit,
        nowcast_tbl     = nowcast_tbl,
        oos_preds       = oos_preds,
        oos_start       = oos_start,
        oos_end         = oos_end
      ),
      file = rds_path
    )
    if (verbose) message("  Saved: ", basename(rds_path))
  }

  if (verbose) message("=== XGBoost model complete. ===\n")

  invisible(list(
    full_fit_xgb    = full_fit_xgb,
    full_fit        = full_fit,
    optimal_params  = optimal_params,
    importance_df   = importance_df,
    shap_importance = shap_importance,
    nowcast_tbl     = nowcast_tbl,
    oos_preds       = oos_preds,
    dta_trans       = dta_trans
  ))
}
