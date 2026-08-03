# ============================================================
# model_randomforest.R
# Random Forest (randomForest) model estimation, nowcasting,
# SHAP decomposition, and output generation. ntree/mtry/nodesize
# tuned once via grid-search cross-validation on the full sample;
# the forest is refit at each OOS step with those fixed
# hyperparameters.
# ============================================================


# ---- Internal: filter candidate RHS vars for one Random Forest fit --
# Mirrors filter_glmnet_rhs() in model_elasticnet.R / filter_xgb_rhs() in
# model_xgboost.R and the reference vintage scripts' defensive checks: a
# candidate must be observed in the target row, have enough overlapping
# history with dep_var, and have non-zero variance in that history.
#' @noRd
filter_rf_rhs <- function(candidates, target_row, train_df, dep_var, min_history = 8L) {
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


#' Estimate a Random Forest model and generate nowcasts
#'
#' Runs the full Random Forest pipeline on a transformed dataset: tunes
#' \code{ntree}, \code{mtry}, and \code{nodesize} via grid-search \code{k}-fold
#' cross-validation on the full sample, fits a final
#' \code{\link[randomForest]{randomForest}} model with those hyperparameters,
#' refits an interpretable OLS on all candidate variables for reporting,
#' applies ARIMA imputation to fill predictor trailing NAs, generates an
#' expanding-window OOS nowcast, optionally computes a SHAP decomposition of
#' the full-sample fit, and saves all outputs.
#'
#' Like \code{\link{np_model_elasticnet}} and \code{\link{np_model_xgboost}},
#' the forest hyperparameters are tuned \emph{once} on the full sample and
#' held fixed across the OOS window (as in the reference methodology) —
#' re-running a full cross-validated grid search at every period would be far
#' more expensive, and the tuned hyperparameters are generally more stable
#' than the fitted trees themselves across a modest expanding window. The
#' forest itself is still re-estimated at each OOS step using only data
#' strictly before the target period, so no future information leaks into a
#' forecast. Candidate predictors that are missing at the target period, lack
#' sufficient training history, or have zero variance in that history are
#' dropped from that period's fit (mirroring the reference vintage scripts'
#' defensive filtering); \code{mtry} is capped at the number of predictors
#' still available in that period's fit.
#'
#' Unlike a linear model, Random Forest has no coefficients to report zero/
#' non-zero selection for, so the interpretable OLS used for the
#' regression-table output is fit on \emph{all} \code{candidate_vars}
#' (rather than a selected subset as in \code{\link{np_model_elasticnet}}).
#' Feature importance (\%IncMSE / IncNodePurity) and, optionally, SHAP values
#' (via \pkg{fastshap}, since \pkg{randomForest} has no native TreeSHAP
#' support) are used instead to characterise which predictors drive the
#' forest's predictions.
#'
#' @param dta_trans Data frame as returned by \code{\link{np_transform_data}}.
#'   Must contain a \code{date} column of class \code{Date} and all columns
#'   named in \code{candidate_vars}.
#' @param dep_var Character. Dependent variable column name, e.g.
#'   \code{"im_SA"}.
#' @param candidate_vars Character vector of predictor column names to feed
#'   into Random Forest (AR lags, dummies, and a wide indicator pool can all
#'   be combined here). These must already be present in \code{dta_trans}.
#'   Names of the form \code{X_lagN} whose base \code{X} exists in
#'   \code{dta_trans} are auto-created. At least two are required.
#' @param param_grid Data frame of hyperparameter combinations to
#'   cross-validate over, typically built with \code{\link{expand.grid}}
#'   with columns \code{ntree}, \code{mtry}, \code{nodesize}. Defaults to
#'   \code{NULL}, which builds \code{ntree = c(100, 300, 500)},
#'   \code{mtry} spanning \code{c(2, floor(sqrt(p)), floor(p/3))} (\code{p}
#'   = number of \code{candidate_vars}, after auto-lag creation), and
#'   \code{nodesize = c(1, 3, 5)}.
#' @param nfolds Integer. Number of cross-validation folds. Defaults to
#'   \code{10L}.
#' @param min_history Integer. Minimum number of overlapping non-NA
#'   observations (with \code{dep_var}) a candidate needs in the training
#'   window to be eligible for that OOS period's fit. Defaults to \code{8L}.
#' @param oos_start Date. Start of the OOS nowcast window.
#' @param oos_end Date. End of the OOS nowcast window (may extend beyond the
#'   last observed actual, producing true nowcasts).
#' @param out_dir Character. Directory where output files are written.
#'   Created if it does not exist. Defaults to \code{NULL} (no files saved).
#' @param model_label Character. Short label used in chart titles and file
#'   names. Defaults to \code{"RandomForest"}.
#' @param compute_shap Logical. If \code{TRUE} (default) and \code{out_dir}
#'   is supplied, compute a SHAP decomposition of the full-sample fit
#'   (beeswarm summary, bar importance, latest-observation waterfall,
#'   top-5 dependence plots, and a time series of top-5 signed
#'   contributions) using the \pkg{fastshap} package. Silently skipped with
#'   a message if \pkg{fastshap} is not installed.
#' @param shap_nsim Integer. Number of Monte Carlo simulations per SHAP value
#'   passed to \code{\link[fastshap]{explain}}. Defaults to \code{100L}.
#' @param seed Integer. Random seed for CV fold assignment and SHAP
#'   simulation reproducibility. Defaults to \code{123}.
#' @param verbose Logical. Print progress messages. Defaults to \code{TRUE}.
#'
#' @return Invisibly, a named list with elements:
#' \describe{
#'   \item{\code{full_fit_rf}}{The final \code{randomForest} object fit on
#'     all available observations with the tuned hyperparameters.}
#'   \item{\code{full_fit}}{An \code{lm} object on all \code{candidate_vars},
#'     for reporting (t-stats, R²).}
#'   \item{\code{optimal_params}}{One-row data frame with the tuned
#'     \code{ntree}, \code{mtry}, \code{nodesize}, and \code{cv_mse}.}
#'   \item{\code{importance_df}}{Data frame of \code{Feature}, \code{IncMSE},
#'     \code{IncNodePurity}, sorted by \code{IncNodePurity} descending.}
#'   \item{\code{shap_importance}}{Data frame of mean absolute SHAP value
#'     per feature, sorted descending. \code{NULL} if \code{compute_shap}
#'     was \code{FALSE} or \pkg{fastshap} was unavailable.}
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
#' rf_out <- np_model_randomforest(
#'   dta_trans      = dta_trans,
#'   dep_var        = "im_SA",
#'   candidate_vars = c("im_SA_lag1", "im_SA_lag12", "crisis_dummy",
#'                      "cpi_SA", "iip_SA", "petrol_SA", "money_supply_SA"),
#'   oos_start      = as.Date("2024-01-01"),
#'   oos_end        = as.Date("2026-04-01"),
#'   out_dir        = "Outputs/RandomForest",
#'   model_label    = "Baseline"
#' )
#' rf_out$importance_df
#' }
#'
#' @export
np_model_randomforest <- function(dta_trans,
                                  dep_var,
                                  candidate_vars,
                                  param_grid   = NULL,
                                  nfolds       = 10L,
                                  min_history  = 8L,
                                  oos_start,
                                  oos_end,
                                  out_dir      = NULL,
                                  model_label  = "RandomForest",
                                  compute_shap = TRUE,
                                  shap_nsim    = 100L,
                                  seed         = 123,
                                  verbose      = TRUE) {

  if (!requireNamespace("randomForest", quietly = TRUE))
    stop("Package 'randomForest' is required for np_model_randomforest().")

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

  full_train_df <- dta_trans[, c("date", dep_var, candidate_vars), drop = FALSE]
  full_train_df <- stats::na.omit(full_train_df)

  if (nrow(full_train_df) <= length(candidate_vars) + 1)
    stop("Insufficient observations for Random Forest after listwise deletion.")

  X_train <- full_train_df[, candidate_vars, drop = FALSE]
  y_train <- full_train_df[[dep_var]]

  if (is.null(param_grid)) {
    p <- ncol(X_train)
    mtry_vals <- sort(unique(pmax(1L, pmin(p, c(2L, floor(sqrt(p)), floor(p / 3))))))
    param_grid <- expand.grid(
      ntree    = c(100, 300, 500),
      mtry     = mtry_vals,
      nodesize = c(1, 3, 5)
    )
  }

  if (verbose) message("Step 1: Full-sample Random Forest tuning (",
                        nrow(param_grid), " combination(s) x ", nfolds, "-fold CV) ...")

  set.seed(seed)
  fold_indices <- sample(rep(seq_len(nfolds), length.out = nrow(full_train_df)))

  cv_results <- do.call(rbind, lapply(seq_len(nrow(param_grid)), function(i) {
    fold_errors <- vapply(seq_len(nfolds), function(fold) {
      train_idx <- fold_indices != fold
      test_idx  <- fold_indices == fold

      fit <- tryCatch(
        randomForest::randomForest(
          x         = X_train[train_idx, , drop = FALSE],
          y         = y_train[train_idx],
          ntree     = param_grid$ntree[i],
          mtry      = param_grid$mtry[i],
          nodesize  = param_grid$nodesize[i],
          importance = FALSE
        ),
        error = function(e) NULL
      )
      if (is.null(fit)) return(NA_real_)

      preds <- stats::predict(fit, X_train[test_idx, , drop = FALSE])
      mean((preds - y_train[test_idx])^2)
    }, numeric(1))

    data.frame(
      ntree    = param_grid$ntree[i],
      mtry     = param_grid$mtry[i],
      nodesize = param_grid$nodesize[i],
      cv_mse   = mean(fold_errors, na.rm = TRUE)
    )
  }))

  cv_results <- cv_results[is.finite(cv_results$cv_mse), ]
  if (nrow(cv_results) == 0)
    stop("Random Forest cross-validation failed for every combination in param_grid.")

  optimal_params <- cv_results[which.min(cv_results$cv_mse), ]

  if (verbose) {
    message("  Optimal ntree    = ", optimal_params$ntree)
    message("  Optimal mtry     = ", optimal_params$mtry)
    message("  Optimal nodesize = ", optimal_params$nodesize)
    message("  CV MSE           = ", round(optimal_params$cv_mse, 6))
  }

  full_fit_rf <- randomForest::randomForest(
    x          = X_train,
    y          = y_train,
    ntree      = optimal_params$ntree,
    mtry       = optimal_params$mtry,
    nodesize   = optimal_params$nodesize,
    importance = TRUE
  )

  # ---- Step 2: Feature importance -------------------------------------

  importance_mat <- randomForest::importance(full_fit_rf)
  importance_df <- data.frame(
    Feature       = rownames(importance_mat),
    IncMSE        = importance_mat[, "%IncMSE"],
    IncNodePurity = importance_mat[, "IncNodePurity"],
    stringsAsFactors = FALSE
  )
  importance_df <- importance_df[order(-importance_df$IncNodePurity), ]
  rownames(importance_df) <- NULL

  if (verbose) {
    message("  Top feature(s) by IncNodePurity:")
    for (i in seq_len(min(10, nrow(importance_df))))
      message("    ", formatC(importance_df$Feature[i], width = 30, flag = "-"),
              " IncNodePurity=", round(importance_df$IncNodePurity[i], 4))
  }

  if (!is.null(out_dir)) {
    imp_png <- file.path(out_dir, paste0("randomforest_importance_", safe_label, ".png"))
    grDevices::png(imp_png, width = 800, height = 600)
    randomForest::varImpPlot(full_fit_rf, main = paste0("Random Forest Variable Importance — ", model_label),
                             n.var = min(15, nrow(importance_df)))
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
      title     = paste0("Random Forest Model — ", model_label, " (", dep_var,
                         ", ntree=", optimal_params$ntree,
                         ", mtry=", optimal_params$mtry, ")"),
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
  # function docs); the forest is re-estimated each period on data
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

    rhs_vars <- filter_rf_rhs(candidate_vars, target_row, train_base, dep_var, min_history)

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

    X_tr <- train_df[, rhs_vars, drop = FALSE]
    y_tr <- train_df[[dep_var]]

    fit <- tryCatch(
      randomForest::randomForest(
        x          = X_tr,
        y          = y_tr,
        ntree      = optimal_params$ntree,
        mtry       = max(1, min(optimal_params$mtry, length(rhs_vars))),
        nodesize   = optimal_params$nodesize,
        importance = FALSE
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) return(na_row())

    X_new <- target_row[, rhs_vars, drop = FALSE]
    if (any(is.na(X_new))) return(na_row())

    data.frame(
      date      = t_month,
      actual    = actual,
      predicted = as.numeric(stats::predict(fit, newdata = X_new)),
      stringsAsFactors = FALSE
    )
  }))

  # Restore Date class lost by do.call(rbind) on data frames with Date columns
  oos_preds$date <- as.Date(oos_preds$date, origin = "1970-01-01")

  # ---- Step 6: Combine in-sample + OOS into nowcast_tbl --------------

  full_train_fitted <- as.numeric(stats::predict(full_fit_rf, newdata = X_train))

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
      ggplot2::geom_line(ggplot2::aes(y = nowcast, colour = "Random Forest",
                                      linetype = type),
                         linewidth = 0.8, na.rm = TRUE) +
      ggplot2::geom_point(data = oos_pts,
                          ggplot2::aes(y = predicted, colour = "Random Forest"),
                          size = 2) +
      ggplot2::geom_vline(xintercept = as.numeric(oos_start),
                          linetype = "dotdash", colour = "grey50",
                          linewidth = 0.6) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dotted",
                          colour = "grey60") +
      ggplot2::scale_colour_manual(
        values = c("Actual" = "#0C4550", "Random Forest" = "#2A9D8F")
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
        title    = paste0(dep_var, " — Random Forest (", model_label,
                          "): Actual vs Nowcast"),
        subtitle = paste0("Solid = in-sample fit  |  Dashed + points = OOS expanding-window\n",
                          "ntree=", optimal_params$ntree,
                          "  mtry=", optimal_params$mtry,
                          "  nodesize=", optimal_params$nodesize,
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
                            paste0("randomforest_nowcast_chart_", safe_label, ".png"))
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
      ggplot2::geom_line(ggplot2::aes(y = predicted, colour = "Random Forest"),
                         linewidth = 0.9, linetype = "dashed", na.rm = TRUE) +
      ggplot2::geom_point(ggplot2::aes(y = predicted, colour = "Random Forest"),
                          size = 2.5, na.rm = TRUE) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dotted",
                          colour = "grey60") +
      ggplot2::scale_colour_manual(
        values = c("Actual" = "#0C4550", "Random Forest" = "#2A9D8F")
      ) +
      ggplot2::scale_x_date(date_breaks = "3 months",
                            date_labels = "%Y-%m") +
      ggplot2::labs(
        title    = paste0(dep_var, " — Random Forest (",
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
                                paste0("randomforest_oos_eval_", safe_label, ".png"))
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
    if (!requireNamespace("fastshap", quietly = TRUE)) {
      if (verbose) message("  SHAP decomposition skipped: install 'fastshap' to enable.")
    } else if (!requireNamespace("ggplot2", quietly = TRUE)) {
      if (verbose) message("  SHAP decomposition skipped: install 'ggplot2' to enable.")
    } else {
      if (verbose) message("Step 4: SHAP decomposition (fastshap, nsim=", shap_nsim, ") ...")

      pred_wrapper <- function(object, newdata) stats::predict(object, newdata = newdata)

      set.seed(seed)
      shap_mat <- fastshap::explain(
        object       = full_fit_rf,
        X            = X_train,
        pred_wrapper = pred_wrapper,
        nsim         = shap_nsim
      )
      shap_mat <- as.matrix(shap_mat)

      shap_importance <- data.frame(
        Feature     = colnames(shap_mat),
        MeanAbsSHAP = colMeans(abs(shap_mat), na.rm = TRUE),
        stringsAsFactors = FALSE
      )
      shap_importance <- shap_importance[order(-shap_importance$MeanAbsSHAP), ]
      rownames(shap_importance) <- NULL

      top5_shap <- shap_importance$Feature[seq_len(min(5, nrow(shap_importance)))]

      # 1) Beeswarm summary plot (mean absolute SHAP per feature)
      shap_long <- data.frame(
        Feature    = rep(colnames(shap_mat), each = nrow(shap_mat)),
        SHAP_value = as.vector(shap_mat),
        stringsAsFactors = FALSE
      )
      top15 <- shap_importance$Feature[seq_len(min(15, nrow(shap_importance)))]

      p_bee <- ggplot2::ggplot(
        shap_long[shap_long$Feature %in% top15, ],
        ggplot2::aes(x = stats::reorder(Feature, SHAP_value, FUN = function(x) mean(abs(x))),
                    y = SHAP_value)
      ) +
        ggplot2::geom_jitter(width = 0.2, alpha = 0.5, size = 2, colour = "steelblue") +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
        ggplot2::coord_flip() +
        ggplot2::labs(
          title    = paste0("Random Forest SHAP Beeswarm — ", model_label),
          subtitle = "Each point is one observation; magnitude = feature impact on prediction",
          x = NULL, y = paste0("SHAP value (impact on ", dep_var, ")")
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
      ggplot2::ggsave(file.path(out_dir, paste0("randomforest_SHAP_beeswarm_", safe_label, ".png")),
                      plot = p_bee, width = 9, height = 6, dpi = 150)

      # 2) Bar chart — mean |SHAP|
      p_bar <- ggplot2::ggplot(
        shap_importance[seq_len(min(15, nrow(shap_importance))), ],
        ggplot2::aes(x = stats::reorder(Feature, MeanAbsSHAP), y = MeanAbsSHAP)
      ) +
        ggplot2::geom_col(fill = "steelblue") +
        ggplot2::coord_flip() +
        ggplot2::labs(
          title    = paste0("Random Forest Mean |SHAP| — ", model_label),
          subtitle = "Average magnitude of Shapley contribution across all training observations",
          x = NULL, y = "Mean |SHAP value|"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
      ggplot2::ggsave(file.path(out_dir, paste0("randomforest_SHAP_importance_bar_", safe_label, ".png")),
                      plot = p_bar, width = 8, height = 6, dpi = 150)

      # 3) Waterfall for the most recent training observation
      last_idx  <- nrow(X_train)
      last_shap <- shap_mat[last_idx, ]
      keep_n    <- min(15, length(last_shap))
      last_shap_sorted <- sort(abs(last_shap), decreasing = TRUE)[seq_len(keep_n)]

      waterfall_df <- data.frame(
        Feature = names(last_shap_sorted),
        SHAP    = as.numeric(last_shap[names(last_shap_sorted)]),
        stringsAsFactors = FALSE
      )
      waterfall_df <- waterfall_df[order(-waterfall_df$SHAP), ]
      waterfall_df$Feature <- factor(waterfall_df$Feature, levels = waterfall_df$Feature)

      p_wf <- ggplot2::ggplot(waterfall_df, ggplot2::aes(x = Feature, y = SHAP, fill = SHAP > 0)) +
        ggplot2::geom_col() +
        ggplot2::scale_fill_manual(values = c("FALSE" = "#d7191c", "TRUE" = "#2b83ba"), guide = "none") +
        ggplot2::coord_flip() +
        ggplot2::geom_hline(yintercept = 0, colour = "black", linewidth = 0.5) +
        ggplot2::labs(
          title    = paste0("Random Forest SHAP Waterfall — ", model_label),
          subtitle = sprintf("Most recent training observation: %s",
                             format(full_train_df$date[last_idx], "%Y-%m-%d")),
          x = NULL, y = "SHAP value (contribution to prediction)"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
      ggplot2::ggsave(file.path(out_dir, paste0("randomforest_SHAP_waterfall_", safe_label, ".png")),
                      plot = p_wf, width = 8, height = 6, dpi = 150)

      # 4) Dependence plots — top 5 features
      for (feat in top5_shap) {
        dep_df <- data.frame(
          feature_value = X_train[[feat]],
          shap_value    = shap_mat[, feat]
        )
        p_dep <- ggplot2::ggplot(dep_df, ggplot2::aes(x = feature_value, y = shap_value)) +
          ggplot2::geom_point(alpha = 0.6, colour = "steelblue", size = 2) +
          ggplot2::geom_smooth(method = "loess", colour = "red", se = FALSE, linewidth = 1,
                               formula = y ~ x) +
          ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
          ggplot2::labs(
            title    = sprintf("Random Forest SHAP Dependence: %s (%s)", feat, model_label),
            subtitle = "Red line shows LOESS trend of SHAP impact vs feature value",
            x = feat, y = sprintf("SHAP value for %s", feat)
          ) +
          ggplot2::theme_minimal() +
          ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))

        safe_feat <- gsub("[^A-Za-z0-9_]", "_", feat)
        ggplot2::ggsave(
          file.path(out_dir, sprintf("randomforest_SHAP_dependence_%s_%s.png", safe_feat, safe_label)),
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
        ggplot2::geom_point(size = 1.5, alpha = 0.7) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
        ggplot2::labs(
          title = paste0("Random Forest SHAP Contributions Over Time — ", model_label),
          x = NULL, y = "SHAP value", colour = "Feature"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "bottom")
      ggplot2::ggsave(file.path(out_dir, paste0("randomforest_SHAP_timeseries_", safe_label, ".png")),
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
                         paste0("randomforest_nowcast_", safe_label, ".xlsx"))
    openxlsx::saveWorkbook(wb, xl_path, overwrite = TRUE)
    if (verbose) message("  Saved: ", basename(xl_path))
  }

  # ---- Step 10: Save RDS -----------------------------------------------

  if (!is.null(out_dir)) {
    rds_path <- file.path(out_dir,
                          paste0("randomforest_model_data_", safe_label, ".rds"))
    saveRDS(
      list(
        dta_trans       = dta_trans,
        candidate_vars  = candidate_vars,
        dep_var         = dep_var,
        optimal_params  = optimal_params,
        cv_results      = cv_results,
        importance_df   = importance_df,
        shap_importance = shap_importance,
        full_fit_rf     = full_fit_rf,
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

  if (verbose) message("=== Random Forest model complete. ===\n")

  invisible(list(
    full_fit_rf     = full_fit_rf,
    full_fit        = full_fit,
    optimal_params  = optimal_params,
    importance_df   = importance_df,
    shap_importance = shap_importance,
    nowcast_tbl     = nowcast_tbl,
    oos_preds       = oos_preds,
    dta_trans       = dta_trans
  ))
}
