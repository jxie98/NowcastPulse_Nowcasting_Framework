# ============================================================
# model_summary.R
# Cross-model comparison: consolidates OOS predictions, the
# latest nowcast, comparison charts, and an evaluation-metrics
# table across any set of np_model_*() outputs.
# ============================================================


# ---- Internal: OOS accuracy metrics for one model's oos_preds -------
#' @noRd
summary_oos_metrics <- function(oos_preds) {
  ev <- oos_preds[!is.na(oos_preds$predicted) & !is.na(oos_preds$actual), ]
  ev <- ev[order(ev$date), ]
  n  <- nrow(ev)

  if (n == 0) {
    return(data.frame(N_Obs = 0L, MAE = NA_real_, RMSE = NA_real_,
                      MAPE = NA_real_, Bias = NA_real_,
                      Correlation = NA_real_, HitRate = NA_real_,
                      DirAccuracy = NA_real_, TheilU2 = NA_real_))
  }

  err  <- ev$predicted - ev$actual
  mae  <- mean(abs(err))
  rmse <- sqrt(mean(err^2))
  mape <- if (all(ev$actual != 0)) mean(abs(err / ev$actual)) * 100 else NA_real_
  bias <- mean(err)
  corr <- if (n >= 2 && stats::sd(ev$actual) > 0 && stats::sd(ev$predicted) > 0) {
    stats::cor(ev$actual, ev$predicted)
  } else NA_real_

  hit_rate <- NA_real_
  if (n >= 2) {
    actual_diff <- diff(ev$actual)
    pred_diff   <- diff(ev$predicted)
    valid <- !is.na(actual_diff) & !is.na(pred_diff) & actual_diff != 0
    if (any(valid)) hit_rate <- mean(sign(actual_diff[valid]) == sign(pred_diff[valid])) * 100
  }

  # Directional accuracy: does the level itself (e.g. a growth/change rate)
  # have the same sign as actual (expansion vs. contraction call), as
  # opposed to HitRate above (period-over-period change in that rate).
  # Periods where actual == 0 are excluded (sign is ambiguous).
  dir_accuracy <- NA_real_
  dir_valid <- !is.na(ev$actual) & !is.na(ev$predicted) & ev$actual != 0
  if (any(dir_valid))
    dir_accuracy <- mean(sign(ev$predicted[dir_valid]) == sign(ev$actual[dir_valid])) * 100

  # Theil U2: RMSE(model) / RMSE(naive random-walk); <1 beats no-change forecast
  theil_u2 <- NA_real_
  if (n >= 2) {
    naive_err <- diff(ev$actual)          # naive: predict actual[t-1] for actual[t]
    model_err <- ev$predicted[-1] - ev$actual[-1]
    denom <- sqrt(mean(naive_err^2))
    if (!is.na(denom) && denom > 0)
      theil_u2 <- sqrt(mean(model_err^2)) / denom
  }

  data.frame(N_Obs = n, MAE = mae, RMSE = rmse, MAPE = mape,
            Bias = bias, Correlation = corr, HitRate = hit_rate,
            DirAccuracy = dir_accuracy, TheilU2 = theil_u2)
}


#' Summarise and compare OOS nowcasts across multiple models
#'
#' Consolidates the out-of-sample (OOS) output of any number of
#' \code{np_model_*()} results (Bridge, PCA, DFM, 3PRF, Elastic Net,
#' XGBoost, Random Forest, ...) into a single cross-model comparison:
#' a wide/long Excel export of actual vs. nowcast values for every OOS
#' period and model, an Excel export of each model's nowcast for the
#' next period immediately after the last observed actual, stacked and
#' per-model comparison charts, and an evaluation-metrics table (MAE,
#' RMSE, MAPE, Bias, Correlation, HitRate, DirAccuracy) with the best
#' model per metric highlighted. \code{HitRate} and \code{DirAccuracy}
#' are both directional-accuracy metrics but answer different questions:
#' \code{HitRate} is the \% of OOS periods where the period-over-period
#' \emph{change} in the predicted value moves the same way as the actual
#' change (turning points), while \code{DirAccuracy} is the \% of OOS
#' periods where \code{predicted} and \code{actual} simply have the same
#' sign (e.g. both call growth vs. both call contraction) — the more
#' relevant one when \code{dep_var} is already a growth/change rate.
#'
#' Every \code{np_model_*()} function in this package returns
#' \code{oos_preds} with the same three columns (\code{date},
#' \code{actual}, \code{predicted}), which is the only thing this
#' function requires from each entry of \code{model_outputs} — so any
#' mix of models fit on the same \code{dep_var} and any (possibly
#' different) OOS windows can be compared directly.
#'
#' @param model_outputs Named list of \code{np_model_*()} return values
#'   (e.g. \code{list(Bridge = bridge_out, ElasticNet = elasticnet_out,
#'   XGBoost = xgboost_out, RandomForest = randomforest_out)}). Names
#'   are used as model labels throughout (columns, legends, filenames)
#'   and must be unique. Each element must contain an \code{oos_preds}
#'   data frame with columns \code{date}, \code{actual}, \code{predicted}.
#'   Any element that is \code{NULL} (e.g. a model you haven't run yet) or
#'   lacks a valid \code{oos_preds} is skipped with a message rather than
#'   raising an error \u2014 the comparison proceeds with whichever models
#'   remain, as long as at least one is valid.
#' @param dep_var Character. Dependent variable name, used for axis/
#'   sheet labelling only. Defaults to \code{"dep_var"}.
#' @param out_dir Character. Directory where output files are written.
#'   Created if it does not exist. Defaults to \code{NULL} (no files
#'   saved; the function still returns its computed objects).
#' @param model_label Character. Short label used in chart titles and
#'   file names. Defaults to \code{"Comparison"}.
#' @param verbose Logical. Print progress messages. Defaults to \code{TRUE}.
#'
#' @return Invisibly, a named list with elements:
#' \describe{
#'   \item{\code{oos_wide}}{Data frame: \code{date}, \code{actual}, and
#'     one column per model with that model's OOS predicted value.}
#'   \item{\code{oos_long}}{Data frame: \code{date}, \code{model},
#'     \code{actual}, \code{predicted} — one row per model per OOS period.}
#'   \item{\code{latest_nowcast}}{Data frame: \code{model}, \code{date}
#'     (the first period after the last observed actual), \code{nowcast}.}
#'   \item{\code{last_actual_date}, \code{target_date}}{The last date
#'     with an observed actual (across all models) and the immediate
#'     next period being nowcast.}
#'   \item{\code{eval_metrics}}{Data frame of \code{Model}, \code{N_Obs},
#'     \code{MAE}, \code{RMSE}, \code{MAPE}, \code{Bias},
#'     \code{Correlation}, \code{HitRate}, \code{DirAccuracy},
#'     \code{TheilU2}, one row per model.}
#' }
#'
#' @examples
#' \dontrun{
#' summary_out <- np_model_summary(
#'   model_outputs = list(
#'     Bridge        = bridge_out,
#'     PCA           = pca_out,
#'     DFM           = dfm_out,
#'     ThreePRF      = tprf_out,
#'     ElasticNet    = elasticnet_out,
#'     XGBoost       = xgboost_out,
#'     RandomForest  = randomforest_out
#'   ),
#'   dep_var     = "im_SA",
#'   out_dir     = "Outputs/9_Summary",
#'   model_label = "Baseline"
#' )
#' summary_out$eval_metrics
#' }
#'
#' @export
np_model_summary <- function(model_outputs,
                             dep_var     = "dep_var",
                             out_dir     = NULL,
                             model_label = "Comparison",
                             verbose     = TRUE) {

  if (!is.list(model_outputs) || length(model_outputs) == 0)
    stop("model_outputs must be a non-empty named list of np_model_*() outputs.")

  model_names_all <- names(model_outputs)
  if (is.null(model_names_all) || any(model_names_all == "") || anyDuplicated(model_names_all))
    stop("model_outputs must be a named list with unique, non-empty names.")

  # Skip any model that hasn't been run yet (NULL) or lacks a valid
  # $oos_preds, rather than failing the whole comparison outright.
  is_valid <- vapply(model_names_all, function(m) {
    op <- model_outputs[[m]]$oos_preds
    !is.null(op) && all(c("date", "actual", "predicted") %in% names(op))
  }, logical(1))

  if (!all(is_valid)) {
    skipped <- model_names_all[!is_valid]
    if (verbose) message("  Skipping model(s) with no valid $oos_preds (need date/actual/predicted): ",
                          paste(skipped, collapse = ", "))
    model_outputs <- model_outputs[is_valid]
  }

  model_names <- names(model_outputs)
  if (length(model_names) == 0)
    stop("None of the supplied model_outputs have a valid $oos_preds (need date/actual/predicted).")

  if (!is.null(out_dir) && !dir.exists(out_dir))
    dir.create(out_dir, recursive = TRUE)

  safe_label <- gsub("[^A-Za-z0-9_]", "_", model_label)

  # ---- Step 1: Build long + wide OOS comparison tables ----------------

  if (verbose) message("Step 1: Consolidating OOS predictions across ",
                        length(model_names), " model(s) ...")

  oos_long <- do.call(rbind, lapply(model_names, function(m) {
    op <- model_outputs[[m]]$oos_preds
    data.frame(
      date      = as.Date(op$date, origin = "1970-01-01"),
      model     = m,
      actual    = as.numeric(op$actual),
      predicted = as.numeric(op$predicted),
      stringsAsFactors = FALSE
    )
  }))

  all_dates <- sort(unique(oos_long$date))
  actual_vec <- vapply(all_dates, function(d) {
    vals <- oos_long$actual[oos_long$date == d]
    vals <- vals[!is.na(vals)]
    if (length(vals) == 0) NA_real_ else vals[1]
  }, numeric(1))
  actual_by_date <- data.frame(date = all_dates, actual = actual_vec)

  pred_wide <- tidyr::pivot_wider(
    oos_long[, c("date", "model", "predicted")],
    names_from = "model", values_from = "predicted"
  )
  pred_wide <- as.data.frame(pred_wide)

  oos_wide <- merge(actual_by_date, pred_wide, by = "date", all = TRUE)
  oos_wide <- oos_wide[order(oos_wide$date), ]
  oos_wide <- oos_wide[, c("date", "actual", model_names)]
  rownames(oos_wide) <- NULL

  if (verbose) message("  OOS comparison spans ", format(min(all_dates), "%Y-%m"),
                        " to ", format(max(all_dates), "%Y-%m"),
                        "  (", length(all_dates), " period(s))")

  # ---- Step 2: Latest nowcast (first period after last observed actual) --

  if (verbose) message("Step 2: Identifying latest nowcast period ...")

  actual_dates <- actual_by_date$date[!is.na(actual_by_date$actual)]
  last_actual_date <- if (length(actual_dates) > 0) max(actual_dates) else NA

  future_dates <- if (!is.na(last_actual_date)) all_dates[all_dates > last_actual_date] else all_dates
  target_date <- if (length(future_dates) > 0) min(future_dates) else max(all_dates)

  latest_nowcast <- do.call(rbind, lapply(model_names, function(m) {
    val <- oos_long$predicted[oos_long$model == m & oos_long$date == target_date]
    data.frame(
      model   = m,
      date    = target_date,
      nowcast = if (length(val) > 0) val[1] else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  latest_nowcast$date <- as.Date(latest_nowcast$date, origin = "1970-01-01")

  if (verbose) {
    message("  Last observed actual: ",
            if (is.na(last_actual_date)) "none" else format(last_actual_date, "%Y-%m"))
    message("  Nowcast target period: ", format(target_date, "%Y-%m"))
    message("  Nowcasts:\n",
            paste(utils::capture.output(print(latest_nowcast, row.names = FALSE)), collapse = "\n"))
  }

  # ---- Step 3: Evaluation metrics table --------------------------------

  if (verbose) message("Step 3: Computing OOS evaluation metrics ...")

  eval_metrics <- do.call(rbind, lapply(model_names, function(m) {
    cbind(Model = m, summary_oos_metrics(model_outputs[[m]]$oos_preds))
  }))
  rownames(eval_metrics) <- NULL

  if (verbose) message("  ", paste(utils::capture.output(
    print(eval_metrics, row.names = FALSE)
  ), collapse = "\n  "))

  # ---- Step 4: Charts ----------------------------------------------------

  if (!is.null(out_dir) && requireNamespace("ggplot2", quietly = TRUE)) {

    if (verbose) message("Step 4: Building comparison charts ...")

    plot_long <- rbind(
      data.frame(date = oos_long$date, model = oos_long$model,
                series = "Actual", value = oos_long$actual, stringsAsFactors = FALSE),
      data.frame(date = oos_long$date, model = oos_long$model,
                series = "Nowcast", value = oos_long$predicted, stringsAsFactors = FALSE)
    )
    plot_long$model <- factor(plot_long$model, levels = model_names)

    series_colours    <- c("Actual" = "#0C4550", "Nowcast" = "#800000")
    series_linetypes  <- c("Actual" = "solid",   "Nowcast" = "dashed")

    # ---- 4a) Stacked (faceted) comparison chart, one row per model -----
    p_stacked <- ggplot2::ggplot(plot_long, ggplot2::aes(x = date, y = value,
                                                          colour = series, linetype = series)) +
      ggplot2::geom_line(linewidth = 0.8, na.rm = TRUE) +
      ggplot2::geom_point(data = plot_long[plot_long$series == "Nowcast", ],
                          size = 1.6, na.rm = TRUE) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
      ggplot2::facet_wrap(~model, ncol = 1, scales = "free_y") +
      ggplot2::scale_colour_manual(values = series_colours) +
      ggplot2::scale_linetype_manual(values = series_linetypes) +
      ggplot2::scale_x_date(date_breaks = "6 months", date_labels = "%Y-%m") +
      ggplot2::labs(
        title = paste0(dep_var, " — OOS Actual vs Nowcast Across Models (", model_label, ")"),
        x = NULL, y = paste0("Annual Growth Rate (", dep_var, ")"), colour = NULL, linetype = NULL
      ) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        legend.position = "bottom",
        axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1),
        strip.text      = ggplot2::element_text(face = "bold")
      )

    stacked_path <- file.path(out_dir, paste0("summary_oos_stacked_", safe_label, ".png"))
    ggplot2::ggsave(stacked_path, p_stacked, width = 10,
                    height = max(4, 2.2 * length(model_names)), dpi = 150, limitsize = FALSE)
    if (verbose) message("  Saved: ", basename(stacked_path))

    # ---- 4b) Individual per-model comparison charts ---------------------
    for (m in model_names) {
      m_data <- plot_long[plot_long$model == m, ]
      m_metrics <- eval_metrics[eval_metrics$Model == m, ]
      metrics_label <- if (!is.na(m_metrics$MAE)) {
        sprintf("OOS MAE = %.4f  |  RMSE = %.4f  (n = %d)", m_metrics$MAE, m_metrics$RMSE, m_metrics$N_Obs)
      } else "No evaluated OOS periods yet (actuals not yet available)"

      p_m <- ggplot2::ggplot(m_data, ggplot2::aes(x = date, y = value,
                                                    colour = series, linetype = series)) +
        ggplot2::geom_line(linewidth = 0.9, na.rm = TRUE) +
        ggplot2::geom_point(data = m_data[m_data$series == "Nowcast", ], size = 2.2, na.rm = TRUE) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dotted", colour = "grey60") +
        ggplot2::scale_colour_manual(values = series_colours) +
        ggplot2::scale_linetype_manual(values = series_linetypes) +
        ggplot2::scale_x_date(date_breaks = "3 months", date_labels = "%Y-%m") +
        ggplot2::labs(
          title = paste0(dep_var, " — ", m, " (", model_label, "): OOS Actual vs Nowcast"),
          subtitle = metrics_label,
          x = NULL, y = paste0("Annual Growth Rate (", dep_var, ")"), colour = NULL, linetype = NULL
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          legend.position = "bottom",
          axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1),
          plot.subtitle   = ggplot2::element_text(size = 9, colour = "grey30", family = "mono")
        )

      safe_m <- gsub("[^A-Za-z0-9_]", "_", m)
      m_path <- file.path(out_dir, paste0("summary_oos_", safe_m, "_", safe_label, ".png"))
      ggplot2::ggsave(m_path, p_m, width = 10, height = 5, dpi = 150)
    }
    if (verbose) message("  Saved: ", length(model_names), " individual per-model chart(s)")
  }

  # ---- Step 5: Excel export — OOS comparison ---------------------------

  if (!is.null(out_dir) && requireNamespace("openxlsx", quietly = TRUE)) {

    if (verbose) message("Step 5: Writing Excel outputs ...")

    wb1 <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb1, "OOS_Wide")
    openxlsx::writeDataTable(wb1, "OOS_Wide", oos_wide, tableStyle = "TableStyleMedium9")
    openxlsx::setColWidths(wb1, "OOS_Wide", cols = seq_len(ncol(oos_wide)), widths = "auto")

    openxlsx::addWorksheet(wb1, "OOS_Long")
    openxlsx::writeDataTable(wb1, "OOS_Long", oos_long, tableStyle = "TableStyleMedium9")
    openxlsx::setColWidths(wb1, "OOS_Long", cols = seq_len(ncol(oos_long)), widths = "auto")

    oos_xl_path <- file.path(out_dir, paste0("summary_oos_all_models_", safe_label, ".xlsx"))
    openxlsx::saveWorkbook(wb1, oos_xl_path, overwrite = TRUE)
    if (verbose) message("  Saved: ", basename(oos_xl_path))

    # ---- Step 6: Excel export — latest nowcast --------------------------

    wb2 <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb2, "Latest_Nowcast")
    openxlsx::writeDataTable(wb2, "Latest_Nowcast", latest_nowcast, tableStyle = "TableStyleMedium2")
    openxlsx::setColWidths(wb2, "Latest_Nowcast", cols = seq_len(ncol(latest_nowcast)), widths = "auto")

    info_df <- data.frame(
      Field = c("dep_var", "last_actual_date", "target_nowcast_date", "n_models"),
      Value = c(dep_var,
               if (is.na(last_actual_date)) "NA" else format(last_actual_date, "%Y-%m-%d"),
               format(target_date, "%Y-%m-%d"),
               as.character(length(model_names))),
      stringsAsFactors = FALSE
    )
    openxlsx::addWorksheet(wb2, "Info")
    openxlsx::writeDataTable(wb2, "Info", info_df, tableStyle = "TableStyleMedium2")

    nowcast_xl_path <- file.path(out_dir, paste0("summary_nowcast_latest_", safe_label, ".xlsx"))
    openxlsx::saveWorkbook(wb2, nowcast_xl_path, overwrite = TRUE)
    if (verbose) message("  Saved: ", basename(nowcast_xl_path))

    # ---- Step 7: Excel export — evaluation metrics (best highlighted) ---

    wb3 <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb3, "Evaluation_Metrics")
    openxlsx::writeDataTable(wb3, "Evaluation_Metrics", eval_metrics, tableStyle = "TableStyleMedium9")
    openxlsx::setColWidths(wb3, "Evaluation_Metrics", cols = seq_len(ncol(eval_metrics)), widths = "auto")

    # lower-is-best / higher-is-best / closest-to-zero-is-best per metric
    metric_direction <- list(
      MAE = "min", RMSE = "min", MAPE = "min",
      Bias = "min_abs", Correlation = "max", HitRate = "max",
      DirAccuracy = "max", TheilU2 = "min"
    )
    best_style <- openxlsx::createStyle(fgFill = "#C6EFCE", textDecoration = "bold")

    for (metric in names(metric_direction)) {
      if (!metric %in% names(eval_metrics)) next
      direction <- metric_direction[[metric]]
      col_vals  <- eval_metrics[[metric]]
      cmp_vals  <- if (direction == "min_abs") abs(col_vals) else col_vals
      if (all(is.na(cmp_vals))) next

      best_rows <- if (direction %in% c("min", "min_abs")) {
        which(cmp_vals == min(cmp_vals, na.rm = TRUE))
      } else {
        which(cmp_vals == max(cmp_vals, na.rm = TRUE))
      }
      col_idx <- which(names(eval_metrics) == metric)
      openxlsx::addStyle(wb3, "Evaluation_Metrics", style = best_style,
                         rows = best_rows + 1, cols = col_idx,
                         gridExpand = TRUE, stack = TRUE)
    }

    eval_xl_path <- file.path(out_dir, paste0("summary_evaluation_metrics_", safe_label, ".xlsx"))
    openxlsx::saveWorkbook(wb3, eval_xl_path, overwrite = TRUE)
    if (verbose) message("  Saved: ", basename(eval_xl_path))
  }

  # ---- Step 8: Save RDS -------------------------------------------------

  if (!is.null(out_dir)) {
    rds_path <- file.path(out_dir, paste0("summary_model_data_", safe_label, ".rds"))
    saveRDS(
      list(
        dep_var          = dep_var,
        oos_wide         = oos_wide,
        oos_long         = oos_long,
        latest_nowcast   = latest_nowcast,
        last_actual_date = last_actual_date,
        target_date      = target_date,
        eval_metrics     = eval_metrics
      ),
      file = rds_path
    )
    if (verbose) message("  Saved: ", basename(rds_path))
  }

  if (verbose) message("=== Model summary complete. ===\n")

  invisible(list(
    oos_wide         = oos_wide,
    oos_long         = oos_long,
    latest_nowcast   = latest_nowcast,
    last_actual_date = last_actual_date,
    target_date      = target_date,
    eval_metrics     = eval_metrics
  ))
}

