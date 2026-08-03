# ============================================================
# model_bridge.R
# Bridge model estimation, nowcasting, and output generation.
# ============================================================


#' Estimate a Bridge model and generate nowcasts
#'
#' Runs the full Bridge model pipeline on a transformed monthly dataset:
#' estimates a full-sample OLS regression using the supplied predictor
#' variables, applies ARIMA imputation to fill predictor trailing NAs,
#' generates an expanding-window OOS nowcast, and saves all outputs.
#'
#' Users specify their baseline model by passing a character vector of
#' predictor column names (which must already exist in \code{dta_trans},
#' e.g. as produced by \code{\link{np_transform_data}}) to
#' \code{model_vars}.
#'
#' @param dta_trans Data frame as returned by \code{\link{np_transform_data}}.
#'   Must contain a \code{date} column of class \code{Date} and all columns
#'   named in \code{model_vars}.
#' @param dep_var Character. Dependent variable column name, e.g.
#'   \code{"im_SA"}.
#' @param model_vars Character vector of predictor column names to include in
#'   the Bridge model.  These must already be present in \code{dta_trans}.
#'   Example:
#'   \code{c("im_SA_lag12", "loan_transport_SA", "crisis_dummy")}.
#' @param oos_start Date. Start of the OOS nowcast window.
#' @param oos_end Date. End of the OOS nowcast window (may extend beyond the
#'   last observed actual, producing true nowcasts).
#' @param out_dir Character. Directory where output files are written.
#'   Created if it does not exist.  Defaults to \code{NULL} (no files saved).
#' @param model_label Character. Short label used in chart titles and file
#'   names.  Defaults to \code{"Baseline"}.
#' @param verbose Logical. Print progress messages.  Defaults to \code{TRUE}.
#'
#' @return Invisibly, a named list with elements:
#' \describe{
#'   \item{\code{full_fit}}{The \code{lm} object estimated on all available
#'     observations.}
#'   \item{\code{nowcast_tbl}}{Data frame combining in-sample fitted values
#'     and OOS expanding-window nowcasts.  Columns: \code{date},
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
#' # Add any extra lags / dummies to dta_trans here if needed, then:
#' bridge_out <- np_model_bridge(
#'   dta_trans   = dta_trans,
#'   dep_var     = "im_SA",
#'   model_vars  = c("im_SA_lag12", "loan_transport_SA",
#'                   "crisis_dummy", "recovery_dummy"),
#'   oos_start   = as.Date("2024-01-01"),
#'   oos_end     = as.Date("2026-04-01"),
#'   out_dir     = "Outputs/Bridge",
#'   model_label = "Baseline"
#' )
#' summary(bridge_out$full_fit)
#' }
#'
#' @export
np_model_bridge <- function(dta_trans,
                            dep_var,
                            model_vars,
                            oos_start,
                            oos_end,
                            out_dir     = NULL,
                            model_label = "Baseline",
                            verbose     = TRUE) {

  oos_start <- as.Date(oos_start)
  oos_end   <- as.Date(oos_end)

  if (!dep_var %in% names(dta_trans))
    stop("'", dep_var, "' not found in dta_trans.")

  lag_res    <- auto_create_lags(model_vars, dta_trans, verbose = verbose)
  dta_trans  <- lag_res$dta_trans
  model_vars <- lag_res$vars

  if (!is.null(out_dir) && !dir.exists(out_dir))
    dir.create(out_dir, recursive = TRUE)

  safe_label <- gsub("[^A-Za-z0-9_]", "_", model_label)

  # ---- Step 1: Full-sample regression ---------------------------

  if (verbose) message("Step 1: Full-sample OLS regression ...")

  full_train_df <- dta_trans[, c("date", dep_var, model_vars), drop = FALSE]
  full_train_df <- stats::na.omit(full_train_df)

  if (nrow(full_train_df) <= length(model_vars) + 1)
    stop("Insufficient observations for full-sample regression after listwise deletion.")

  bridge_fml <- stats::as.formula(
    paste(dep_var, "~", paste(model_vars, collapse = " + "))
  )
  full_fit <- stats::lm(bridge_fml, data = full_train_df)

  if (verbose) {
    message("  Sample: ", format(min(full_train_df$date), "%Y-%m"),
            " to ",       format(max(full_train_df$date), "%Y-%m"),
            "  (n = ",    nrow(full_train_df), ")")
    message("  Adj R\u00b2 = ", round(summary(full_fit)$adj.r.squared, 4))
  }

  # ---- Step 2: Save regression table PNG ------------------------

  if (!is.null(out_dir)) {
    sg_html <- file.path(out_dir, paste0("regression_", safe_label, ".html"))
    sg_png  <- file.path(out_dir, paste0("regression_", safe_label, ".png"))

    html_content <- make_reg_html(
      full_fit,
      title     = paste0("Bridge Model \u2014 ", model_label, " (", dep_var, ")"),
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

  # ---- Step 3: ARIMA imputation for predictor trailing NAs ------

  if (verbose) message("Step 2: ARIMA imputation for predictor trailing NAs ...")
  n_imputed <- 0L
  for (col in model_vars) {
    if (any(is.na(dta_trans[[col]]))) {
      dta_trans[[col]] <- arima_impute_col(dta_trans[[col]])
      n_imputed <- n_imputed + 1L
    }
  }
  if (verbose) message("  Imputed ", n_imputed, " column(s).")

  # ---- Step 4: OOS expanding-window nowcast ---------------------

  if (verbose) message("Step 3: Generating OOS nowcasts (expanding window) ...")

  oos_months <- dta_trans$date[dta_trans$date >= oos_start &
                                 dta_trans$date <= oos_end]

  oos_preds <- do.call(rbind, lapply(oos_months, function(t_month) {
    train_df <- dta_trans[dta_trans$date < t_month,
                          c(dep_var, model_vars), drop = FALSE]
    train_df <- stats::na.omit(train_df)

    actual <- dta_trans[dta_trans$date == t_month, dep_var, drop = TRUE]
    actual <- if (length(actual) > 0) actual[1] else NA_real_

    pred_row <- dta_trans[dta_trans$date == t_month, model_vars, drop = FALSE]

    if (nrow(train_df) < length(model_vars) + 2 || any(is.na(pred_row))) {
      return(data.frame(date = t_month, actual = actual, predicted = NA_real_,
                        stringsAsFactors = FALSE))
    }

    fit <- tryCatch(stats::lm(bridge_fml, data = train_df),
                    error = function(e) NULL)
    if (is.null(fit))
      return(data.frame(date = t_month, actual = actual, predicted = NA_real_,
                        stringsAsFactors = FALSE))

    data.frame(
      date      = t_month,
      actual    = actual,
      predicted = as.numeric(stats::predict(fit, newdata = pred_row)),
      stringsAsFactors = FALSE
    )
  }))

  # Restore Date class lost by do.call(rbind) on data frames with Date columns
  oos_preds$date <- as.Date(oos_preds$date, origin = "1970-01-01")

  # ---- Step 5: Combine in-sample + OOS into nowcast_tbl ---------

  insample_df <- data.frame(
    date    = full_train_df$date,
    actual  = full_train_df[[dep_var]],
    nowcast = as.numeric(stats::fitted(full_fit)),
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

  # ---- Step 6: Nowcast chart ------------------------------------

  if (!is.null(out_dir) && requireNamespace("ggplot2", quietly = TRUE)) {
    chart_data <- nowcast_tbl[!is.na(nowcast_tbl$nowcast), ]
    oos_pts    <- oos_preds[!is.na(oos_preds$predicted), ]

    p <- ggplot2::ggplot(chart_data, ggplot2::aes(x = date)) +
      ggplot2::geom_line(ggplot2::aes(y = actual,  colour = "Actual"),
                         linewidth = 0.8, na.rm = TRUE) +
      ggplot2::geom_line(ggplot2::aes(y = nowcast, colour = "Bridge Model",
                                      linetype = type),
                         linewidth = 0.8, na.rm = TRUE) +
      ggplot2::geom_point(data = oos_pts,
                          ggplot2::aes(y = predicted, colour = "Bridge Model"),
                          size = 2) +
      ggplot2::geom_vline(xintercept = as.numeric(oos_start),
                          linetype = "dotdash", colour = "grey50",
                          linewidth = 0.6) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dotted",
                          colour = "grey60") +
      ggplot2::scale_colour_manual(
        values = c("Actual" = "#0C4550", "Bridge Model" = "#E8823C")
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
        title    = paste0(dep_var, " \u2014 Bridge Model (", model_label,
                          "): Actual vs Nowcast"),
        subtitle = paste0("Solid = in-sample fit  |  Dashed + points = OOS expanding-window\n",
                          "Predictors: ",
                          paste(model_vars, collapse = ", ")),
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
                            paste0("bridge_nowcast_chart_", safe_label, ".png"))
    ggplot2::ggsave(chart_path, p, width = 12, height = 6, dpi = 150)
    if (verbose) message("  Saved: ", basename(chart_path))

    # ---- OOS-only chart with evaluation metrics -------------------
    oos_eval <- oos_preds[!is.na(oos_preds$predicted) &
                            !is.na(oos_preds$actual), ]

    # Compute metrics when actuals are available; otherwise label as pending
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

    oos_plot_data       <- oos_preds
    oos_plot_data$date  <- as.Date(oos_plot_data$date, origin = "1970-01-01")

    p_oos <- ggplot2::ggplot(oos_plot_data,
                             ggplot2::aes(x = date)) +
      ggplot2::geom_line(ggplot2::aes(y = actual,    colour = "Actual"),
                         linewidth = 0.9, na.rm = TRUE) +
      ggplot2::geom_line(ggplot2::aes(y = predicted, colour = "Bridge Model"),
                         linewidth = 0.9, linetype = "dashed", na.rm = TRUE) +
      ggplot2::geom_point(ggplot2::aes(y = predicted, colour = "Bridge Model"),
                          size = 2.5, na.rm = TRUE) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dotted",
                          colour = "grey60") +
      ggplot2::scale_colour_manual(
        values = c("Actual" = "#0C4550", "Bridge Model" = "#E8823C")
      ) +
      ggplot2::scale_x_date(date_breaks = "3 months",
                            date_labels = "%Y-%m") +
      ggplot2::labs(
        title    = paste0(dep_var, " \u2014 Bridge Model (",
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
                                paste0("bridge_oos_eval_", safe_label, ".png"))
    ggplot2::ggsave(oos_chart_path, p_oos, width = 10, height = 5, dpi = 150)
    if (verbose) message("  Saved: ", basename(oos_chart_path),
                         if (!is.na(oos_mae))
                           paste0("  (MAE=", round(oos_mae, 4),
                                  ", RMSE=", round(oos_rmse, 4), ")")
                         else "  (no actuals yet)")
  }

  # ---- Step 7: Excel export -------------------------------------

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

    openxlsx::addWorksheet(wb, "Model_Spec")
    spec_df <- data.frame(
      predictor    = model_vars,
      in_model     = TRUE,
      stringsAsFactors = FALSE
    )
    openxlsx::writeDataTable(wb, "Model_Spec", spec_df,
                             tableStyle = "TableStyleMedium2")

    xl_path <- file.path(out_dir,
                         paste0("bridge_nowcast_", safe_label, ".xlsx"))
    openxlsx::saveWorkbook(wb, xl_path, overwrite = TRUE)
    if (verbose) message("  Saved: ", basename(xl_path))
  }

  # ---- Step 8: Save RDS -----------------------------------------

  if (!is.null(out_dir)) {
    rds_path <- file.path(out_dir,
                          paste0("bridge_model_data_", safe_label, ".rds"))
    saveRDS(
      list(
        dta_trans   = dta_trans,
        model_vars  = model_vars,
        dep_var     = dep_var,
        bridge_fml  = bridge_fml,
        full_fit    = full_fit,
        nowcast_tbl = nowcast_tbl,
        oos_preds   = oos_preds,
        oos_start   = oos_start,
        oos_end     = oos_end
      ),
      file = rds_path
    )
    if (verbose) message("  Saved: ", basename(rds_path))
  }

  if (verbose) message("=== Bridge model complete. ===\n")

  invisible(list(
    full_fit    = full_fit,
    nowcast_tbl = nowcast_tbl,
    oos_preds   = oos_preds,
    dta_trans   = dta_trans
  ))
}
