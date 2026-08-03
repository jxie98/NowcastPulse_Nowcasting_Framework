# ============================================================
# model_3prf.R
# Three-Pass Regression Filter (3PRF) model estimation,
# nowcasting, and output generation. Single proxy = dep_var.
# ============================================================


# ---- Internal: 3PRF Pass 1 — indicator loadings on the proxy --------
#' @noRd
tprf_pass1_loadings <- function(indicators_matrix, proxy) {
  tprf_vars <- colnames(indicators_matrix)
  loadings <- vapply(tprf_vars, function(v) {
    df  <- data.frame(y = indicators_matrix[, v], x = as.numeric(proxy))
    fit <- tryCatch(stats::lm(y ~ x, data = df), error = function(e) NULL)
    if (is.null(fit) || !"x" %in% names(stats::coef(fit))) return(NA_real_)
    as.numeric(stats::coef(fit)[["x"]])
  }, numeric(1))
  names(loadings) <- tprf_vars
  loadings
}


# ---- Internal: 3PRF Pass 2 — cross-sectional factor from loadings ---
#' @noRd
tprf_pass2_factor <- function(indicator_row, loadings) {
  ok <- is.finite(as.numeric(loadings)) & is.finite(as.numeric(indicator_row))
  if (sum(ok) < 2) return(NA_real_)
  df  <- data.frame(y = as.numeric(indicator_row)[ok], x = as.numeric(loadings)[ok])
  fit <- tryCatch(stats::lm(y ~ x - 1, data = df), error = function(e) NULL)
  if (is.null(fit)) return(NA_real_)
  as.numeric(stats::coef(fit)[1])
}


#' Estimate a Three-Pass Regression Filter (3PRF) model and generate nowcasts
#'
#' Runs the full 3PRF pipeline on a transformed dataset: extracts a single
#' target-proxy factor from a pool of predictor variables (\code{tprf_vars}),
#' using \code{dep_var} itself as the proxy (Kelly & Pruitt, 2015), estimates
#' a full-sample OLS regression of \code{dep_var} on that factor plus any
#' fixed variables (e.g. AR lags, dummies), applies ARIMA imputation to fill
#' predictor trailing NAs, generates an expanding-window OOS nowcast, and
#' saves all outputs.
#'
#' The factor is built in two passes:
#' \enumerate{
#'   \item \strong{Pass 1} (time-series): each indicator in \code{tprf_vars}
#'     is regressed on the proxy (\code{dep_var}) across periods where both
#'     are observed, giving one loading per indicator.
#'   \item \strong{Pass 2} (cross-sectional): for a given period, that
#'     period's indicator readings are regressed on the Pass-1 loadings
#'     (across indicators, no intercept), giving the factor value for that
#'     period.
#' }
#' Pass 2 only needs that period's indicator readings and the already-fitted
#' loadings — not \code{dep_var} — so a factor can be produced for a period
#' whose \code{dep_var} is still unknown (the whole point of nowcasting).
#' At each OOS step, Pass 1 is re-estimated using only data strictly before
#' the target period (so it never sees the value it is trying to predict),
#' and Pass 2 is then applied to the target period's own indicator readings
#' to get that period's factor; the regression on top is likewise
#' re-estimated on data strictly before the target period. This differs from
#' a literal per-vintage replication of the reference methodology (which
#' would drop the target period from the factor step entirely, since it
#' listwise-deletes on \code{dep_var} before extracting the factor) — here
#' the factor is deliberately available for the very period being nowcast.
#'
#' @param dta_trans Data frame as returned by \code{\link{np_transform_data}}.
#'   Must contain a \code{date} column of class \code{Date} and all columns
#'   named in \code{tprf_vars} and \code{fixed_vars}.
#' @param dep_var Character. Dependent variable column name, e.g.
#'   \code{"im_SA"}. Also used as the 3PRF proxy.
#' @param tprf_vars Character vector of predictor column names from which the
#'   factor is extracted. These must already be present in \code{dta_trans}.
#'   At least two are required (Pass 2 is a cross-sectional regression across
#'   indicators).
#' @param fixed_vars Character vector of predictor column names to include
#'   directly in the regression alongside the factor (e.g. AR lags,
#'   structural dummies). Names of the form \code{X_lagN} whose base
#'   \code{X} exists in \code{dta_trans} are auto-created. Defaults to
#'   \code{NULL}.
#' @param oos_start Date. Start of the OOS nowcast window.
#' @param oos_end Date. End of the OOS nowcast window (may extend beyond the
#'   last observed actual, producing true nowcasts).
#' @param out_dir Character. Directory where output files are written.
#'   Created if it does not exist. Defaults to \code{NULL} (no files saved).
#' @param model_label Character. Short label used in chart titles and file
#'   names. Defaults to \code{"3PRF"}.
#' @param verbose Logical. Print progress messages. Defaults to \code{TRUE}.
#'
#' @return Invisibly, a named list with elements:
#' \describe{
#'   \item{\code{full_fit}}{The \code{lm} object estimated on all available
#'     observations, using the full-sample factor.}
#'   \item{\code{loadings}}{Named numeric vector of full-sample Pass-1
#'     loadings, one per \code{tprf_vars}.}
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
#' tprf_out <- np_model_3prf(
#'   dta_trans   = dta_trans,
#'   dep_var     = "im_SA",
#'   tprf_vars   = c("cpi_SA", "iip_SA", "petrol_SA", "money_supply_SA"),
#'   fixed_vars  = c("im_SA_lag1", "im_SA_lag12", "crisis_dummy"),
#'   oos_start   = as.Date("2024-01-01"),
#'   oos_end     = as.Date("2026-04-01"),
#'   out_dir     = "Outputs/3PRF",
#'   model_label = "Baseline"
#' )
#' summary(tprf_out$full_fit)
#' }
#'
#' @export
np_model_3prf <- function(dta_trans,
                          dep_var,
                          tprf_vars,
                          fixed_vars  = NULL,
                          oos_start,
                          oos_end,
                          out_dir     = NULL,
                          model_label = "3PRF",
                          verbose     = TRUE) {

  oos_start <- as.Date(oos_start)
  oos_end   <- as.Date(oos_end)

  if (!dep_var %in% names(dta_trans))
    stop("'", dep_var, "' not found in dta_trans.")

  if (length(tprf_vars) < 2)
    stop("tprf_vars must contain at least 2 variables for the Pass 2 ",
         "cross-sectional regression.")

  missing_tprf <- setdiff(tprf_vars, names(dta_trans))
  if (length(missing_tprf) > 0)
    stop("tprf_vars not found in dta_trans: ", paste(missing_tprf, collapse = ", "))

  lag_res     <- auto_create_lags(fixed_vars, dta_trans, verbose = verbose)
  dta_trans   <- lag_res$dta_trans
  fixed_vars  <- lag_res$vars

  if (!is.null(out_dir) && !dir.exists(out_dir))
    dir.create(out_dir, recursive = TRUE)

  safe_label <- gsub("[^A-Za-z0-9_]", "_", model_label)

  # ---- Step 1: Full-sample 3PRF factor extraction ------------------

  if (verbose) message("Step 1: Full-sample 3PRF factor extraction on ",
                        length(tprf_vars), " indicator(s), proxy = ", dep_var, " ...")

  tprf_input <- dta_trans[, c("date", tprf_vars, dep_var), drop = FALSE]
  tprf_input <- stats::na.omit(tprf_input)

  if (nrow(tprf_input) < length(tprf_vars) + 2)
    stop("Insufficient observations for 3PRF after listwise deletion.")

  indicators_matrix <- as.matrix(tprf_input[, tprf_vars, drop = FALSE])
  proxy <- tprf_input[[dep_var]]

  loadings <- tprf_pass1_loadings(indicators_matrix, proxy)
  if (all(!is.finite(loadings)))
    stop("3PRF Pass 1 failed to estimate any indicator loadings on the proxy.")

  factors   <- apply(indicators_matrix, 1, tprf_pass2_factor, loadings = loadings)
  factor_df <- data.frame(date = tprf_input$date, factor_3prf = as.numeric(factors),
                          stringsAsFactors = FALSE)

  if (verbose) {
    message("  Pass 1 loadings (proxy = ", dep_var, "):")
    for (v in tprf_vars)
      message("    ", formatC(v, width = 30, flag = "-"), " ", round(loadings[[v]], 4))
  }

  # ---- Step 2: Full-sample regression -----------------------------

  if (verbose) message("Step 2: Full-sample OLS regression on factor ...")

  if (length(fixed_vars) == 0) {
    full_train_df <- dta_trans[, c("date", dep_var), drop = FALSE]
  } else {
    full_train_df <- dta_trans[, c("date", dep_var, fixed_vars), drop = FALSE]
  }
  full_train_df <- merge(full_train_df, factor_df, by = "date")
  full_train_df <- stats::na.omit(full_train_df[, c("date", dep_var, fixed_vars, "factor_3prf"), drop = FALSE])

  if (nrow(full_train_df) <= length(fixed_vars) + 2)
    stop("Insufficient observations for full-sample regression after listwise deletion.")

  tprf_fml <- stats::as.formula(
    paste(dep_var, "~", paste(c(fixed_vars, "factor_3prf"), collapse = " + "))
  )
  full_fit <- stats::lm(tprf_fml, data = full_train_df)

  if (verbose) {
    message("  Sample: ", format(min(full_train_df$date), "%Y-%m"),
            " to ",       format(max(full_train_df$date), "%Y-%m"),
            "  (n = ",    nrow(full_train_df), ")")
    message("  Adj R² = ", round(summary(full_fit)$adj.r.squared, 4))
  }

  # ---- Step 3: Save regression table PNG ---------------------------

  if (!is.null(out_dir)) {
    sg_html <- file.path(out_dir, paste0("regression_", safe_label, ".html"))
    sg_png  <- file.path(out_dir, paste0("regression_", safe_label, ".png"))

    html_content <- make_reg_html(
      full_fit,
      title     = paste0("3PRF Model — ", model_label, " (", dep_var, ")"),
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

  # ---- Step 4: ARIMA imputation for predictor trailing NAs ---------

  if (verbose) message("Step 3: ARIMA imputation for predictor trailing NAs ...")
  predictor_cols <- unique(c(tprf_vars, fixed_vars))
  n_imputed <- 0L
  for (col in predictor_cols) {
    if (any(is.na(dta_trans[[col]]))) {
      dta_trans[[col]] <- arima_impute_col(dta_trans[[col]])
      n_imputed <- n_imputed + 1L
    }
  }
  if (verbose) message("  Imputed ", n_imputed, " column(s).")

  # ---- Step 5: OOS expanding-window nowcast -------------------------
  # Pass 1 loadings are re-estimated each period from data strictly before
  # it (dep_var known there); Pass 2 is then applied to that period's own
  # indicator readings using those loadings, so a factor is available for
  # the very period being nowcast without needing its (unknown) dep_var.
  # The regression on top is re-estimated on data strictly before the
  # period too, so nothing here uses future information.

  if (verbose) message("Step 4: Generating OOS nowcasts (expanding window) ...")

  oos_months <- dta_trans$date[dta_trans$date >= oos_start &
                                 dta_trans$date <= oos_end]

  oos_preds <- do.call(rbind, lapply(oos_months, function(t_month) {

    actual <- dta_trans[dta_trans$date == t_month, dep_var, drop = TRUE]
    actual <- if (length(actual) > 0) actual[1] else NA_real_

    na_row <- function() data.frame(date = t_month, actual = actual,
                                    predicted = NA_real_, stringsAsFactors = FALSE)

    # Pass 1: loadings from history only (proxy = dep_var, known there)
    train_tprf <- dta_trans[dta_trans$date < t_month, c("date", tprf_vars, dep_var), drop = FALSE]
    train_tprf <- stats::na.omit(train_tprf)
    if (nrow(train_tprf) < length(tprf_vars) + 2) return(na_row())

    train_matrix <- as.matrix(train_tprf[, tprf_vars, drop = FALSE])
    win_loadings <- tprf_pass1_loadings(train_matrix, train_tprf[[dep_var]])
    if (all(!is.finite(win_loadings))) return(na_row())

    # Pass 2: factor for every training date plus the target date itself
    win_dates_df <- dta_trans[dta_trans$date <= t_month, c("date", tprf_vars), drop = FALSE]
    win_dates_df <- stats::na.omit(win_dates_df)
    if (!(t_month %in% win_dates_df$date)) return(na_row())

    win_matrix   <- as.matrix(win_dates_df[, tprf_vars, drop = FALSE])
    win_factors  <- apply(win_matrix, 1, tprf_pass2_factor, loadings = win_loadings)
    win_factor_df <- data.frame(date = win_dates_df$date, factor_3prf = as.numeric(win_factors),
                                stringsAsFactors = FALSE)

    if (length(fixed_vars) == 0) {
      train_df <- dta_trans[dta_trans$date < t_month, c("date", dep_var), drop = FALSE]
    } else {
      train_df <- dta_trans[dta_trans$date < t_month, c("date", dep_var, fixed_vars), drop = FALSE]
    }
    train_df <- merge(train_df, win_factor_df, by = "date")
    train_df <- stats::na.omit(train_df[, c(dep_var, fixed_vars, "factor_3prf"), drop = FALSE])

    min_obs <- max(5, length(fixed_vars) + 3)
    if (nrow(train_df) < min_obs) return(na_row())

    if (length(fixed_vars) == 0) {
      fixed_row <- dta_trans[dta_trans$date == t_month, "date", drop = FALSE]
    } else {
      fixed_row <- dta_trans[dta_trans$date == t_month, c("date", fixed_vars), drop = FALSE]
    }
    factor_row <- win_factor_df[win_factor_df$date == t_month, c("date", "factor_3prf"), drop = FALSE]
    pred_row <- merge(fixed_row, factor_row, by = "date")
    if (nrow(pred_row) == 0) return(na_row())
    pred_row <- pred_row[, c(fixed_vars, "factor_3prf"), drop = FALSE]
    if (any(is.na(pred_row))) return(na_row())

    win_fml <- stats::as.formula(
      paste(dep_var, "~", paste(c(fixed_vars, "factor_3prf"), collapse = " + "))
    )
    fit <- tryCatch(stats::lm(win_fml, data = train_df), error = function(e) NULL)
    if (is.null(fit)) return(na_row())

    data.frame(
      date      = t_month,
      actual    = actual,
      predicted = as.numeric(stats::predict(fit, newdata = pred_row)),
      stringsAsFactors = FALSE
    )
  }))

  # Restore Date class lost by do.call(rbind) on data frames with Date columns
  oos_preds$date <- as.Date(oos_preds$date, origin = "1970-01-01")

  # ---- Step 6: Combine in-sample + OOS into nowcast_tbl -------------

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

  # ---- Step 7: Nowcast chart -----------------------------------------

  if (!is.null(out_dir) && requireNamespace("ggplot2", quietly = TRUE)) {
    chart_data <- nowcast_tbl[!is.na(nowcast_tbl$nowcast), ]
    oos_pts    <- oos_preds[!is.na(oos_preds$predicted), ]

    p <- ggplot2::ggplot(chart_data, ggplot2::aes(x = date)) +
      ggplot2::geom_line(ggplot2::aes(y = actual,  colour = "Actual"),
                         linewidth = 0.8, na.rm = TRUE) +
      ggplot2::geom_line(ggplot2::aes(y = nowcast, colour = "3PRF Model",
                                      linetype = type),
                         linewidth = 0.8, na.rm = TRUE) +
      ggplot2::geom_point(data = oos_pts,
                          ggplot2::aes(y = predicted, colour = "3PRF Model"),
                          size = 2) +
      ggplot2::geom_vline(xintercept = as.numeric(oos_start),
                          linetype = "dotdash", colour = "grey50",
                          linewidth = 0.6) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dotted",
                          colour = "grey60") +
      ggplot2::scale_colour_manual(
        values = c("Actual" = "#0C4550", "3PRF Model" = "#2A9D8F")
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
        title    = paste0(dep_var, " — 3PRF Model (", model_label,
                          "): Actual vs Nowcast"),
        subtitle = paste0("Solid = in-sample fit  |  Dashed + points = OOS expanding-window\n",
                          "Factor from: ", paste(tprf_vars, collapse = ", "),
                          "  |  Fixed vars: ",
                          if (length(fixed_vars) > 0) paste(fixed_vars, collapse = ", ") else "none"),
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
                            paste0("tprf_nowcast_chart_", safe_label, ".png"))
    ggplot2::ggsave(chart_path, p, width = 12, height = 6, dpi = 150)
    if (verbose) message("  Saved: ", basename(chart_path))

    # ---- OOS-only chart with evaluation metrics -------------------
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
      ggplot2::geom_line(ggplot2::aes(y = predicted, colour = "3PRF Model"),
                         linewidth = 0.9, linetype = "dashed", na.rm = TRUE) +
      ggplot2::geom_point(ggplot2::aes(y = predicted, colour = "3PRF Model"),
                          size = 2.5, na.rm = TRUE) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dotted",
                          colour = "grey60") +
      ggplot2::scale_colour_manual(
        values = c("Actual" = "#0C4550", "3PRF Model" = "#2A9D8F")
      ) +
      ggplot2::scale_x_date(date_breaks = "3 months",
                            date_labels = "%Y-%m") +
      ggplot2::labs(
        title    = paste0(dep_var, " — 3PRF Model (",
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
                                paste0("tprf_oos_eval_", safe_label, ".png"))
    ggplot2::ggsave(oos_chart_path, p_oos, width = 10, height = 5, dpi = 150)
    if (verbose) message("  Saved: ", basename(oos_chart_path),
                         if (!is.na(oos_mae))
                           paste0("  (MAE=", round(oos_mae, 4),
                                  ", RMSE=", round(oos_rmse, 4), ")")
                         else "  (no actuals yet)")
  }

  # ---- Step 8: Excel export -------------------------------------------

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
      predictor = c(fixed_vars, "factor_3prf"),
      type      = c(rep("fixed", length(fixed_vars)), "factor"),
      stringsAsFactors = FALSE
    )
    openxlsx::writeDataTable(wb, "Model_Spec", spec_df,
                             tableStyle = "TableStyleMedium2")

    openxlsx::addWorksheet(wb, "TPRF_Loadings")
    loadings_df <- data.frame(
      variable = names(loadings),
      loading  = as.numeric(loadings),
      stringsAsFactors = FALSE
    )
    openxlsx::writeDataTable(wb, "TPRF_Loadings", loadings_df,
                             tableStyle = "TableStyleMedium2")

    xl_path <- file.path(out_dir,
                         paste0("tprf_nowcast_", safe_label, ".xlsx"))
    openxlsx::saveWorkbook(wb, xl_path, overwrite = TRUE)
    if (verbose) message("  Saved: ", basename(xl_path))
  }

  # ---- Step 9: Save RDS ------------------------------------------------

  if (!is.null(out_dir)) {
    rds_path <- file.path(out_dir,
                          paste0("tprf_model_data_", safe_label, ".rds"))
    saveRDS(
      list(
        dta_trans   = dta_trans,
        tprf_vars   = tprf_vars,
        fixed_vars  = fixed_vars,
        dep_var     = dep_var,
        loadings    = loadings,
        tprf_fml    = tprf_fml,
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

  if (verbose) message("=== 3PRF model complete. ===\n")

  invisible(list(
    full_fit    = full_fit,
    loadings    = loadings,
    nowcast_tbl = nowcast_tbl,
    oos_preds   = oos_preds,
    dta_trans   = dta_trans
  ))
}
