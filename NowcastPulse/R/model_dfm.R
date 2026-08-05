# ============================================================
# model_dfm.R
# Dynamic Factor Model estimation, nowcasting, and output
# generation. Single common factor extracted via Kalman filter
# (KFAS), initialised from PCA.
# ============================================================


# ---- Internal: extract a single common factor via Kalman filter ----
#' @noRd
extract_dfm_factor <- function(std_matrix) {
  n_indicators <- ncol(std_matrix)

  pca_init      <- stats::prcomp(std_matrix, center = FALSE, scale. = FALSE)
  loadings_init <- pca_init$rotation[, 1]

  # State equation:       f_t = phi * f_{t-1} + eta_t      (AR(1) factor)
  # Observation equation: Y_t = Lambda * f_t + epsilon_t
  Zt <- matrix(loadings_init, nrow = n_indicators, ncol = 1)
  Ht <- diag(0.1, n_indicators)
  Tt <- matrix(0.9, nrow = 1, ncol = 1)
  Rt <- matrix(1,   nrow = 1, ncol = 1)
  Qt <- matrix(1,   nrow = 1, ncol = 1)
  a1 <- matrix(0,   nrow = 1, ncol = 1)
  P1 <- matrix(1,   nrow = 1, ncol = 1)

  # KFAS::SSModel's formula parser only recognises the bare "SSMcustom" term
  # name (a namespace-qualified call is not detected as a model term), and
  # this KFAS version requires an explicit `index` rather than defaulting it
  # to all p series.
  SSMcustom <- KFAS::SSMcustom
  ss_model <- KFAS::SSModel(
    std_matrix ~ -1 +
      SSMcustom(Z = Zt, T = Tt, R = Rt, Q = Qt, a1 = a1, P1 = P1,
               index = seq_len(n_indicators)),
    H = Ht
  )

  fit <- tryCatch(
    KFAS::fitSSM(ss_model, inits = c(0.5, rep(0.1, n_indicators), 0.5),
                method = "BFGS", control = list(maxit = 500)),
    error = function(e) list(model = ss_model, optim.out = list(convergence = 1L))
  )

  smoothed <- tryCatch(
    KFAS::KFS(fit$model, smoothing = c("state", "mean")),
    error = function(e) NULL
  )
  if (is.null(smoothed)) return(NULL)

  list(
    factor      = as.numeric(smoothed$alphahat[, 1]),
    loadings    = stats::setNames(fit$model$Z[, 1, 1], colnames(std_matrix)),
    ar_coef     = fit$model$T[1, 1, 1],
    convergence = fit$optim.out$convergence
  )
}


#' Estimate a Dynamic Factor Model and generate nowcasts
#'
#' Runs the full DFM pipeline on a transformed dataset: extracts a single
#' common factor from a pool of predictor variables (\code{dfm_vars}) via
#' Kalman filtering of a state-space model (state equation: AR(1) factor;
#' observation equation: factor loadings), estimates a full-sample OLS
#' regression of \code{dep_var} on that factor plus any fixed variables
#' (e.g. AR lags, dummies), applies ARIMA imputation to fill predictor
#' trailing NAs, generates an expanding-window OOS nowcast, and saves all
#' outputs.
#'
#' The factor's loadings and AR(1) coefficient are estimated by maximum
#' likelihood (\code{KFAS::fitSSM}, BFGS), initialised from the first
#' principal component of the standardised predictor matrix, and the factor
#' itself is the Kalman-smoothed state. At each OOS step the factor is
#' re-estimated from scratch using only data available up to and including
#' that period (an expanding window, standardised on that window's own mean
#' and SD), and the regression is re-estimated using only data strictly
#' before that period — mirroring how the model would actually be built in
#' real time and avoiding look-ahead bias. Because each step involves a
#' numerical MLE fit, this is markedly slower than \code{\link{np_model_pca}}
#' or \code{\link{np_model_bridge}} over long OOS windows.
#'
#' @param dta_trans Data frame as returned by \code{\link{np_transform_data}}.
#'   Must contain a \code{date} column of class \code{Date} and all columns
#'   named in \code{dfm_vars} and \code{fixed_vars}.
#' @param dep_var Character. Dependent variable column name, e.g.
#'   \code{"im_SA"}.
#' @param dfm_vars Character vector of predictor column names from which the
#'   common factor is extracted. These must already be present in
#'   \code{dta_trans}. At least two are required.
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
#'   names. Defaults to \code{"DFM"}.
#' @param verbose Logical. Print progress messages. Defaults to \code{TRUE}.
#'
#' @return Invisibly, a named list with elements:
#' \describe{
#'   \item{\code{full_fit}}{The \code{lm} object estimated on all available
#'     observations, using the full-sample factor.}
#'   \item{\code{dfm_fit}}{List with the full-sample factor extraction
#'     result: \code{factor}, \code{loadings}, \code{ar_coef},
#'     \code{convergence} (see \code{KFAS::fitSSM}; \code{0} = converged).}
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
#' dfm_out <- np_model_dfm(
#'   dta_trans   = dta_trans,
#'   dep_var     = "im_SA",
#'   dfm_vars    = c("cpi_SA", "iip_SA", "petrol_SA", "money_supply_SA"),
#'   fixed_vars  = c("im_SA_lag1", "im_SA_lag12", "crisis_dummy"),
#'   oos_start   = as.Date("2024-01-01"),
#'   oos_end     = as.Date("2026-04-01"),
#'   out_dir     = "Outputs/DFM",
#'   model_label = "Baseline"
#' )
#' summary(dfm_out$full_fit)
#' }
#'
#' @export
np_model_dfm <- function(dta_trans,
                         dep_var,
                         dfm_vars,
                         fixed_vars  = NULL,
                         oos_start,
                         oos_end,
                         out_dir     = NULL,
                         model_label = "DFM",
                         verbose     = TRUE) {

  oos_start <- as.Date(oos_start)
  oos_end   <- as.Date(oos_end)

  if (!dep_var %in% names(dta_trans))
    stop("'", dep_var, "' not found in dta_trans.")

  if (length(dfm_vars) < 2)
    stop("dfm_vars must contain at least 2 variables to extract a common factor.")

  lag_res_dfm <- auto_create_lags(dfm_vars, dta_trans, verbose = verbose)
  dta_trans   <- lag_res_dfm$dta_trans
  missing_dfm <- setdiff(dfm_vars, lag_res_dfm$vars)
  if (length(missing_dfm) > 0)
    stop("dfm_vars not found in dta_trans: ", paste(missing_dfm, collapse = ", "))
  dfm_vars    <- lag_res_dfm$vars

  lag_res     <- auto_create_lags(fixed_vars, dta_trans, verbose = verbose)
  dta_trans   <- lag_res$dta_trans
  fixed_vars  <- lag_res$vars

  if (!is.null(out_dir) && !dir.exists(out_dir))
    dir.create(out_dir, recursive = TRUE)

  safe_label <- gsub("[^A-Za-z0-9_]", "_", model_label)

  # ---- Step 1: Full-sample factor extraction via Kalman filter ---

  if (verbose) message("Step 1: Full-sample DFM factor extraction (Kalman filter) on ",
                        length(dfm_vars), " indicator(s) ...")

  dfm_input <- dta_trans[, c("date", dfm_vars), drop = FALSE]
  dfm_input <- stats::na.omit(dfm_input)

  if (nrow(dfm_input) < length(dfm_vars) + 2)
    stop("Insufficient observations for DFM after listwise deletion.")

  raw_matrix <- as.matrix(dfm_input[, dfm_vars, drop = FALSE])
  std_matrix <- scale(raw_matrix)

  dfm_fit <- extract_dfm_factor(std_matrix)
  if (is.null(dfm_fit))
    stop("Kalman filter factor extraction failed on the full sample.")

  factor_df <- data.frame(date = dfm_input$date, factor_dfm = dfm_fit$factor,
                          stringsAsFactors = FALSE)

  if (verbose) {
    message("  Convergence status: ", dfm_fit$convergence, "  (0 = success)")
    message("  AR(1) coefficient: ", round(dfm_fit$ar_coef, 4))
  }

  # ---- Step 2: Full-sample regression -----------------------------

  if (verbose) message("Step 2: Full-sample OLS regression on factor ...")

  if (length(fixed_vars) == 0) {
    full_train_df <- dta_trans[, c("date", dep_var), drop = FALSE]
  } else {
    full_train_df <- dta_trans[, c("date", dep_var, fixed_vars), drop = FALSE]
  }
  full_train_df <- merge(full_train_df, factor_df, by = "date")
  full_train_df <- stats::na.omit(full_train_df[, c("date", dep_var, fixed_vars, "factor_dfm"), drop = FALSE])

  if (nrow(full_train_df) <= length(fixed_vars) + 2)
    stop("Insufficient observations for full-sample regression after listwise deletion.")

  dfm_fml <- stats::as.formula(
    paste(dep_var, "~", paste(c(fixed_vars, "factor_dfm"), collapse = " + "))
  )
  full_fit <- stats::lm(dfm_fml, data = full_train_df)

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
      title     = paste0("DFM Model — ", model_label, " (", dep_var, ")"),
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
  predictor_cols <- unique(c(dfm_vars, fixed_vars))
  n_imputed <- 0L
  for (col in predictor_cols) {
    if (any(is.na(dta_trans[[col]]))) {
      dta_trans[[col]] <- arima_impute_col(dta_trans[[col]])
      n_imputed <- n_imputed + 1L
    }
  }
  if (verbose) message("  Imputed ", n_imputed, " column(s).")

  # ---- Step 5: OOS expanding-window nowcast -------------------------
  # The factor is re-extracted each period from data available up to and
  # including that period (standardised on that window's own mean/SD); the
  # regression is re-estimated on data strictly before that period, so
  # neither step uses future information.

  if (verbose) message("Step 4: Generating OOS nowcasts (expanding window) ...")

  oos_months <- dta_trans$date[dta_trans$date >= oos_start &
                                 dta_trans$date <= oos_end]

  oos_preds <- do.call(rbind, lapply(oos_months, function(t_month) {

    actual <- dta_trans[dta_trans$date == t_month, dep_var, drop = TRUE]
    actual <- if (length(actual) > 0) actual[1] else NA_real_

    na_row <- function() data.frame(date = t_month, actual = actual,
                                    predicted = NA_real_, stringsAsFactors = FALSE)

    window_df <- dta_trans[dta_trans$date <= t_month, c("date", dfm_vars), drop = FALSE]
    window_df <- stats::na.omit(window_df)

    if (nrow(window_df) < length(dfm_vars) + 2 || !(t_month %in% window_df$date))
      return(na_row())

    win_std_matrix <- scale(as.matrix(window_df[, dfm_vars, drop = FALSE]))
    win_dfm_fit <- tryCatch(extract_dfm_factor(win_std_matrix), error = function(e) NULL)
    if (is.null(win_dfm_fit)) return(na_row())

    win_factor_df <- data.frame(date = window_df$date, factor_dfm = win_dfm_fit$factor,
                                stringsAsFactors = FALSE)

    if (length(fixed_vars) == 0) {
      train_df <- dta_trans[dta_trans$date < t_month, c("date", dep_var), drop = FALSE]
    } else {
      train_df <- dta_trans[dta_trans$date < t_month, c("date", dep_var, fixed_vars), drop = FALSE]
    }
    train_df <- merge(train_df, win_factor_df, by = "date")
    train_df <- stats::na.omit(train_df[, c(dep_var, fixed_vars, "factor_dfm"), drop = FALSE])

    min_obs <- max(5, length(fixed_vars) + 3)
    if (nrow(train_df) < min_obs) return(na_row())

    if (length(fixed_vars) == 0) {
      fixed_row <- dta_trans[dta_trans$date == t_month, "date", drop = FALSE]
    } else {
      fixed_row <- dta_trans[dta_trans$date == t_month, c("date", fixed_vars), drop = FALSE]
    }
    factor_row <- win_factor_df[win_factor_df$date == t_month, c("date", "factor_dfm"), drop = FALSE]
    pred_row <- merge(fixed_row, factor_row, by = "date")
    if (nrow(pred_row) == 0) return(na_row())
    pred_row <- pred_row[, c(fixed_vars, "factor_dfm"), drop = FALSE]
    if (any(is.na(pred_row))) return(na_row())

    win_fml <- stats::as.formula(
      paste(dep_var, "~", paste(c(fixed_vars, "factor_dfm"), collapse = " + "))
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
      ggplot2::geom_line(ggplot2::aes(y = nowcast, colour = "DFM Model",
                                      linetype = type),
                         linewidth = 0.8, na.rm = TRUE) +
      ggplot2::geom_point(data = oos_pts,
                          ggplot2::aes(y = predicted, colour = "DFM Model"),
                          size = 2) +
      ggplot2::geom_vline(xintercept = as.numeric(oos_start),
                          linetype = "dotdash", colour = "grey50",
                          linewidth = 0.6) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dotted",
                          colour = "grey60") +
      ggplot2::scale_colour_manual(
        values = c("Actual" = "#0C4550", "DFM Model" = "#B23A48")
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
        title    = paste0(dep_var, " — DFM Model (", model_label,
                          "): Actual vs Nowcast"),
        subtitle = paste0("Solid = in-sample fit  |  Dashed + points = OOS expanding-window\n",
                          "Factor from: ", paste(dfm_vars, collapse = ", "),
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
                            paste0("dfm_nowcast_chart_", safe_label, ".png"))
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
      ggplot2::geom_line(ggplot2::aes(y = predicted, colour = "DFM Model"),
                         linewidth = 0.9, linetype = "dashed", na.rm = TRUE) +
      ggplot2::geom_point(ggplot2::aes(y = predicted, colour = "DFM Model"),
                          size = 2.5, na.rm = TRUE) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dotted",
                          colour = "grey60") +
      ggplot2::scale_colour_manual(
        values = c("Actual" = "#0C4550", "DFM Model" = "#B23A48")
      ) +
      ggplot2::scale_x_date(date_breaks = "3 months",
                            date_labels = "%Y-%m") +
      ggplot2::labs(
        title    = paste0(dep_var, " — DFM Model (",
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
                                paste0("dfm_oos_eval_", safe_label, ".png"))
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
      predictor = c(fixed_vars, "factor_dfm"),
      type      = c(rep("fixed", length(fixed_vars)), "factor"),
      stringsAsFactors = FALSE
    )
    openxlsx::writeDataTable(wb, "Model_Spec", spec_df,
                             tableStyle = "TableStyleMedium2")

    openxlsx::addWorksheet(wb, "DFM_Loadings")
    loadings_df <- data.frame(
      variable = names(dfm_fit$loadings),
      loading  = as.numeric(dfm_fit$loadings),
      stringsAsFactors = FALSE
    )
    openxlsx::writeDataTable(wb, "DFM_Loadings", loadings_df,
                             tableStyle = "TableStyleMedium2")

    openxlsx::addWorksheet(wb, "DFM_Diagnostics")
    diag_df <- data.frame(
      ar_coefficient = dfm_fit$ar_coef,
      convergence    = dfm_fit$convergence,
      n_indicators   = length(dfm_vars),
      stringsAsFactors = FALSE
    )
    openxlsx::writeDataTable(wb, "DFM_Diagnostics", diag_df,
                             tableStyle = "TableStyleMedium2")

    xl_path <- file.path(out_dir,
                         paste0("dfm_nowcast_", safe_label, ".xlsx"))
    openxlsx::saveWorkbook(wb, xl_path, overwrite = TRUE)
    if (verbose) message("  Saved: ", basename(xl_path))
  }

  # ---- Step 9: Save RDS ------------------------------------------------

  if (!is.null(out_dir)) {
    rds_path <- file.path(out_dir,
                          paste0("dfm_model_data_", safe_label, ".rds"))
    saveRDS(
      list(
        dta_trans   = dta_trans,
        dfm_vars    = dfm_vars,
        fixed_vars  = fixed_vars,
        dep_var     = dep_var,
        dfm_fit     = dfm_fit,
        dfm_fml     = dfm_fml,
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

  if (verbose) message("=== DFM model complete. ===\n")

  invisible(list(
    full_fit    = full_fit,
    dfm_fit     = dfm_fit,
    nowcast_tbl = nowcast_tbl,
    oos_preds   = oos_preds,
    dta_trans   = dta_trans
  ))
}
