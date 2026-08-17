# ============================================================
# variable_selection.R
# Variable transformation, correlation ranking, and
# forward-backward stepwise OOS MAE variable selection
# for the monthly Bridge model.
# ============================================================


# ---- Internal: build regression table HTML without stargazer --
#' @noRd
make_reg_html <- function(model, title = "Bridge Model", dep_label = NULL) {
  sm       <- summary(model)
  cf       <- as.data.frame(coef(sm))
  names(cf) <- c("Estimate", "Std_Error", "t_value", "p_value")

  stars <- function(p) {
    ifelse(p < 0.01, "***", ifelse(p < 0.05, "**", ifelse(p < 0.1, "*", "")))
  }

  if (is.null(dep_label)) dep_label <- as.character(formula(model)[[2]])

  rows <- mapply(function(nm, est, se, pv) {
    st <- stars(pv)
    sprintf(
      "<tr><td style='text-align:left'>%s</td>
       <td style='text-align:center'>%.3f%s<br><span style='font-size:0.9em'>(%.3f)</span></td></tr>",
      nm, est, st, se
    )
  }, rownames(cf), cf$Estimate, cf$Std_Error, cf$p_value)

  n_obs  <- nrow(model$model)
  r2     <- sm$r.squared
  adj_r2 <- sm$adj.r.squared
  fstat  <- sm$fstatistic
  f_str  <- if (!is.null(fstat))
    sprintf("%.3f<sup>***</sup> (df = %.0f; %.0f)", fstat[1], fstat[2], fstat[3])
  else ""
  rse    <- sprintf("%.3f (df = %.0f)", sm$sigma, sm$df[2])

  sprintf('<!DOCTYPE html>
<html><head><meta charset="utf-8">
<style>
  body { font-family: "Times New Roman", serif; font-size: 13px; margin: 20px; }
  table { border-collapse: collapse; width: 420px; }
  th, td { padding: 4px 10px; }
  .title { font-weight: bold; text-align: center; font-size: 14px; padding-bottom: 4px; }
  .depvar { text-align: center; font-style: italic; border-bottom: 1px solid #000; padding-bottom: 2px; }
  .colhead { text-align: center; border-bottom: 1px solid #000; padding: 2px 10px; }
  .divider { border-top: 1px solid #000; }
  .note { font-size: 0.85em; }
  sup { font-size: 0.75em; color: #b8860b; }
</style></head><body>
<table>
  <tr><td colspan="2" class="title">%s</td></tr>
  <tr><td colspan="2" class="depvar">Dependent variable:</td></tr>
  <tr><td></td><td class="colhead">%s<br>Baseline</td></tr>
  %s
  <tr class="divider"><td style="text-align:left">Observations</td><td style="text-align:center">%d</td></tr>
  <tr><td style="text-align:left">R<sup>2</sup></td><td style="text-align:center">%.3f</td></tr>
  <tr><td style="text-align:left">Adjusted R<sup>2</sup></td><td style="text-align:center">%.3f</td></tr>
  <tr><td style="text-align:left">Residual Std. Error</td><td style="text-align:center">%s</td></tr>
  <tr><td style="text-align:left">F Statistic</td><td style="text-align:center">%s</td></tr>
  <tr class="divider"><td colspan="2" class="note"><em>Note:</em>
    &nbsp;<sup>*</sup>p&lt;0.1; <sup>**</sup>p&lt;0.05; <sup>***</sup>p&lt;0.01</td></tr>
</table>
</body></html>',
    title, dep_label,
    paste(rows, collapse = "\n"),
    n_obs, r2, adj_r2, rse, f_str
  )
}


# ---- Internal transformation helpers --------------------------

#' @noRd
dlog <- function(x) c(NA_real_, diff(log(as.numeric(x))))

#' @noRd
log_level <- function(x) log(as.numeric(x))

#' @noRd
dfirst <- function(x) c(NA_real_, diff(as.numeric(x)))

#' @noRd
pchy <- function(x, n_periods = 12L) {
  x <- as.numeric(x)
  n <- length(x)
  if (n <= n_periods) return(rep(NA_real_, n))
  c(rep(NA_real_, n_periods), x[(n_periods + 1):n] / x[1:(n - n_periods)] - 1)
}

#' @noRd
saar <- function(x) (exp(4 * dlog(x)) - 1) * 100

#' @noRd
apply_trans <- function(x, trans, n_periods = 12L) {
  trans <- trimws(trans)
  if (grepl("^none$",  trans, ignore.case = TRUE)) return(as.numeric(x))
  if (grepl("^PCHY",   trans, ignore.case = TRUE)) return(pchy(x, n_periods))
  if (grepl("^SAAR",   trans, ignore.case = TRUE)) return(saar(x))
  if (grepl("^DLOG",   trans, ignore.case = TRUE)) return(dlog(x))
  if (grepl("^LOG",    trans, ignore.case = TRUE)) return(log_level(x))
  if (grepl("^D\\(",   trans, ignore.case = TRUE)) return(dfirst(x))
  pchy(x, n_periods)   # default fallback
}

# ---- Internal: ARIMA imputation for trailing NAs --------------
#' @noRd
arima_impute_col <- function(x) {
  if (!any(is.na(x))) return(x)
  last_obs <- max(which(!is.na(x)))
  if (last_obs >= length(x)) return(x)
  n_ahead  <- length(x) - last_obs
  obs      <- x[1:last_obs]
  if (sum(!is.na(obs)) < 12) return(x)
  obs[is.na(obs)] <- mean(obs, na.rm = TRUE)
  tryCatch({
    fit <- forecast::auto.arima(obs, stepwise = TRUE, approximation = TRUE)
    fc  <- forecast::forecast(fit, h = n_ahead)
    x[(last_obs + 1):length(x)] <- as.numeric(fc$mean)
    x
  }, error = function(e) x)
}

# ---- Internal: auto-create X_lagN columns from base variable --
#' @noRd
auto_create_lags <- function(vars, dta_trans, verbose = TRUE) {
  missing_v <- setdiff(vars, names(dta_trans))

  for (v in missing_v) {
    lag_match <- regmatches(v, regexpr("_lag(\\d+)$", v))
    if (length(lag_match) == 1) {
      lag_n    <- as.integer(sub("_lag", "", lag_match))
      base_var <- sub("_lag\\d+$", "", v)
      if (base_var %in% names(dta_trans)) {
        dta_trans[[v]] <- dplyr::lag(dta_trans[[base_var]], lag_n)
        if (verbose) message("  Auto-created lag: ", v,
                              " = lag(", base_var, ", ", lag_n, ")")
        missing_v <- setdiff(missing_v, v)
      }
    }
  }

  if (length(missing_v) > 0)
    warning("Variable(s) not found in dta_trans — excluded from model: ",
            paste(missing_v, collapse = ", "), call. = FALSE)

  list(dta_trans = dta_trans, vars = intersect(vars, names(dta_trans)))
}

# ---- Internal: OOS MAE for a fixed variable set ---------------
# base_controls (e.g. AR lags, structural dummies) are forced into every
# trial model in addition to hf_vars; they are never ADD/DROP candidates.
#' @noRd
eval_oos_mae <- function(hf_vars, bridge_full, eval_months, dep_var, base_controls = NULL) {
  if (length(hf_vars) == 0 && length(base_controls) == 0) return(Inf)
  model_vars <- unique(c(base_controls, hf_vars))
  errors <- vapply(eval_months, function(t_month) {
    train_df <- bridge_full[bridge_full$date < t_month, c(dep_var, model_vars), drop = FALSE]
    train_df <- stats::na.omit(train_df)
    if (nrow(train_df) < length(model_vars) + 2) return(NA_real_)
    fml <- stats::as.formula(paste(dep_var, "~", paste(model_vars, collapse = " + ")))
    fit <- tryCatch(stats::lm(fml, data = train_df), error = function(e) NULL)
    if (is.null(fit)) return(NA_real_)
    pred_row <- bridge_full[bridge_full$date == t_month, model_vars, drop = FALSE]
    if (any(is.na(pred_row))) return(NA_real_)
    pred_val <- as.numeric(stats::predict(fit, newdata = pred_row))
    actual   <- bridge_full[bridge_full$date == t_month, dep_var, drop = TRUE]
    pred_val - actual
  }, numeric(1))
  mean(abs(errors), na.rm = TRUE)
}

# ---- Internal: forward-backward stepwise ----------------------
# base_controls stay in every trial/final model; only `candidates` are
# subject to the forward (ADD) / backward (DROP) search.
#' @noRd
run_stepwise_oos_mae <- function(candidates, bridge_full, eval_months,
                                 dep_var, corr_tbl, min_improve = 0.005,
                                 base_controls = NULL) {
  candidates <- setdiff(candidates, base_controls)
  init_var <- corr_tbl$variable[corr_tbl$variable %in% candidates][1]
  if (is.na(init_var)) stop("No valid candidate to initialise stepwise selection.")

  selected    <- init_var
  current_mae <- eval_oos_mae(selected, bridge_full, eval_months, dep_var, base_controls)
  message("  Init: [", init_var, "]  |  MAE = ", round(current_mae, 6))

  # Track history: one row per accepted step
  history <- data.frame(
    step     = 1L,
    action   = "INIT",
    variable = init_var,
    n_vars   = 1L,
    oos_mae  = current_mae,
    stringsAsFactors = FALSE
  )

  improved <- TRUE
  step     <- 1L
  while (improved) {
    improved <- FALSE

    # Forward pass
    remaining <- setdiff(candidates, selected)
    best_add  <- NULL
    best_mae  <- current_mae
    for (v in remaining) {
      trial_mae <- eval_oos_mae(c(selected, v), bridge_full, eval_months, dep_var, base_controls)
      if (!is.na(trial_mae) && trial_mae < best_mae * (1 - min_improve)) {
        best_mae <- trial_mae
        best_add <- v
      }
    }
    if (!is.null(best_add)) {
      step        <- step + 1L
      selected    <- c(selected, best_add)
      current_mae <- best_mae
      message("    ADD  ", formatC(best_add, width = 30, flag = "-"),
              "  OOS MAE = ", round(current_mae, 6))
      history <- rbind(history, data.frame(
        step     = step, action = "ADD", variable = best_add,
        n_vars   = length(selected), oos_mae = current_mae,
        stringsAsFactors = FALSE
      ))
      improved <- TRUE
    }

    # Backward pass
    if (length(selected) > 1) {
      best_drop     <- NULL
      best_mae_drop <- current_mae
      for (v in selected) {
        trial_mae <- eval_oos_mae(setdiff(selected, v), bridge_full, eval_months, dep_var, base_controls)
        if (!is.na(trial_mae) && trial_mae < best_mae_drop * (1 - min_improve)) {
          best_mae_drop <- trial_mae
          best_drop     <- v
        }
      }
      if (!is.null(best_drop)) {
        step        <- step + 1L
        selected    <- setdiff(selected, best_drop)
        current_mae <- best_mae_drop
        message("    DROP ", formatC(best_drop, width = 30, flag = "-"),
                "  OOS MAE = ", round(current_mae, 6))
        history <- rbind(history, data.frame(
          step     = step, action = "DROP", variable = best_drop,
          n_vars   = length(selected), oos_mae = current_mae,
          stringsAsFactors = FALSE
        ))
        improved <- TRUE
      }
    }
  }

  list(vars = unique(selected), oos_mae = current_mae, history = history)
}


# ---- Exported function: np_transform_data ---------------------

#' Transform monthly data and add AR lags and structural dummies
#'
#' Applies per-variable transformations (log-difference, log level,
#' first-difference, year-on-year growth rate, or annualised quarterly
#' growth rate) to the combined monthly SA dataset, then appends AR lags of
#' the dependent variable and optional structural dummy columns.
#'
#' The transformation for each variable is read from \code{trans_map}:
#' \itemize{
#'   \item \code{"DLOG"} or \code{"DLOG(...)"} — log-difference
#'     \eqn{\Delta\ln x_t}
#'   \item \code{"LOG"} or \code{"LOG(...)"} — log level \eqn{\ln x_t}
#'     (no differencing)
#'   \item \code{"D(...)"} — first difference \eqn{\Delta x_t}
#'   \item \code{"PCHY"} — year-on-year growth rate
#'     \eqn{x_t / x_{t-12} - 1}
#'   \item \code{"SAAR"} or \code{"SAAR(...)"} — annualised quarterly growth
#'     rate, percent: \eqn{(\exp(4 \cdot \Delta\ln x_t) - 1) \times 100},
#'     the standard seasonally-adjusted-annualised-rate convention (what the
#'     annual growth rate would be if this quarter's log-growth persisted
#'     for 4 quarters). Intended for quarterly level series (e.g. real GDP);
#'     applying it to data at other frequencies still computes the formula
#'     literally on whatever period-to-period log-difference the series has.
#'   \item \code{"none"} — no transformation
#'   \item anything else / missing — defaults to \code{"PCHY"}
#' }
#'
#' @param combined Data frame as returned by \code{\link{np_merge_monthly}},
#'   \code{\link{np_aggregate_to_target}}, or \code{np_process_data()$combined}.
#'   Must have a \code{date} column giving the anchor month of each period in
#'   \code{"YYYY-MM"} format (e.g. \code{"2020-04"} for 2020Q2 when
#'   \code{target_freq = "quarterly"}).
#' @param trans_map Named character vector mapping variable names (without
#'   \code{_SA} suffix) to transformation codes.  Typically
#'   \code{raw_data$trans_map}.
#' @param dep_var Character. Name of the dependent variable column as it
#'   appears in \code{combined} (including \code{_SA} suffix if present),
#'   e.g. \code{"im_SA"}.
#' @param n_ar_lags Integer. Number of AR lags of \code{dep_var} to add.
#'   Lags 1 through \code{n_ar_lags} plus the seasonal lag implied by
#'   \code{target_freq} (52 for weekly, 12 for monthly, 4 for quarterly;
#'   omitted for annual, where it would duplicate lag 1) are always added.
#'   Defaults to \code{4}.
#' @param target_freq Character. Frequency of \code{combined}: one of
#'   \code{"monthly"}, \code{"quarterly"}, \code{"annual"}, \code{"weekly"}.
#'   Determines the number of periods per year used by the \code{"PCHY"}
#'   (year-on-year) transformation and the seasonal AR lag above. Defaults
#'   to \code{"monthly"}.
#' @param out_dir Character. If supplied, the full transformed dataset
#'   (all candidate variables, AR lags, and date column) is written to
#'   \code{<out_dir>/transformed_data_<dep_var>.csv}.  Defaults to
#'   \code{NULL} (no file written).
#'
#' @return A data frame with:
#' \describe{
#'   \item{\code{date}}{\code{Date} objects, one per period at
#'     \code{target_freq}.}
#'   \item{Transformed value columns}{One \code{_SA}-suffixed column per
#'     variable in \code{combined}, each transformed per \code{trans_map}.}
#'   \item{AR lag columns}{\code{<dep_var>_lag1} … \code{<dep_var>_lag<n>}
#'     and the seasonal lag (see \code{n_ar_lags} above).}
#' }
#'
#' @examples
#' \dontrun{
#' raw    <- np_load_data("Data/", "Fiji")
#' result <- np_process_data("Data/", "Fiji", save_outputs = FALSE)
#'
#' dta_trans <- np_transform_data(
#'   combined   = result$combined,
#'   trans_map  = raw$trans_map,
#'   dep_var    = "im_SA",
#'   n_ar_lags  = 4
#' )
#' head(dta_trans[, 1:6])
#' }
#'
#' @export
np_transform_data <- function(combined,
                              trans_map,
                              dep_var,
                              n_ar_lags   = 4L,
                              target_freq = c("monthly", "quarterly", "annual", "weekly"),
                              out_dir     = NULL) {

  target_freq <- match.arg(target_freq)
  periods_per_year <- switch(target_freq,
    monthly = 12L, quarterly = 4L, annual = 1L, weekly = 52L
  )

  if (!dep_var %in% names(combined))
    stop("'", dep_var, "' not found in combined data frame.")

  # Convert date to Date class. "weekly" combined data is already a full
  # "YYYY-MM-DD" anchor (see np_aggregate_weekly()); other frequencies use
  # a "YYYY-MM" anchor month, so append the first-of-month day.
  dta <- combined
  dta$date <- if (all(nchar(as.character(dta$date)) == 10)) {
    as.Date(dta$date)
  } else {
    as.Date(paste0(dta$date, "-01"))
  }

  # Apply transformations
  value_cols <- setdiff(names(dta), "date")
  dta_trans  <- dta["date"]

  for (col in value_cols) {
    base_name  <- sub("_SA$", "", col)
    trans_code <- if (base_name %in% names(trans_map)) trans_map[[base_name]] else "PCHY"
    dta_trans[[col]] <- apply_trans(dta[[col]], trans_code, n_periods = periods_per_year)
  }

  if (!dep_var %in% names(dta_trans))
    stop("Dependent variable '", dep_var, "' not found after transformation.")

  # AR lags
  for (k in seq_len(n_ar_lags)) {
    dta_trans[[paste0(dep_var, "_lag", k)]] <- dplyr::lag(dta_trans[[dep_var]], k)
  }
  if (periods_per_year > 1) {
    dta_trans[[paste0(dep_var, "_lag", periods_per_year)]] <-
      dplyr::lag(dta_trans[[dep_var]], periods_per_year)
  }

  # Save if out_dir supplied
  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    safe_name <- gsub("[^A-Za-z0-9_]", "_", dep_var)
    out_path  <- file.path(out_dir, paste0("transformed_data_", safe_name, ".csv"))
    utils::write.csv(dta_trans, out_path, row.names = FALSE)
    message("  Saved: ", basename(out_path),
            " (", nrow(dta_trans), " rows x ", ncol(dta_trans), " cols)")
  }

  dta_trans
}


# ---- Exported function: np_rank_correlations ------------------

#' Rank candidate predictors by absolute correlation with the dependent variable
#'
#' Computes Pearson correlation between each candidate predictor and the
#' dependent variable over the full available sample (pairwise complete
#' observations), then returns a data frame ranked by \eqn{|r|} in
#' descending order.
#'
#' @param dta_trans Data frame as returned by \code{\link{np_transform_data}}.
#' @param dep_var Character. Name of the dependent variable column.
#' @param candidates Character vector of candidate predictor column names to
#'   evaluate.  If \code{NULL} (default), all columns except \code{date} and
#'   \code{dep_var} are used.
#'
#' @return A data frame with columns:
#' \describe{
#'   \item{\code{variable}}{Candidate predictor name.}
#'   \item{\code{correlation}}{Pearson \eqn{r} with \code{dep_var}
#'     (signed, full sample).}
#' }
#' Rows are sorted by \eqn{|r|} descending.
#'
#' @examples
#' \dontrun{
#' corr_tbl <- np_rank_correlations(dta_trans, dep_var = "im_SA")
#' head(corr_tbl, 10)
#' }
#'
#' @export
np_rank_correlations <- function(dta_trans, dep_var, candidates = NULL) {

  if (!dep_var %in% names(dta_trans))
    stop("'", dep_var, "' not found in dta_trans.")

  if (is.null(candidates))
    candidates <- setdiff(names(dta_trans), c("date", dep_var))

  y <- dta_trans[[dep_var]]

  corr_vals <- vapply(candidates, function(v) {
    if (!v %in% names(dta_trans)) return(NA_real_)
    r <- stats::cor(y, dta_trans[[v]], use = "pairwise.complete.obs")
    if (!is.finite(r)) 0 else r
  }, numeric(1))

  result <- data.frame(
    variable    = candidates,
    correlation = corr_vals,
    stringsAsFactors = FALSE
  )
  result[order(-abs(result$correlation)), ]
}


# ---- Exported function: np_select_variables -------------------

#' Select bridge model variables via forward-backward stepwise OOS MAE
#'
#' Implements the full variable selection pipeline for the monthly Bridge model:
#'
#' \enumerate{
#'   \item Rank all candidate predictors by \eqn{|r|} with \code{dep_var}.
#'   \item Filter to candidates above \code{corr_threshold}.
#'   \item Add lag-1 and lag-2 of each HF candidate to the pool.
#'   \item Impute trailing NAs in predictors via ARIMA.
#'   \item Run forward-backward stepwise selection minimising OOS MAE on an
#'         expanding estimation window over \code{[oos_start, oos_end]}.
#'   \item Re-estimate the selected model on all available observations.
#'   \item Prune iteratively: drop the highest-\eqn{p} variable if doing so
#'         does not reduce Adjusted \eqn{R^2}.
#' }
#'
#' @param dta_trans Data frame as returned by \code{\link{np_transform_data}}.
#' @param dep_var Character. Dependent variable column name (e.g.
#'   \code{"im_SA"}).
#' @param oos_start Date or character (\code{"YYYY-MM-DD"}). Start of the OOS
#'   evaluation window.
#' @param oos_end Date or character (\code{"YYYY-MM-DD"}). End of the OOS
#'   evaluation window.
#' @param corr_threshold Numeric. Minimum \eqn{|r|} to enter the candidate
#'   pool.  Defaults to \code{0.4}.  All variables in \code{dta_trans}
#'   (including AR lags added by \code{\link{np_transform_data}}) are ranked
#'   and filtered by this threshold.
#' @param min_improve Numeric. Minimum fractional MAE improvement required to
#'   accept an ADD or DROP step.  Defaults to \code{0.005} (0.5\%).
#' @param dummy_vars Character vector of dummy column names to always include
#'   in the candidate pool regardless of correlation.  Defaults to
#'   \code{NULL}.
#' @param base_controls Character vector of column names (typically AR lags
#'   of \code{dep_var} and/or structural break dummies, e.g.
#'   \code{c("im_SA_lag1", "im_SA_lag2", "d_covid")}) that are forced into
#'   every trial and final model rather than being treated as ADD/DROP
#'   candidates. They are excluded from the correlation-ranked candidate
#'   pool and from post-estimation pruning, so they are never dropped.
#'   Defaults to \code{NULL} (no forced regressors; matches prior behaviour).
#' @param out_dir Character. Path to a directory where output files will be
#'   saved.  If \code{NULL} (default), no files are written.  When supplied,
#'   the following files are created:
#'   \itemize{
#'     \item \code{correlation_ranking_<dep_var>.csv} — all candidates ranked
#'       by \eqn{|r|}.
#'     \item \code{oos_mae_history_<dep_var>.csv} — stepwise ADD/DROP history
#'       with OOS MAE at each accepted step.
#'     \item \code{bridge_model_final_<dep_var>.png} — regression table image.
#'     \item \code{bridge_model_final_<dep_var>.xlsx} — model coefficients,
#'       standard errors, p-values, and fit statistics.
#'   }
#' @param verbose Logical. If \code{TRUE} (default), prints stepwise progress
#'   to the console.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{\code{final_model}}{The pruned \code{lm} object estimated on all
#'     available observations.}
#'   \item{\code{selected_vars}}{Character vector of selected predictor names
#'     after pruning.}
#'   \item{\code{oos_mae}}{OOS MAE of the stepwise-selected model (before
#'     post-estimation pruning).}
#'   \item{\code{corr_tbl}}{Data frame of all candidates ranked by \eqn{|r|}.}
#'   \item{\code{oos_history}}{Data frame of the stepwise ADD/DROP history with
#'     OOS MAE at each accepted step.}
#'   \item{\code{spec_tbl}}{Data frame of final model coefficients, standard
#'     errors, p-values, sample info, and fit statistics — ready to save as
#'     CSV.}
#'   \item{\code{train_df}}{The training data frame used to estimate the
#'     final model (listwise deletion on selected variables).}
#' }
#'
#' @examples
#' \dontrun{
#' raw    <- np_load_data("Data/", "Fiji")
#' result <- np_process_data("Data/", "Fiji", save_outputs = FALSE)
#'
#' dta_trans <- np_transform_data(
#'   combined  = result$combined,
#'   trans_map = raw$trans_map,
#'   dep_var   = "im_SA"
#' )
#'
#' sel <- np_select_variables(
#'   dta_trans      = dta_trans,
#'   dep_var        = "im_SA",
#'   oos_start      = as.Date("2020-01-01"),
#'   oos_end        = as.Date("2025-12-01")
#' )
#'
#' summary(sel$final_model)
#' sel$spec_tbl
#' }
#'
#' @export
np_select_variables <- function(dta_trans,
                                dep_var,
                                oos_start,
                                oos_end,
                                corr_threshold = 0.4,
                                min_improve    = 0.005,
                                dummy_vars     = NULL,
                                base_controls  = NULL,
                                out_dir        = NULL,
                                verbose        = TRUE) {

  oos_start <- as.Date(oos_start)
  oos_end   <- as.Date(oos_end)

  if (!dep_var %in% names(dta_trans))
    stop("'", dep_var, "' not found in dta_trans.")

  missing_bc <- setdiff(base_controls, names(dta_trans))
  if (length(missing_bc) > 0)
    stop("base_controls not found in dta_trans: ", paste(missing_bc, collapse = ", "))

  # ---- Step 1: All columns except date, dep_var, and base_controls ----
  all_candidates <- setdiff(names(dta_trans), c("date", dep_var, base_controls))

  # ---- Step 2: Correlation ranking ----
  if (verbose) message("Step 1: Ranking candidates by |r| ...")
  corr_tbl <- np_rank_correlations(dta_trans, dep_var,
                                   candidates = intersect(all_candidates, names(dta_trans)))

  # ---- Step 3: Filter by threshold; always include dummy_vars ----
  pass_threshold  <- corr_tbl$variable[abs(corr_tbl$correlation) > corr_threshold]
  candidates_filtered <- unique(c(pass_threshold, dummy_vars))
  candidates_filtered <- intersect(candidates_filtered, names(dta_trans))

  if (verbose) message("  Candidate pool: ", length(candidates_filtered),
                        " variables (|r| > ", corr_threshold,
                        if (!is.null(dummy_vars)) " + dummies" else "", ")")

  # ---- Step 4: Add lag-1 and lag-2 of non-lag HF candidates only ----
  # AR lags from np_transform_data (e.g. im_SA_lag1) are already in the pool;
  # skip them to avoid creating double-lag columns like im_SA_lag1_lag1.
  hf_base <- candidates_filtered[!grepl("_lag\\d+$", candidates_filtered) &
                                  !candidates_filtered %in% dummy_vars]
  lag_vars <- character(0)
  for (v in hf_base) {
    l1 <- paste0(v, "_lag1"); l2 <- paste0(v, "_lag2")
    dta_trans[[l1]] <- dplyr::lag(dta_trans[[v]], 1)
    dta_trans[[l2]] <- dplyr::lag(dta_trans[[v]], 2)
    lag_vars <- c(lag_vars, l1, l2)
  }
  candidates_filtered <- unique(c(candidates_filtered, lag_vars))
  candidates_filtered <- intersect(candidates_filtered, names(dta_trans))

  if (verbose) message("  Pool after HF lags (lag-1 & lag-2): ",
                        length(candidates_filtered), " variables")

  # ---- Step 5: Build bridge_full with ARIMA imputation ----
  keep_cols   <- intersect(c("date", dep_var, base_controls, candidates_filtered), names(dta_trans))
  bridge_full <- dta_trans[, keep_cols, drop = FALSE]
  bridge_full <- bridge_full[order(bridge_full$date), ]

  pred_cols <- setdiff(keep_cols, c("date", dep_var))
  n_imputed <- 0L
  for (col in pred_cols) {
    if (any(is.na(bridge_full[[col]]))) {
      bridge_full[[col]] <- arima_impute_col(bridge_full[[col]])
      n_imputed <- n_imputed + 1L
    }
  }
  if (verbose) message("  ARIMA imputation applied to ", n_imputed, " column(s).")

  # ---- Step 6: OOS evaluation months ----
  eval_months <- bridge_full$date[
    bridge_full$date >= oos_start &
    bridge_full$date <= oos_end   &
    !is.na(bridge_full[[dep_var]])
  ]
  if (length(eval_months) == 0)
    stop("No OOS evaluation months found between ", oos_start, " and ", oos_end, ".")

  if (verbose) message("Step 2: Forward-backward stepwise OOS MAE (",
                        length(eval_months), " OOS months) ...")

  # ---- Step 7: Stepwise ----
  sw_result <- run_stepwise_oos_mae(
    candidates    = candidates_filtered,
    bridge_full   = bridge_full,
    eval_months   = eval_months,
    dep_var       = dep_var,
    corr_tbl      = corr_tbl,
    min_improve   = min_improve,
    base_controls = base_controls
  )

  if (verbose) {
    message("\n  Best: ", length(sw_result$vars), " variable(s) | OOS MAE = ",
            round(sw_result$oos_mae, 6))
    for (v in sw_result$vars) {
      r <- corr_tbl$correlation[match(v, corr_tbl$variable)]
      message("    ", formatC(v, width = 32, flag = "-"),
              "  r = ", sprintf("%+.3f", ifelse(is.na(r), 0, r)))
    }
  }

  # ---- Step 8: Final model on ALL available observations ----
  if (verbose) message("Step 3: Estimating final model on all observations ...")

  final_vars <- intersect(unique(c(base_controls, sw_result$vars)), names(dta_trans))
  train_df   <- dta_trans[, c("date", dep_var, final_vars), drop = FALSE]
  train_df   <- stats::na.omit(train_df)

  fml              <- stats::as.formula(paste(dep_var, "~", paste(final_vars, collapse = " + ")))
  final_model      <- stats::lm(fml, data = train_df)

  # ---- Step 9: Post-estimation pruning (maximise Adj R²) ----
  # base_controls are exempt: they are forced regressors, never dropped here.
  if (verbose) message("Step 4: Post-estimation pruning (Adj R²) ...")

  pruned <- TRUE
  while (pruned) {
    pruned       <- FALSE
    sm           <- summary(final_model)
    current_adjr <- sm$adj.r.squared
    pvals        <- stats::coef(sm)[, "Pr(>|t|)"]
    droppable    <- pvals[!names(pvals) %in% c("(Intercept)", base_controls)]
    if (length(droppable) == 0) break
    ns_vars      <- droppable[droppable > 0.05]
    if (length(ns_vars) == 0) break
    worst_var    <- names(which.max(ns_vars))
    trial_vars   <- setdiff(names(stats::coef(final_model))[-1], worst_var)
    trial_fml    <- stats::as.formula(paste(dep_var, "~", paste(trial_vars, collapse = " + ")))
    trial_fit    <- tryCatch(stats::lm(trial_fml, data = train_df), error = function(e) NULL)
    if (is.null(trial_fit)) break
    trial_adjr   <- summary(trial_fit)$adj.r.squared
    if (trial_adjr >= current_adjr) {
      if (verbose)
        message("    DROP ", formatC(worst_var, width = 32, flag = "-"),
                "  p = ", sprintf("%.4f", ns_vars[[worst_var]]),
                "  Adj R\u00b2: ", round(current_adjr, 4),
                " -> ", round(trial_adjr, 4))
      final_model <- trial_fit
      pruned      <- TRUE
    } else {
      break
    }
  }

  selected_vars <- names(stats::coef(final_model))[-1]
  final_sm      <- summary(final_model)

  if (verbose) message("  Final Adj R\u00b2 = ", round(final_sm$adj.r.squared, 4),
                        "  (", length(selected_vars), " variable(s))")

  # ---- Step 10: Build spec_tbl ----
  coef_vals <- stats::coef(final_model)
  spec_tbl  <- data.frame(
    term        = names(coef_vals),
    coefficient = as.numeric(coef_vals),
    std_error   = as.numeric(stats::coef(final_sm)[, "Std. Error"]),
    p_value     = as.numeric(stats::coef(final_sm)[, "Pr(>|t|)"]),
    dep_var     = dep_var,
    sample_start = format(min(train_df$date), "%Y-%m"),
    sample_end   = format(max(train_df$date), "%Y-%m"),
    n_obs        = nrow(train_df),
    adj_r2       = final_sm$adj.r.squared,
    oos_mae      = sw_result$oos_mae,
    oos_start    = format(min(eval_months), "%Y-%m"),
    oos_end      = format(max(eval_months), "%Y-%m"),
    stringsAsFactors = FALSE
  )

  if (verbose) message("=== Variable selection complete. ===\n")

  # ---- Step 11: Save outputs if out_dir supplied ----
  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
    safe_name <- gsub("[^A-Za-z0-9_]", "_", dep_var)

    # 1) Correlation table
    corr_path <- file.path(out_dir, paste0("correlation_ranking_", safe_name, ".csv"))
    utils::write.csv(corr_tbl, corr_path, row.names = FALSE)
    if (verbose) message("  Saved: ", basename(corr_path))

    # 2) OOS MAE stepwise history
    hist_path <- file.path(out_dir, paste0("oos_mae_history_", safe_name, ".csv"))
    utils::write.csv(sw_result$history, hist_path, row.names = FALSE)
    if (verbose) message("  Saved: ", basename(hist_path))

    # 3) Regression table as PNG — built with make_reg_html() (no stargazer needed)
    #    webshot2 screenshots the HTML; only webshot2 is an optional dependency.
    sg_html <- file.path(out_dir, paste0("bridge_model_final_", safe_name, ".html"))
    sg_png  <- file.path(out_dir, paste0("bridge_model_final_", safe_name, ".png"))

    html_content <- make_reg_html(
      final_model,
      title     = paste0("Bridge Model \u2014 ", dep_var),
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

    # 4) Specification table as XLSX (requires openxlsx)
    xlsx_path <- file.path(out_dir, paste0("bridge_model_final_", safe_name, ".xlsx"))
    if (requireNamespace("openxlsx", quietly = TRUE)) {
      wb <- openxlsx::createWorkbook()

      # Sheet 1: Coefficients
      openxlsx::addWorksheet(wb, "Coefficients")
      openxlsx::writeDataTable(wb, "Coefficients", spec_tbl, tableStyle = "TableStyleMedium9")
      openxlsx::setColWidths(wb, "Coefficients", cols = seq_len(ncol(spec_tbl)), widths = "auto")

      # Sheet 2: OOS MAE history
      openxlsx::addWorksheet(wb, "OOS_MAE_History")
      openxlsx::writeDataTable(wb, "OOS_MAE_History", sw_result$history,
                               tableStyle = "TableStyleMedium2")
      openxlsx::setColWidths(wb, "OOS_MAE_History",
                             cols = seq_len(ncol(sw_result$history)), widths = "auto")

      # Sheet 3: Correlation ranking
      openxlsx::addWorksheet(wb, "Correlation_Ranking")
      openxlsx::writeDataTable(wb, "Correlation_Ranking", corr_tbl,
                               tableStyle = "TableStyleMedium2")
      openxlsx::setColWidths(wb, "Correlation_Ranking",
                             cols = seq_len(ncol(corr_tbl)), widths = "auto")

      openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)
      if (verbose) message("  Saved: ", basename(xlsx_path))
    } else {
      warning("XLSX output skipped: install 'openxlsx' to enable.", call. = FALSE)
    }
  }

  list(
    final_model   = final_model,
    selected_vars = selected_vars,
    oos_mae       = sw_result$oos_mae,
    corr_tbl      = corr_tbl,
    oos_history   = sw_result$history,
    spec_tbl      = spec_tbl,
    train_df      = train_df
  )
}


# ---- Internal: variable_names -> descriptions lookup ----------
# Reads the 'descriptions' column (2nd column) from the Index file and each
# raw data file that defines it; first source to define a given variable
# name wins. Used to annotate np_baseline_selection()'s saved outputs.
#' @noRd
read_variable_descriptions <- function(data_dir, prefix) {
  raw_dir  <- file.path(data_dir, "Raw")
  desc_map <- character(0)

  add_source <- function(path) {
    if (!file.exists(path)) return(invisible(NULL))
    dta <- tryCatch(read.csv(path, check.names = FALSE, stringsAsFactors = FALSE),
                     error = function(e) NULL)
    if (is.null(dta) || !all(c("variable_names", "descriptions") %in% names(dta)))
      return(invisible(NULL))
    vn   <- trimws(as.character(dta[["variable_names"]]))
    ds   <- as.character(dta[["descriptions"]])
    keep <- !is.na(vn) & nchar(vn) > 0 &
            !is.na(ds) & nchar(trimws(ds)) > 0 &
            !vn %in% names(desc_map)
    if (any(keep)) desc_map[vn[keep]] <<- ds[keep]
    invisible(NULL)
  }

  add_source(file.path(raw_dir, paste0(prefix, "_Index.csv")))
  for (suffix in c("_Daily_Data.csv", "_Daily_Data_Oil.csv", "_Weekly_Data.csv",
                   "_Monthly_Data.csv", "_Quarterly_Data.csv", "_Annual_Data.csv")) {
    add_source(file.path(raw_dir, paste0(prefix, suffix)))
  }

  desc_map
}

# ---- Internal: strip _lagN / _SA suffixes to recover the base variable ----
#' @noRd
strip_variable_suffixes <- function(x) {
  x <- sub("_lag\\d+$", "", x)
  x <- sub("_SA$", "", x)
  x
}


# ---- Exported function: np_baseline_selection -----------------

#' Run the full baseline variable selection pipeline in one call
#'
#' Convenience wrapper that sequentially calls:
#' \enumerate{
#'   \item \code{\link{np_load_data}} — load raw data files and index.
#'   \item \code{\link{np_process_data}} — seasonal adjustment, aggregation,
#'     and merging into a combined monthly dataset.
#'   \item \code{\link{np_transform_data}} — apply per-variable
#'     transformations and add AR lags.
#'   \item \code{\link{np_select_variables}} — correlation ranking,
#'     forward-backward stepwise OOS MAE selection, final OLS, and
#'     post-estimation pruning.
#' }
#'
#' @param data_dir Character. Root data directory containing a \code{Raw/}
#'   sub-folder.  Defaults to \code{"Data/"}.
#' @param prefix Character. File-name prefix (e.g. \code{"Fiji"}).
#'   Defaults to \code{"Fiji"}.
#' @param dep_var Character. Dependent variable column name in the combined
#'   monthly dataset (including \code{_SA} suffix), e.g. \code{"im_SA"}.
#' @param start_date Character. Earliest \code{"YYYY-MM"} month to retain in
#'   the combined dataset.  Defaults to \code{"2000-01"}.
#' @param n_ar_lags Integer. Number of AR lags to add in the transformation
#'   step.  Defaults to \code{4L}.
#' @param target_freq Character. Frequency of \code{dep_var} and the final
#'   modelling dataset: one of \code{"monthly"}, \code{"quarterly"},
#'   \code{"annual"}, \code{"weekly"}. Defaults to \code{"monthly"}. When set
#'   to \code{"quarterly"} or \code{"annual"}, all higher-frequency
#'   predictors are rolled up to that frequency (see
#'   \code{\link{np_aggregate_to_target}}), and the year-on-year transform /
#'   seasonal AR lag adjust accordingly (see \code{\link{np_transform_data}}).
#'   When set to \code{"weekly"}, daily/daily oil series are aggregated
#'   \emph{down} to weekly instead (see \code{\link{np_aggregate_weekly}});
#'   monthly/quarterly/annual source data cannot be split into weeks, so it
#'   is not used for this target.
#' @param oos_start Date or character (\code{"YYYY-MM-DD"}). Start of the OOS
#'   evaluation window for variable selection.
#' @param oos_end Date or character (\code{"YYYY-MM-DD"}). End of the OOS
#'   evaluation window.
#' @param corr_threshold Numeric. Minimum \eqn{|r|} for a variable to enter
#'   the candidate pool.  Defaults to \code{0.4}.
#' @param min_improve Numeric. Minimum fractional OOS MAE improvement to
#'   accept an ADD or DROP step.  Defaults to \code{0.005}.
#' @param dummy_vars Character vector of dummy column names to always include
#'   in the candidate pool.  Defaults to \code{NULL}.
#' @param base_controls Character vector of column names forced into every
#'   trial and final model (e.g. AR lags of \code{dep_var}, structural break
#'   dummies). Excluded from the candidate pool and from post-estimation
#'   pruning, so they are never dropped. Defaults to \code{NULL}. See
#'   \code{\link{np_select_variables}}.
#' @param interpolate_vars Character vector of variable names to interpolate
#'   right after the combined dataset is built (see
#'   \code{\link{np_interpolate_data}}), before transformation. Values of
#'   \code{0} are first treated as missing (see \code{zero_as_na}), then
#'   interior \code{NA} gaps are filled via natural cubic spline. Defaults to
#'   \code{NULL} (no interpolation).
#' @param zero_as_na Logical. If \code{TRUE} (default), values of exactly
#'   \code{0} are treated as missing before interpolating
#'   \code{interpolate_vars}. Ignored if \code{interpolate_vars} is
#'   \code{NULL}.
#' @param out_dir Character. If supplied, all output files from
#'   \code{np_transform_data} and \code{np_select_variables} are written
#'   here, plus:
#'   \itemize{
#'     \item \code{baseline_selection_<dep_var>.rds} — this function's
#'       entire return value. Reload it in a later session with
#'       \code{\link{np_load_baseline_selection}} to run other
#'       \code{np_model_*()} functions on \code{dta_trans} without
#'       re-running this pipeline.
#'     \item \code{variable_descriptions_<dep_var>.csv} — every column in
#'       \code{dta_trans} (\code{variable}), its base variable name with
#'       \code{_lagN}/\code{_SA} suffixes stripped (\code{base_variable}),
#'       and the human-readable \code{description} looked up from the
#'       \code{descriptions} column of the index file / raw data files.
#'     \item \code{baseline_vars_<dep_var>.R} — the final
#'       \code{selection$selected_vars} formatted as a ready-to-source
#'       \code{baseline_vars <- c(...)} script, for quick reuse when
#'       setting up other models.
#'   }
#'   Defaults to \code{NULL} (nothing saved to disk).
#' @param save_processed Logical. If \code{TRUE} (default), processed SA CSVs
#'   are saved to \code{data_dir/Processed/} by \code{np_process_data}.
#' @param verbose Logical. If \code{TRUE} (default), progress messages are
#'   printed throughout the pipeline.
#'
#' @return Invisibly, a named list with elements:
#' \describe{
#'   \item{\code{raw}}{Output of \code{np_load_data}.}
#'   \item{\code{processed}}{Output of \code{np_process_data}.}
#'   \item{\code{dta_trans}}{Transformed data frame from
#'     \code{np_transform_data}.}
#'   \item{\code{selection}}{Output of \code{np_select_variables} — includes
#'     \code{final_model}, \code{selected_vars}, \code{spec_tbl}, etc.}
#' }
#'
#' @examples
#' \dontrun{
#' result <- np_baseline_selection(
#'   data_dir       = "Data/",
#'   prefix         = "Fiji",
#'   dep_var        = "im_SA",
#'   oos_start      = as.Date("2020-01-01"),
#'   oos_end        = as.Date("2025-12-01"),
#'   out_dir        = "Outputs/1_Baseline"
#' )
#' summary(result$selection$final_model)
#' }
#'
#' @export
np_baseline_selection <- function(data_dir       = "Data/",
                                  prefix         = "Fiji",
                                  dep_var        = "im_SA",
                                  start_date     = "2000-01",
                                  n_ar_lags      = 4L,
                                  target_freq    = c("monthly", "quarterly", "annual", "weekly"),
                                  oos_start,
                                  oos_end,
                                  corr_threshold = 0.4,
                                  min_improve    = 0.005,
                                  dummy_vars     = NULL,
                                  base_controls  = NULL,
                                  interpolate_vars = NULL,
                                  zero_as_na       = TRUE,
                                  out_dir        = NULL,
                                  save_processed = TRUE,
                                  verbose        = TRUE) {

  target_freq <- match.arg(target_freq)

  if (!is.null(out_dir) && !dir.exists(out_dir))
    dir.create(out_dir, recursive = TRUE)

  # Step 1: Load raw data
  if (verbose) message("=== NowcastPulse: Baseline Selection Pipeline ===")
  if (verbose) message("Step 1/6: Loading raw data ...")
  raw <- np_load_data(data_dir = data_dir, prefix = prefix)

  # Step 2: Process data (SA + aggregate + merge)
  if (verbose) message("Step 2/6: Processing data (X-13, aggregation, merge) ...")
  processed <- np_process_data(
    data_dir     = data_dir,
    prefix       = prefix,
    start_date   = start_date,
    target_freq  = target_freq,
    save_outputs = save_processed
  )

  # Step 3: Interpolate specified variables (treat 0 as NA, cubic spline)
  if (!is.null(interpolate_vars)) {
    if (verbose) message("Step 3/6: Interpolating specified variables ...")
    processed$combined <- np_interpolate_data(
      combined   = processed$combined,
      vars       = interpolate_vars,
      zero_as_na = zero_as_na,
      verbose    = verbose
    )
  }

  # Step 4: Transform
  if (verbose) message("Step 4/6: Transforming data ...")
  dta_trans <- np_transform_data(
    combined    = processed$combined,
    trans_map   = raw$trans_map,
    dep_var     = dep_var,
    n_ar_lags   = n_ar_lags,
    target_freq = target_freq,
    out_dir     = out_dir
  )

  # Step 5: Variable selection
  if (verbose) message("Step 5/6: Variable selection ...")
  selection <- np_select_variables(
    dta_trans      = dta_trans,
    dep_var        = dep_var,
    oos_start      = oos_start,
    oos_end        = oos_end,
    corr_threshold = corr_threshold,
    min_improve    = min_improve,
    dummy_vars     = dummy_vars,
    base_controls  = base_controls,
    out_dir        = out_dir,
    verbose        = verbose
  )

  # Step 6: Variable descriptions + a ready-to-source baseline_vars script
  safe_name <- gsub("[^A-Za-z0-9_]", "_", dep_var)
  if (!is.null(out_dir)) {
    if (verbose) message("Step 6/6: Saving variable descriptions and baseline_vars script ...")

    desc_map  <- read_variable_descriptions(data_dir, prefix)
    dta_vars  <- setdiff(names(dta_trans), "date")
    desc_tbl  <- data.frame(
      variable      = dta_vars,
      base_variable = strip_variable_suffixes(dta_vars),
      stringsAsFactors = FALSE
    )
    desc_tbl$description <- unname(desc_map[desc_tbl$base_variable])

    desc_path <- file.path(out_dir, paste0("variable_descriptions_", safe_name, ".csv"))
    utils::write.csv(desc_tbl, desc_path, row.names = FALSE, na = "")
    if (verbose) message("  Saved: ", basename(desc_path))

    n_sel <- length(selection$selected_vars)
    baseline_vars_lines <- paste0(
      '  "', selection$selected_vars, '"',
      ifelse(seq_len(n_sel) < n_sel, ",", "")
    )
    script_path <- file.path(out_dir, paste0("baseline_vars_", safe_name, ".R"))
    writeLines(c("baseline_vars <- c(", baseline_vars_lines, ")"), script_path)
    if (verbose) message("  Saved: ", basename(script_path))
  }

  if (verbose) message("=== Baseline selection complete. ===")

  result <- list(
    raw       = raw,
    processed = processed,
    dta_trans = dta_trans,
    selection = selection
  )

  # Persist everything needed to run other np_model_*() functions later
  # without re-running this pipeline (see np_load_baseline_selection()).
  if (!is.null(out_dir)) {
    rds_path <- file.path(out_dir, paste0("baseline_selection_", safe_name, ".rds"))
    saveRDS(result, rds_path)
    if (verbose) message("  Saved: ", basename(rds_path))
  }

  invisible(result)
}


# ---- Exported function: np_load_baseline_selection -------------

#' Reload a previously saved \code{np_baseline_selection()} result
#'
#' Reads back the RDS file written by \code{\link{np_baseline_selection}}
#' (when called with a non-\code{NULL} \code{out_dir}), so that
#' \code{dta_trans} and the variable-selection results can be reused in a
#' later session — e.g. to run \code{\link{np_model_bridge}},
#' \code{\link{np_model_pca}}, \code{\link{np_model_dfm}}, or any other
#' \code{np_model_*()} function — without re-running data loading,
#' processing, transformation, and selection from scratch.
#'
#' @param out_dir Character. The \code{out_dir} that was passed to the
#'   original \code{\link{np_baseline_selection}} call.
#' @param dep_var Character. The \code{dep_var} that was passed to the
#'   original \code{\link{np_baseline_selection}} call (used to build the
#'   saved file name).
#'
#' @return The same named list returned (invisibly) by
#'   \code{\link{np_baseline_selection}}: \code{raw}, \code{processed},
#'   \code{dta_trans}, \code{selection}.
#'
#' @examples
#' \dontrun{
#' # Session 1
#' np_baseline_selection(
#'   data_dir  = "Data/", prefix = "Fiji", dep_var = "im_SA",
#'   oos_start = as.Date("2020-01-01"), oos_end = as.Date("2025-12-01"),
#'   out_dir   = "Outputs/1_Baseline"
#' )
#'
#' # Session 2 (later) — no need to re-run np_baseline_selection()
#' result    <- np_load_baseline_selection("Outputs/1_Baseline", dep_var = "im_SA")
#' dta_trans <- result$dta_trans
#'
#' bridge_out <- np_model_bridge(
#'   dta_trans  = dta_trans,
#'   dep_var    = "im_SA",
#'   model_vars = result$selection$selected_vars,
#'   oos_start  = as.Date("2024-01-01"),
#'   oos_end    = as.Date("2026-04-01")
#' )
#' }
#'
#' @export
np_load_baseline_selection <- function(out_dir, dep_var) {
  safe_name <- gsub("[^A-Za-z0-9_]", "_", dep_var)
  rds_path  <- file.path(out_dir, paste0("baseline_selection_", safe_name, ".rds"))
  if (!file.exists(rds_path))
    stop("No saved baseline selection found at '", rds_path, "'. ",
         "Run np_baseline_selection(..., out_dir = '", out_dir, "') first.")
  readRDS(rds_path)
}
