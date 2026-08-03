# ============================================================
# data_processing.R
# Seasonal adjustment (X-13 ARIMA-SEATS), temporal aggregation
# to monthly, merging, and the all-in-one pipeline wrapper.
# ============================================================


# ---- Internal helper: parse_date_by_type ----------------------
#' @noRd
parse_date_by_type <- function(date_vec, dataset_type) {
  x <- as.character(date_vec)
  switch(dataset_type,
    daily     = as.Date(x, format = "%Y%m%d"),
    daily_oil = as.Date(x, format = "%m/%d/%Y"),
    # Weekly dates are stored as "M/YY/CCYY"
    # (e.g. "1/16/2026" = January 2016).
    weekly    = {
      parts     <- strsplit(x, "/")
      month_val <- as.integer(vapply(parts, `[`, character(1L), 1L))
      year2     <- as.integer(vapply(parts, `[`, character(1L), 2L))
      as.Date(sprintf("%04d-%02d-01", 2000L + year2, month_val))
    },
    monthly   = as.Date(paste0(gsub("-", "", x), "01"), format = "%Y%m%d"),
    quarterly = {
      yr <- as.integer(substr(x, 1, 4))
      q  <- as.integer(substr(x, 5, 5))
      as.Date(sprintf("%04d-%02d-01", yr, c(1L, 4L, 7L, 10L)[q]))
    },
    # Annual dates are stored as plain 4-digit years, e.g. "2020".
    annual    = as.Date(paste0(x, "-01-01")),
    stop("Unknown dataset_type: '", dataset_type,
         "'. Must be one of: daily, daily_oil, weekly, monthly, quarterly, annual.")
  )
}


# ---- Internal helper: adjust_seasonal_align -------------------
#' @noRd
adjust_seasonal_align <- function(x, dates, freq) {
  x <- as.numeric(x)
  n <- length(x)

  first_valid <- which(!is.na(x))
  if (length(first_valid) == 0) return(rep(NA_real_, n))
  first_valid <- min(first_valid)
  last_valid  <- max(which(!is.na(x)))

  x_trim <- x[first_valid:last_valid]

  if (sum(!is.na(x_trim)) < max(8, 2 * freq)) return(rep(NA_real_, n))
  if (isTRUE(stats::sd(x_trim, na.rm = TRUE) == 0)) return(x)

  start_date <- dates[first_valid]
  yr  <- as.integer(format(start_date, "%Y"))
  per <- if (freq == 12) {
    as.integer(format(start_date, "%m"))
  } else {
    as.integer((as.integer(format(start_date, "%m")) - 1) / 3) + 1L
  }

  x_ts <- stats::ts(x_trim, start = c(yr, per), frequency = freq)

  fit <- tryCatch(
    seasonal::seas(x_ts, x11 = ""),
    error = function(e) tryCatch(seasonal::seas(x_ts), error = function(e2) NULL)
  )

  if (is.null(fit)) return(x)

  sa_trim <- as.numeric(seasonal::final(fit))

  out <- rep(NA_real_, n)
  out[first_valid:(first_valid + length(sa_trim) - 1L)] <- sa_trim
  out
}


# ---- Internal helper: aggregate a data frame by an arbitrary period key ---
#' @noRd
aggregate_df_by_period <- function(df, period_key, agg_map = NULL) {
  out        <- data.frame(date = sort(unique(period_key)), stringsAsFactors = FALSE)
  value_cols <- setdiff(names(df), "date")

  for (nm in value_cols) {
    v       <- as.numeric(df[[nm]])
    base_nm <- sub("_SA$", "", nm)
    method  <- if (!is.null(agg_map) && base_nm %in% names(agg_map))
                 trimws(agg_map[[base_nm]])
               else "Average"

    agg_fn <- switch(toupper(method),
      "SUM"     = function(z) if (all(is.na(z))) NA_real_ else sum(z,  na.rm = TRUE),
      "AVERAGE" = function(z) if (all(is.na(z))) NA_real_ else mean(z, na.rm = TRUE),
      "LAST"    = function(z) {
                    valid <- z[!is.na(z)]
                    if (length(valid) == 0) NA_real_ else valid[length(valid)]
                  },
                  function(z) if (all(is.na(z))) NA_real_ else mean(z, na.rm = TRUE)
    )

    agg <- tapply(v, period_key, agg_fn)
    out[[nm]] <- as.numeric(agg[out$date])
  }
  out
}


# ---- Internal helper: anchor-month key for a target frequency -------------
# Every target frequency is keyed by the "YYYY-MM" of its period's first
# month, so downstream code (np_transform_data's as.Date(paste0(date,"-01")))
# does not need to know which frequency it is looking at.
#' @noRd
period_anchor_key <- function(dates, target_freq) {
  yr <- as.integer(format(dates, "%Y"))
  if (target_freq == "quarterly") {
    m  <- as.integer(format(dates, "%m"))
    qm <- c(1L, 1L, 1L, 4L, 4L, 4L, 7L, 7L, 7L, 10L, 10L, 10L)[m]
    sprintf("%04d-%02d", yr, qm)
  } else if (target_freq == "annual") {
    sprintf("%04d-01", yr)
  } else {
    format(dates, "%Y-%m")
  }
}


# ---- Internal helper: aggregate_to_monthly_internal -----------
#' @noRd
aggregate_to_monthly_internal <- function(df, dataset_type, agg_map = NULL) {
  dt        <- parse_date_by_type(df$date, dataset_type)
  month_key <- format(dt, "%Y-%m")
  aggregate_df_by_period(df, month_key, agg_map)
}


# ---- Internal helper: apply_x13 -------------------------------
#' @noRd
apply_x13 <- function(df, dataset_type, freq, dataset_label, sa_map = NULL) {
  dates      <- parse_date_by_type(df$date, dataset_type)
  value_cols <- setdiff(names(df), "date")

  message("  [X-13] ", dataset_label, " (", length(value_cols), " variables)")

  for (var in value_cols) {
    base_var <- sub("_SA$", "", var)
    needs_sa <- if (!is.null(sa_map) && base_var %in% names(sa_map))
                  trimws(toupper(sa_map[[base_var]])) == "NSA"
                else TRUE

    if (needs_sa) {
      message("    Adjusting: ", var)
      df[[paste0(var, "_SA")]] <- adjust_seasonal_align(df[[var]], dates = dates, freq = freq)
    } else {
      message("    Pass-through (already SA): ", var)
      df[[paste0(var, "_SA")]] <- df[[var]]
    }
  }

  df[, c("date", paste0(value_cols, "_SA"))]
}


# ---- Internal helper: passthrough_sa (annual data, no seasonality) --------
#' @noRd
passthrough_sa <- function(df, dataset_label) {
  value_cols <- setdiff(names(df), "date")

  message("  [Pass-through] ", dataset_label, " (", length(value_cols),
          " variable(s), annual — no seasonal adjustment)")

  for (var in value_cols) df[[paste0(var, "_SA")]] <- df[[var]]
  df[, c("date", paste0(value_cols, "_SA"))]
}


# ---- Exported function: np_seasonal_adjust --------------------

#' Apply X-13 ARIMA-SEATS seasonal adjustment to all datasets
#'
#' Runs X-13 ARIMA-SEATS (via the \pkg{seasonal} package) on each data frame
#' returned by \code{\link{np_load_data}}.  Variables flagged \code{"NSA"} in
#' \code{sa_map} are adjusted; all others are passed through unchanged.
#' The adjusted columns are renamed with an \code{_SA} suffix. Annual data
#' has no sub-annual seasonality to remove, so it is always passed through
#' unchanged (\code{sa_map} is not consulted for it).
#'
#' Non-standard frequencies (daily \code{freq = 365}, weekly \code{freq = 52})
#' are attempted; if X-13 fails the original series is returned as-is.
#'
#' @param raw_data Named list as returned by \code{\link{np_load_data}}.
#' @param sa_map Named character vector mapping variable names to SA flags.
#'   Typically \code{raw_data$sa_map}. Variables set to \code{"NSA"} will be
#'   seasonally adjusted; all others are passed through.
#'
#' @return A named list with elements \code{daily}, \code{daily_oil},
#'   \code{weekly}, \code{monthly}, \code{quarterly}, \code{annual} — each a
#'   data frame with an \code{_SA}-suffixed column for every variable.
#'
#' @examples
#' \dontrun{
#' raw <- np_load_data(data_dir = "Data/", prefix = "Fiji")
#' sa  <- np_seasonal_adjust(raw, sa_map = raw$sa_map)
#' head(sa$monthly)
#' }
#'
#' @export
np_seasonal_adjust <- function(raw_data, sa_map = raw_data$sa_map) {

  message("Applying X-13 seasonal adjustment ...")

  # Skip any dataset that was not loaded (NULL = file was missing)
  maybe_x13 <- function(df, type, freq, label)
    if (!is.null(df)) apply_x13(df, type, freq, label, sa_map) else NULL
  maybe_passthrough <- function(df, label)
    if (!is.null(df)) passthrough_sa(df, label) else NULL

  list(
    daily     = maybe_x13(raw_data$daily,     "daily",     365, "Daily"),
    daily_oil = maybe_x13(raw_data$daily_oil, "daily_oil", 365, "Daily Oil"),
    weekly    = maybe_x13(raw_data$weekly,    "weekly",    52,  "Weekly"),
    monthly   = maybe_x13(raw_data$monthly,   "monthly",   12,  "Monthly"),
    quarterly = maybe_x13(raw_data$quarterly, "quarterly", 4,   "Quarterly"),
    annual    = maybe_passthrough(raw_data$annual, "Annual")
  )
}


# ---- Exported function: np_aggregate_monthly ------------------

#' Aggregate seasonally adjusted series to monthly frequency
#'
#' Converts daily, daily oil, and weekly SA data frames from
#' \code{\link{np_seasonal_adjust}} to monthly frequency using the per-variable
#' aggregation rules in \code{agg_map} (\code{"Sum"}, \code{"Average"}, or
#' \code{"Last"}).  The monthly SA data frame is standardised to \code{"YYYY-MM"}
#' date format but is otherwise returned unchanged (already monthly).
#'
#' @param sa_data Named list as returned by \code{\link{np_seasonal_adjust}}.
#' @param agg_map Named character vector mapping variable names to aggregation
#'   methods.  Typically \code{raw_data$agg_map}.  Defaults to
#'   \code{"Average"} for any variable not found in the map.
#'
#' @return A named list with elements \code{daily}, \code{daily_oil},
#'   \code{weekly}, \code{monthly} — each a monthly-frequency data frame.
#'   (Quarterly is excluded from the monthly merge workflow.)
#'
#' @examples
#' \dontrun{
#' raw     <- np_load_data(data_dir = "Data/", prefix = "Fiji")
#' sa      <- np_seasonal_adjust(raw)
#' monthly <- np_aggregate_monthly(sa, agg_map = raw$agg_map)
#' head(monthly$daily)
#' }
#'
#' @export
np_aggregate_monthly <- function(sa_data, agg_map = NULL) {

  message("Aggregating SA series to monthly frequency ...")

  # Standardise monthly date format; skip if monthly was not loaded
  monthly_sa_m <- if (!is.null(sa_data$monthly)) {
    m        <- sa_data$monthly
    m$date   <- format(as.Date(paste0(m$date, "01"), format = "%Y%m%d"), "%Y-%m")
    m
  } else NULL

  # Skip NULL datasets (file was missing)
  maybe_agg <- function(df, type)
    if (!is.null(df)) aggregate_to_monthly_internal(df, type, agg_map) else NULL

  list(
    daily     = maybe_agg(sa_data$daily,     "daily"),
    daily_oil = maybe_agg(sa_data$daily_oil, "daily_oil"),
    weekly    = maybe_agg(sa_data$weekly,    "weekly"),
    monthly   = monthly_sa_m
  )
}


# ---- Exported function: np_merge_monthly ----------------------

#' Merge all monthly SA series into a single combined dataset
#'
#' Performs a full outer join on \code{date} across the daily, daily oil,
#' weekly, and monthly SA data frames produced by
#' \code{\link{np_aggregate_monthly}}.  Rows are trimmed to
#' \code{start_date} and sorted chronologically.
#'
#' @param monthly_list Named list as returned by \code{\link{np_aggregate_monthly}}.
#' @param start_date Character. Earliest \code{"YYYY-MM"} date to retain in
#'   the merged output.  Defaults to \code{"2000-01"}.
#'
#' @return A single data frame with column \code{date} (\code{"YYYY-MM"}) and
#'   one column per variable across all four frequency sources.
#'
#' @examples
#' \dontrun{
#' raw      <- np_load_data(data_dir = "Data/", prefix = "Fiji")
#' sa       <- np_seasonal_adjust(raw)
#' monthly  <- np_aggregate_monthly(sa, agg_map = raw$agg_map)
#' combined <- np_merge_monthly(monthly)
#' dim(combined)
#' }
#'
#' @export
np_merge_monthly <- function(monthly_list, start_date = "2000-01") {

  message("Merging monthly SA series ...")

  # Drop any NULL entries (their source file was missing)
  frames <- Filter(
    Negate(is.null),
    list(
      monthly_list$daily,
      monthly_list$daily_oil,
      monthly_list$weekly,
      monthly_list$monthly
    )
  )

  if (length(frames) == 0)
    stop("No monthly data available to merge. Check that at least one data file was loaded.")

  combined <- if (length(frames) == 1) {
    frames[[1]]
  } else {
    Reduce(function(a, b) merge(a, b, by = "date", all = TRUE), frames)
  }

  combined <- combined[order(combined$date), ]
  combined <- combined[combined$date >= start_date, ]
  rownames(combined) <- NULL
  combined
}


# ---- Exported function: np_aggregate_to_target -----------------

#' Aggregate a combined monthly dataset up to a coarser target frequency
#'
#' Rolls the combined monthly dataset from \code{\link{np_merge_monthly}} up
#' to quarterly or annual frequency, using the same per-variable aggregation
#' rules (\code{agg_map}) as the monthly aggregation step. This lets the
#' dependent variable — and therefore the whole modelling pipeline — target
#' quarterly or annual data (e.g. quarterly GDP) while still drawing on
#' higher-frequency (daily/weekly/monthly) predictors, mirroring how a
#' bridge/PCA/DFM model uses partial within-period data to nowcast the
#' current, not-yet-complete period.
#'
#' If \code{target_freq = "monthly"}, \code{combined} is returned unchanged.
#' If native quarterly source data is supplied via \code{quarterly_sa}
#' (\code{sa_data$quarterly} from \code{\link{np_seasonal_adjust}}), it is
#' merged in directly (for \code{target_freq = "quarterly"}) or rolled up to
#' annual first (for \code{target_freq = "annual"}), rather than being
#' dropped as it is in the monthly-only workflow. Likewise, native annual
#' source data supplied via \code{annual_sa} is merged in directly when
#' \code{target_freq = "annual"} (it is too coarse to use otherwise, so it
#' is ignored for \code{"monthly"} / \code{"quarterly"} targets).
#'
#' @param combined Data frame as returned by \code{\link{np_merge_monthly}}.
#' @param target_freq Character. One of \code{"monthly"}, \code{"quarterly"},
#'   \code{"annual"}. Defaults to \code{"monthly"} (no-op).
#' @param agg_map Named character vector mapping variable names to
#'   aggregation methods (\code{"Sum"}, \code{"Average"}, \code{"Last"}).
#'   Typically \code{raw_data$agg_map}.
#' @param quarterly_sa Optional data frame of seasonally adjusted native
#'   quarterly series (\code{sa_data$quarterly}), with \code{date} in the raw
#'   \code{"<year><quarter>"} format produced by \code{\link{np_load_data}}.
#'   Defaults to \code{NULL} (no native quarterly data to merge in).
#' @param annual_sa Optional data frame of seasonally adjusted native annual
#'   series (\code{sa_data$annual}), with \code{date} in the raw 4-digit-year
#'   format produced by \code{\link{np_load_data}}. Only used when
#'   \code{target_freq = "annual"}. Defaults to \code{NULL}.
#'
#' @return A data frame with column \code{date} — the anchor month
#'   (\code{"YYYY-MM"}) of each period, e.g. \code{"2020-01"},
#'   \code{"2020-04"}, ... for quarterly, or \code{"2020-01"}, \code{"2021-01"},
#'   ... for annual — and one column per variable, at \code{target_freq}.
#'
#' @examples
#' \dontrun{
#' raw      <- np_load_data(data_dir = "Data/", prefix = "Fiji")
#' sa       <- np_seasonal_adjust(raw)
#' monthly  <- np_aggregate_monthly(sa, agg_map = raw$agg_map)
#' combined <- np_merge_monthly(monthly)
#' quarterly_combined <- np_aggregate_to_target(
#'   combined, target_freq = "quarterly",
#'   agg_map = raw$agg_map, quarterly_sa = sa$quarterly
#' )
#' annual_combined <- np_aggregate_to_target(
#'   combined, target_freq = "annual",
#'   agg_map = raw$agg_map, quarterly_sa = sa$quarterly, annual_sa = sa$annual
#' )
#' }
#'
#' @export
np_aggregate_to_target <- function(combined,
                                   target_freq = c("monthly", "quarterly", "annual"),
                                   agg_map      = NULL,
                                   quarterly_sa = NULL,
                                   annual_sa    = NULL) {

  target_freq <- match.arg(target_freq)
  if (target_freq == "monthly") return(combined)

  message("Aggregating combined dataset to ", target_freq, " frequency ...")

  dates  <- as.Date(paste0(combined$date, "-01"))
  period <- period_anchor_key(dates, target_freq)
  out    <- aggregate_df_by_period(combined, period, agg_map)

  if (!is.null(quarterly_sa)) {
    q_dates  <- parse_date_by_type(quarterly_sa$date, "quarterly")
    q_period <- period_anchor_key(q_dates, target_freq)
    q_out    <- aggregate_df_by_period(quarterly_sa, q_period, agg_map)
    out <- merge(out, q_out, by = "date", all = TRUE)
  }

  if (target_freq == "annual" && !is.null(annual_sa)) {
    a_dates  <- parse_date_by_type(annual_sa$date, "annual")
    a_period <- period_anchor_key(a_dates, "annual")
    a_out    <- aggregate_df_by_period(annual_sa, a_period, agg_map)
    out <- merge(out, a_out, by = "date", all = TRUE)
  }

  out <- out[order(out$date), ]
  rownames(out) <- NULL
  out
}


# ---- Exported function: np_process_data -----------------------

#' Run the full data processing pipeline
#'
#' Convenience wrapper that sequentially calls \code{\link{np_load_data}},
#' \code{\link{np_seasonal_adjust}}, \code{\link{np_aggregate_monthly}}, and
#' \code{\link{np_merge_monthly}}, then optionally saves all outputs to the
#' \code{Processed/} sub-folder of \code{data_dir}.
#'
#' @param data_dir Character. Path to the root data directory (must contain a
#'   \code{Raw/} sub-folder).  Defaults to \code{"Data/"}.
#' @param prefix Character. Country / project file-name prefix, e.g.
#'   \code{"Fiji"}.  Defaults to \code{"Fiji"}.
#' @param start_date Character. Earliest \code{"YYYY-MM"} date to retain in
#'   the combined output.  Defaults to \code{"2000-01"}.
#' @param target_freq Character. Frequency of the dependent variable / final
#'   modelling dataset: one of \code{"monthly"}, \code{"quarterly"},
#'   \code{"annual"}. Defaults to \code{"monthly"}. When set to
#'   \code{"quarterly"} or \code{"annual"}, all higher-frequency series are
#'   rolled up to that frequency via \code{\link{np_aggregate_to_target}}
#'   (native quarterly source data is merged in rather than dropped; native
#'   annual source data, from an optional \code{<prefix>_Annual_Data.csv},
#'   is merged in too when \code{target_freq = "annual"}).
#' @param save_outputs Logical. If \code{TRUE} (the default), all processed
#'   data frames are written as CSV files to \code{data_dir/Processed/}.
#'
#' @return Invisibly, a named list with elements:
#' \describe{
#'   \item{\code{sa}}{List of seasonally adjusted data frames by frequency
#'     (\code{daily}, \code{daily_oil}, \code{weekly}, \code{monthly},
#'     \code{quarterly}, \code{annual}).}
#'   \item{\code{monthly}}{List of monthly-aggregated SA data frames.}
#'   \item{\code{combined}}{Single combined data frame, at \code{target_freq},
#'     ready for variable selection and modelling.}
#'   \item{\code{target_freq}}{The frequency \code{combined} is at.}
#' }
#'
#' @examples
#' \dontrun{
#' result <- np_process_data(data_dir = "Data/", prefix = "Fiji")
#' head(result$combined)
#'
#' # Quarterly dependent variable (e.g. quarterly GDP):
#' result_q <- np_process_data(data_dir = "Data/", prefix = "Fiji",
#'                             target_freq = "quarterly")
#' }
#'
#' @export
np_process_data <- function(data_dir    = "Data/",
                            prefix      = "Fiji",
                            start_date  = "2000-01",
                            target_freq = c("monthly", "quarterly", "annual"),
                            save_outputs = TRUE) {

  target_freq <- match.arg(target_freq)

  message("=== NowcastPulse: Data Processing Pipeline ===")

  message("Step 1/4: Loading raw data ...")
  raw <- np_load_data(data_dir = data_dir, prefix = prefix)

  message("Step 2/4: Seasonal adjustment (X-13) ...")
  sa  <- np_seasonal_adjust(raw, sa_map = raw$sa_map)

  message("Step 3/4: Aggregating to monthly ...")
  monthly_list <- np_aggregate_monthly(sa, agg_map = raw$agg_map)

  message("Step 4/4: Merging into combined monthly dataset ...")
  combined <- np_merge_monthly(monthly_list, start_date = start_date)

  if (target_freq != "monthly") {
    combined <- np_aggregate_to_target(
      combined, target_freq = target_freq,
      agg_map = raw$agg_map, quarterly_sa = sa$quarterly, annual_sa = sa$annual
    )
  }

  if (save_outputs) {
    proc_dir <- file.path(data_dir, "Processed")
    if (!dir.exists(proc_dir)) dir.create(proc_dir, recursive = TRUE)

    # Only save datasets that were actually loaded (NULL = file was missing)
    save_if_present <- function(df, filename) {
      if (!is.null(df)) {
        write.csv(df, file.path(proc_dir, filename), row.names = FALSE)
        message("  Saved: ", filename)
      } else {
        message("  Skipped (source file missing): ", filename)
      }
    }

    save_if_present(sa$daily,     paste0(prefix, "_Daily_Data_SA.csv"))
    save_if_present(sa$daily_oil, paste0(prefix, "_Daily_Data_Oil_SA.csv"))
    save_if_present(sa$weekly,    paste0(prefix, "_Weekly_Data_SA.csv"))
    save_if_present(sa$monthly,   paste0(prefix, "_Monthly_Data_SA.csv"))
    save_if_present(sa$quarterly, paste0(prefix, "_Quarterly_Data_SA.csv"))
    save_if_present(sa$annual,    paste0(prefix, "_Annual_Data_SA.csv"))

    freq_label <- switch(target_freq,
      monthly = "Monthly", quarterly = "Quarterly", annual = "Annual"
    )
    combined_filename <- paste0(prefix, "_Combined_", freq_label, "_SA.csv")
    write.csv(combined, file.path(proc_dir, combined_filename), row.names = FALSE)
    message("  Saved: ", combined_filename)
    message("Processed files saved to: ", proc_dir)
  }

  message("=== Pipeline complete. ===")
  invisible(list(sa = sa, monthly = monthly_list, combined = combined,
                 target_freq = target_freq))
}
