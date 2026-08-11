# ============================================================
# data_loading.R
# Functions for reading raw multi-frequency data files and the
# variable-level index into R data frames.
# ============================================================


# ---- Internal helper: wide_to_long ----------------------------
#' @noRd
wide_to_long <- function(dta) {
  var_names <- dta[["variable_names"]]
  dates     <- names(dta)[-(1:2)]

  keep      <- !is.na(var_names) & nchar(trimws(var_names)) > 0
  var_names <- var_names[keep]
  dta       <- dta[keep, ]

  mat <- t(as.matrix(dta[, -(1:2)]))
  df  <- as.data.frame(mat, stringsAsFactors = FALSE)
  colnames(df) <- var_names

  df[] <- lapply(df, function(x) {
    x[x == "#N/A"] <- NA
    suppressWarnings(as.numeric(x))
  })

  cbind(date = dates, df, stringsAsFactors = FALSE)
}


# ---- Exported function: np_load_data --------------------------

#' Load raw multi-frequency data files
#'
#' Reads all raw CSV files (daily, daily oil, weekly, monthly, quarterly,
#' annual) from the \code{raw/} sub-folder of \code{data_dir} and reshapes
#' each from the project-standard wide format (variables × dates) to a long
#' format (dates × variables).  Also reads the variable-level index file
#' which controls seasonal-adjustment and aggregation behaviour.
#'
#' @section Expected file layout:
#' \code{data_dir} must contain a \code{raw/} sub-folder with:
#' \itemize{
#'   \item \code{<prefix>_Daily_Data.csv}
#'   \item \code{<prefix>_Daily_Data_Oil.csv}
#'   \item \code{<prefix>_Weekly_Data.csv}
#'   \item \code{<prefix>_Monthly_Data.csv}
#'   \item \code{<prefix>_Quarterly_Data.csv}
#'   \item \code{<prefix>_Annual_Data.csv} — date labels are plain 4-digit
#'     years, e.g. \code{"2020"}. Only useful when \code{target_freq =
#'     "annual"} in \code{\link{np_process_data}} / \code{\link{np_baseline_selection}}
#'     — it is too coarse to roll down to monthly or quarterly.
#'   \item \code{<prefix>_Index.csv}
#' }
#' Each data CSV must have a \code{variable_names} column in the first
#' position, a \code{descriptions} column in the second position, and
#' date labels as the remaining column headers.
#'
#' The index CSV must have columns \code{variable_names}, \code{agg_index}
#' (\code{"Sum"} / \code{"Average"} / \code{"Last"}), and \code{sa_index}
#' (\code{"NSA"} = needs adjustment, anything else = already SA).
#'
#' @param data_dir Character. Path to the root data directory (must contain a
#'   \code{raw/} sub-folder). Defaults to \code{"Data/"}.
#' @param prefix Character. Country / project prefix used in file names,
#'   e.g. \code{"Fiji"}. Defaults to \code{"Fiji"}.
#'
#' @return A named list with elements:
#' \describe{
#'   \item{\code{daily}}{Data frame — daily series (dates × variables).}
#'   \item{\code{daily_oil}}{Data frame — daily oil price series.}
#'   \item{\code{weekly}}{Data frame — weekly series.}
#'   \item{\code{monthly}}{Data frame — monthly series.}
#'   \item{\code{quarterly}}{Data frame — quarterly series.}
#'   \item{\code{annual}}{Data frame — annual series.}
#'   \item{\code{agg_map}}{Named character vector mapping variable names to
#'     aggregation methods (\code{"Sum"}, \code{"Average"}, \code{"Last"}).}
#'   \item{\code{sa_map}}{Named character vector mapping variable names to SA
#'     flags (\code{"NSA"} or otherwise).}
#'   \item{\code{trans_map}}{Named character vector mapping variable names to
#'     transformation codes (\code{"DLOG"}, \code{"LOG"}, \code{"D(...)"},
#'     \code{"PCHY"}, \code{"none"}).}
#' }
#'
#' @examples
#' \dontrun{
#' raw <- np_load_data(data_dir = "Data/", prefix = "Fiji")
#' head(raw$monthly)
#' raw$agg_map[1:5]
#' }
#'
#' @export
np_load_data <- function(data_dir = "Data/", prefix = "Fiji") {

  raw_dir <- file.path(data_dir, "Raw")

  # ---- Index file is required (no maps = no pipeline) -----------------------
  index_file <- paste0(prefix, "_Index.csv")
  index_path <- file.path(raw_dir, index_file)
  if (!file.exists(index_path))
    stop("Index file not found: ", index_path, "\n",
         "This file is required to supply agg_index, sa_index, and trans_index.")

  fiji_index <- read.csv(index_path, check.names = FALSE, stringsAsFactors = FALSE)

  required_idx_cols <- c("variable_names", "agg_index", "sa_index", "trans_index")
  missing_idx_cols  <- setdiff(required_idx_cols, names(fiji_index))
  if (length(missing_idx_cols) > 0)
    stop("Index file '", index_file, "' is missing column(s): ",
         paste(missing_idx_cols, collapse = ", "))

  # ---- Data files are optional — warn and return NULL for any that are missing
  data_files <- c(
    daily     = paste0(prefix, "_Daily_Data.csv"),
    daily_oil = paste0(prefix, "_Daily_Data_Oil.csv"),
    weekly    = paste0(prefix, "_Weekly_Data.csv"),
    monthly   = paste0(prefix, "_Monthly_Data.csv"),
    quarterly = paste0(prefix, "_Quarterly_Data.csv"),
    annual    = paste0(prefix, "_Annual_Data.csv")
  )

  missing_data <- names(data_files)[!file.exists(file.path(raw_dir, data_files))]
  if (length(missing_data) > 0) {
    warning(
      length(missing_data), " data file(s) not found — skipped:\n",
      paste0("  - ", file.path(raw_dir, data_files[missing_data]), collapse = "\n"),
      "\nThe pipeline will continue with the available files.",
      call. = FALSE
    )
  }

  available <- names(data_files)[file.exists(file.path(raw_dir, data_files))]
  if (length(available) == 0)
    stop("No data files found in '", raw_dir, "'. At least one data CSV is required.")

  read_wide <- function(filename) {
    path <- file.path(raw_dir, filename)
    wide_to_long(read.csv(path, check.names = FALSE, stringsAsFactors = FALSE))
  }

  # Return NULL for missing datasets so downstream functions can skip them
  load_or_null <- function(key) {
    if (key %in% available) read_wide(data_files[[key]]) else NULL
  }

  list(
    daily     = load_or_null("daily"),
    daily_oil = load_or_null("daily_oil"),
    weekly    = load_or_null("weekly"),
    monthly   = load_or_null("monthly"),
    quarterly = load_or_null("quarterly"),
    annual    = load_or_null("annual"),
    agg_map   = stats::setNames(fiji_index[["agg_index"]],   fiji_index[["variable_names"]]),
    sa_map    = stats::setNames(fiji_index[["sa_index"]],    fiji_index[["variable_names"]]),
    trans_map = stats::setNames(fiji_index[["trans_index"]], fiji_index[["variable_names"]])
  )
}
