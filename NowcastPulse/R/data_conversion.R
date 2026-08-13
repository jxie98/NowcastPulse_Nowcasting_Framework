# ============================================================
# data_conversion.R
# Convert an arbitrary raw source file (variable_names,
# descriptions, then date-labelled columns, in any incoming
# date-header format) into the exact wide CSV format that
# np_load_data() expects for a given dataset_type.
# ============================================================


# ---- Internal helper: Excel serial date -> Date --------------------------
#' @noRd
excel_serial_to_date <- function(x) as.Date(as.numeric(x), origin = "1899-12-30")


# ---- Internal helper: parse "2020Q1" / "Q1-2020" style labels ------------
#' @noRd
parse_quarter_label <- function(x) {
  m1 <- regmatches(x, regexec("^(\\d{4})[-_ ]?[Qq]([1-4])$", x))
  yr <- vapply(m1, function(z) if (length(z) == 3) as.integer(z[2]) else NA_integer_, integer(1))
  q  <- vapply(m1, function(z) if (length(z) == 3) as.integer(z[3]) else NA_integer_, integer(1))

  m2  <- regmatches(x, regexec("^[Qq]([1-4])[-_ ]?(\\d{4})$", x))
  yr2 <- vapply(m2, function(z) if (length(z) == 3) as.integer(z[3]) else NA_integer_, integer(1))
  q2  <- vapply(m2, function(z) if (length(z) == 3) as.integer(z[2]) else NA_integer_, integer(1))

  yr[is.na(yr)] <- yr2[is.na(yr)]
  q[is.na(q)]   <- q2[is.na(q)]

  list(year = yr, quarter = q)
}


# ---- Internal helper: flexible date-label parser --------------------------
#' @noRd
parse_flexible_dates <- function(x, dataset_type, date_format = "auto") {
  x <- trimws(as.character(x))
  n <- length(x)

  # Explicit override: caller supplies a strptime format string, or the
  # literal "excel_serial" — skip all auto-detection heuristics.
  if (!identical(date_format, "auto")) {
    if (identical(date_format, "excel_serial")) return(excel_serial_to_date(x))
    parsed <- as.Date(x, format = date_format)
    if (anyNA(parsed))
      stop("date_format = '", date_format, "' failed to parse: ",
           paste(unique(x[is.na(parsed)]), collapse = ", "))
    return(parsed)
  }

  out        <- as.Date(rep(NA, n))
  unresolved <- rep(TRUE, n)
  xn         <- sub("\\.0+$", "", x)  # tolerate "43831.0"-style numerics

  # 1. Quarterly shorthand labels: "2020Q1", "Q1-2020", "2020-Q1" ...
  if (dataset_type == "quarterly") {
    ql  <- parse_quarter_label(x)
    hit <- unresolved & !is.na(ql$year) & !is.na(ql$quarter)
    if (any(hit)) {
      out[hit] <- as.Date(sprintf("%04d-%02d-01", ql$year[hit],
                                   c(1L, 4L, 7L, 10L)[ql$quarter[hit]]))
      unresolved[hit] <- FALSE
    }
  }

  # 2. Plain 4-digit year — only for annual files ("2020")
  hit <- unresolved & dataset_type == "annual" & grepl("^[0-9]{4}$", xn)
  if (any(hit)) {
    out[hit] <- as.Date(paste0(xn[hit], "-01-01"))
    unresolved[hit] <- FALSE
  }

  # 3. Bare "YYYYMM" — only for monthly files ("202001")
  hit <- unresolved & dataset_type == "monthly" & grepl("^[0-9]{6}$", xn)
  if (any(hit)) {
    parsed <- as.Date(paste0(xn[hit], "01"), format = "%Y%m%d")
    ok <- hit; ok[hit] <- !is.na(parsed)
    out[ok] <- parsed[!is.na(parsed)]
    unresolved[ok] <- FALSE
  }

  # 4. Bare "YYYYMMDD" (daily), e.g. "20100103"
  hit <- unresolved & grepl("^[0-9]{8}$", xn)
  if (any(hit)) {
    parsed <- as.Date(xn[hit], format = "%Y%m%d")
    ok <- hit; ok[hit] <- !is.na(parsed)
    out[ok] <- parsed[!is.na(parsed)]
    unresolved[ok] <- FALSE
  }

  # 5. Excel serial date numbers — remaining short pure-numeric labels
  #    (4-5 digits; 6/8-digit numerics were already claimed above)
  hit <- unresolved & grepl("^[0-9]{4,5}$", xn)
  if (any(hit)) {
    parsed <- excel_serial_to_date(xn[hit])
    ok <- hit; ok[hit] <- !is.na(parsed)
    out[ok] <- parsed[!is.na(parsed)]
    unresolved[ok] <- FALSE
  }

  # 6a. Full calendar-date text formats. R's as.Date() does not enforce
  # %Y to be exactly 4 digits, so e.g. "1/1/2020" can silently (and
  # wrongly) match format = "%Y/%m/%d" as well as the intended
  # "%m/%d/%Y" — guard each format with a shape regex so only the
  # correctly-shaped candidates are attempted.
  full_pairs <- list(
    list(regex = "^[0-9]{4}-[0-9]{1,2}-[0-9]{1,2}$", fmt = "%Y-%m-%d"),
    list(regex = "^[0-9]{4}/[0-9]{1,2}/[0-9]{1,2}$", fmt = "%Y/%m/%d"),
    list(regex = "^[0-9]{1,2}/[0-9]{1,2}/[0-9]{4}$", fmt = "%m/%d/%Y"),
    list(regex = "^[0-9]{1,2}-[0-9]{1,2}-[0-9]{4}$", fmt = "%m-%d-%Y")
  )
  for (p in full_pairs) {
    if (!any(unresolved)) break
    cand <- unresolved & grepl(p$regex, x)
    if (!any(cand)) next
    parsed <- suppressWarnings(as.Date(x[cand], format = p$fmt))
    ok <- which(cand)[!is.na(parsed)]
    out[ok] <- parsed[!is.na(parsed)]
    unresolved[ok] <- FALSE
  }

  # 6b. Day-month-year with a month name, e.g. "16-Jan-2020", "16 January 2020"
  dmy_formats <- c("%d-%b-%Y", "%d-%B-%Y", "%d %b %Y", "%d %B %Y")
  for (fmt in dmy_formats) {
    if (!any(unresolved)) break
    parsed <- suppressWarnings(as.Date(x[unresolved], format = fmt))
    ok <- which(unresolved)[!is.na(parsed)]
    out[ok] <- parsed[!is.na(parsed)]
    unresolved[ok] <- FALSE
  }

  # 6c. Month-year only, month name (assumes the 1st of the month), e.g.
  #     "Jan-2020", "January 2020"
  my_formats <- c("%b-%Y", "%B-%Y", "%b %Y", "%B %Y")
  for (fmt in my_formats) {
    if (!any(unresolved)) break
    parsed <- suppressWarnings(as.Date(paste0(x[unresolved], "-01"),
                                        format = paste0(fmt, "-%d")))
    ok <- which(unresolved)[!is.na(parsed)]
    out[ok] <- parsed[!is.na(parsed)]
    unresolved[ok] <- FALSE
  }

  # 6d. Month-year only, numeric ("2020-01", "2020/01") — shape-guarded
  # for the same reason as 6a.
  my_num_pairs <- list(
    list(regex = "^[0-9]{4}-[0-9]{1,2}$", fmt = "%Y-%m"),
    list(regex = "^[0-9]{4}/[0-9]{1,2}$", fmt = "%Y/%m")
  )
  for (p in my_num_pairs) {
    if (!any(unresolved)) break
    cand <- unresolved & grepl(p$regex, x)
    if (!any(cand)) next
    parsed <- suppressWarnings(as.Date(paste0(x[cand], "-01"), format = paste0(p$fmt, "-%d")))
    ok <- which(cand)[!is.na(parsed)]
    out[ok] <- parsed[!is.na(parsed)]
    unresolved[ok] <- FALSE
  }

  if (any(unresolved))
    stop("Could not parse date label(s) for dataset_type = '", dataset_type, "': ",
         paste(unique(x[unresolved]), collapse = ", "),
         "\nPass an explicit `date_format` (a strptime format string, e.g. \"%d.%m.%Y\", ",
         "or the literal \"excel_serial\").")

  out
}


# ---- Internal helper: reformat a Date into the target header string ------
#' @noRd
format_date_header <- function(d, dataset_type) {
  switch(dataset_type,
    daily      = format(d, "%Y%m%d"),
    daily_oil  = ,
    weekly     = format(d, "%m/%d/%Y"),
    monthly    = format(d, "%Y%m"),
    quarterly  = {
      yr <- as.integer(format(d, "%Y"))
      q  <- (as.integer(format(d, "%m")) - 1L) %/% 3L + 1L
      sprintf("%04d%d", yr, q)
    },
    annual     = format(d, "%Y"),
    stop("Unknown dataset_type: '", dataset_type, "'.")
  )
}


# ---- Exported function: np_convert_raw_data -------------------------------

#' Convert a raw source file into the package's required data-file format
#'
#' Many source extracts already have the right shape (\code{variable_names}
#' in the first column, a description in the second, and one column per
#' date) but use whatever date-header convention the source system happens
#' to export — \code{"1/16/2020"}, \code{"2020-01-16"}, \code{"16-Jan-2020"},
#' an Excel serial number, and so on. \code{np_convert_raw_data} rewrites
#' those date headers into the exact format \code{\link{np_load_data}}
#' requires for the given \code{dataset_type} (see \code{vignette("introduction")}
#' for the full table), and writes the result to \code{<data_dir>/Raw/} under
#' the package's standard file-naming convention, ready to be picked up by
#' \code{\link{np_load_data}} / \code{\link{np_process_data}}.
#'
#' Column position — not column name — determines role: whatever the source
#' file calls its first two columns, they are treated as \code{variable_names}
#' and \code{descriptions} and renamed accordingly in the output. Date labels
#' are auto-detected (Excel serial numbers, \code{"YYYYMMDD"}/\code{"YYYYMM"}/
#' \code{"YYYYQ"} shorthand, and common text formats such as
#' \code{"M/D/YYYY"}, \code{"YYYY-MM-DD"}, \code{"DD-Mon-YYYY"},
#' \code{"Mon-YYYY"}, \code{"YYYY-MM"}, and \code{"YYYYQ#"} quarter labels).
#' If auto-detection fails or is ambiguous, pass \code{date_format} explicitly.
#'
#' @param input_file Character. Path to the raw source CSV.
#' @param dataset_type Character. One of \code{"daily"}, \code{"daily_oil"},
#'   \code{"weekly"}, \code{"monthly"}, \code{"quarterly"}, \code{"annual"}.
#'   Determines both the required output date format and the output file name.
#' @param data_dir Character. Root data directory; output is written to its
#'   \code{Raw/} sub-folder (created if missing). Defaults to \code{"Data/"}.
#' @param prefix Character. Country / project file-name prefix used in the
#'   output file name, e.g. \code{"Fiji"} -> \code{Fiji_Monthly_Data.csv}.
#'   Defaults to \code{"Fiji"}.
#' @param date_col_start Integer. Column index where date-labelled columns
#'   begin. Defaults to \code{3} (i.e. columns 1-2 are variable_names /
#'   descriptions, matching the package's standard layout).
#' @param date_format Character. \code{"auto"} (default) tries a series of
#'   common formats/heuristics. Alternatively pass an explicit
#'   \code{\link[base]{strptime}} format string (e.g. \code{"\%d.\%m.\%Y"}), or
#'   the literal \code{"excel_serial"} to force Excel serial-number parsing.
#'   Use this to disambiguate when auto-detection errors or guesses wrong.
#' @param overwrite Logical. If \code{TRUE} (default), replaces an existing
#'   output file. \code{Raw/} is treated as a build output regenerated from
#'   \code{input_file}, so overwriting it is the expected, safe case; set to
#'   \code{FALSE} to refuse instead.
#'
#' @return Invisibly, the path to the written output file.
#'
#' @examples
#' \dontrun{
#' np_convert_raw_data(
#'   input_file   = "Data/Source/Fiji_Monthly_Extract.csv",
#'   dataset_type = "monthly",
#'   data_dir     = "Data/",
#'   prefix       = "Fiji"
#' )
#'
#' # Ambiguous source dates (e.g. Excel export) — force interpretation:
#' np_convert_raw_data(
#'   input_file   = "Data/Source/Fiji_Weekly_Extract.csv",
#'   dataset_type = "weekly",
#'   date_format  = "excel_serial"
#' )
#' }
#'
#' @export
np_convert_raw_data <- function(input_file,
                                dataset_type   = c("daily", "daily_oil", "weekly",
                                                    "monthly", "quarterly", "annual"),
                                data_dir       = "Data/",
                                prefix         = "Fiji",
                                date_col_start = 3,
                                date_format    = "auto",
                                overwrite      = TRUE) {

  dataset_type <- match.arg(dataset_type)

  if (!file.exists(input_file))
    stop("input_file not found: ", input_file)

  df <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)

  if (ncol(df) < date_col_start)
    stop("input_file has ", ncol(df), " column(s) but date_col_start = ",
         date_col_start, " — expected variable_names, descriptions, then dates.")

  date_headers <- names(df)[date_col_start:ncol(df)]
  parsed_dates <- parse_flexible_dates(date_headers, dataset_type, date_format)

  new_headers <- format_date_header(parsed_dates, dataset_type)

  dup <- new_headers[duplicated(new_headers)]
  if (length(dup) > 0)
    warning("Converted date headers contain duplicate period(s) after formatting: ",
            paste(unique(dup), collapse = ", "),
            " — check date_format / dataset_type are correct for this file.",
            call. = FALSE)

  # Order columns chronologically; keep variable_names / descriptions first
  ord <- order(parsed_dates)
  df  <- df[, c(seq_len(date_col_start - 1), (date_col_start:ncol(df))[ord]), drop = FALSE]
  names(df) <- c("variable_names", "descriptions", new_headers[ord])

  file_map <- c(
    daily     = "Daily_Data.csv",
    daily_oil = "Daily_Data_Oil.csv",
    weekly    = "Weekly_Data.csv",
    monthly   = "Monthly_Data.csv",
    quarterly = "Quarterly_Data.csv",
    annual    = "Annual_Data.csv"
  )

  out_dir  <- file.path(data_dir, "Raw")
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  out_file <- file.path(out_dir, paste0(prefix, "_", file_map[[dataset_type]]))

  if (file.exists(out_file) && !overwrite)
    stop("Output file already exists: ", out_file, "\nPass overwrite = TRUE to replace it.")

  write.csv(df, out_file, row.names = FALSE)
  message("Converted '", input_file, "' -> '", out_file, "' (", length(new_headers),
          " date column(s), dataset_type = '", dataset_type, "')")

  invisible(out_file)
}
