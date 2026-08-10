# NowcastPulse 0.1.4

## Enhancements

* `np_baseline_selection()`: when `out_dir` is supplied, its entire return
  value (`raw`, `processed`, `dta_trans`, `selection`) is now also saved as
  a single `baseline_selection_<dep_var>.rds` file.
* Added `np_load_baseline_selection(out_dir, dep_var)` to reload that file
  in a later session. Since every `np_model_*()` function only needs
  `dta_trans` (plus `dep_var`/variable names/OOS dates), users can now run
  `np_baseline_selection()` once, then reload its saved output and go
  straight to `np_model_bridge()`, `np_model_pca()`, `np_model_dfm()`, etc.
  without re-running data loading, processing, transformation, and
  variable selection every time.

# NowcastPulse 0.1.3

## Enhancements

* `np_select_variables()` / `np_baseline_selection()`: added a `base_controls`
  argument. Columns listed there (typically AR lags of `dep_var` and/or
  structural break dummies) are forced into every trial model during the
  forward-backward stepwise search and into the final model, rather than
  being treated as ADD/DROP candidates. They are excluded from the
  correlation-ranked candidate pool and exempt from post-estimation
  pruning, so they can never be dropped. Aligns the package's variable
  selection with the fixed-base-controls bridge-model methodology used in
  the India/US GDP nowcasting demo scripts. Defaults to `NULL`, so existing
  calls are unaffected.

# NowcastPulse 0.1.2

## Enhancements

* Internal `adjust_seasonal_align()`: when X-13/SEATS identifies a degenerate
  non-seasonal `(0 0 0)` ARIMA model, the series is now refit with a
  standard airline model (`(0 1 1)(0 1 1)`) so it still gets a real
  seasonal adjustment instead of silently falling through to the raw,
  unadjusted series. The unadjusted-series fallback is kept as a last
  resort if the refit still fails or produces a mismatched-length result.

# NowcastPulse 0.1.1

## Bug fixes

* `apply_x13()` / internal `adjust_seasonal_align()`: fixed a crash
  ("replacement has length zero") that occurred when X-13/SEATS selected a
  degenerate ARIMA(0,0,0) model and `seasonal::final()` silently returned a
  zero-length vector. The function now falls back to the unadjusted series
  whenever the adjusted output isn't the expected length, instead of
  computing a malformed (descending) assignment index. Reproduced with FRED
  series IPDCONGD / IPNCONGD, but not specific to those series.

* Internal `parse_date_by_type(x, "weekly")`: fixed misparsing of full
  `M/D/YYYY` weekly date strings (e.g. `"1/16/2010"`). The 2-token `"M/YY"`
  parser was reading the day-of-month as a 2-digit year, scattering
  consecutive weekly observations across different decades. The `weekly`
  branch now parses `M/D/YYYY` first (matching `daily_oil`) and falls back
  to the legacy 2-token `"M/YY"` format only when that fails.
