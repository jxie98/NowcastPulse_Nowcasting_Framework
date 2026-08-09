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
