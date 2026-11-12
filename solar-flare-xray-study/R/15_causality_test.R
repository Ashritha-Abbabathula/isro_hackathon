# 15_causality_test.R -- proves features at time t use no data after t.
# For several cut times T, features are recomputed with ALL data after T deleted; features at times <= T must be identical to the full run.
suppressPackageStartupMessages(library(data.table))
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
proc <- Sys.getenv("SOLAR_PROC", file.path(root, "data", "processed_adl1")); outdir <- Sys.getenv("SOLAR_OUT", file.path(root, "outputs_adl1")); FS <- Sys.getenv("FEATURE_SCRIPT", "10_features_adl1.R"); FF <- Sys.getenv("FEATURE_FILE", "features_adl1.csv"); rs <- file.path(Sys.getenv("R_HOME"), "bin", "Rscript.exe")
to_utc <- function(v) if (inherits(v, "POSIXct")) { attr(v, "tzone") <- "UTC"; v } else as.POSIXct(as.character(v), tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ")
full <- fread(file.path(proc, FF)); full[, time := to_utc(time)]; feat <- grep("^(sxr|hxr|euv|nep|hist|mag)_", names(full), value = TRUE)
cuts <- c("2024-10-20 12:00:00", "2024-11-25 06:30:00", "2024-12-15 18:45:00")
res <- list()
for (T in cuts) {
  Sys.setenv(TRUNC_END = T, SOLAR_PROC = proc); system2(rs, shQuote(file.path(root, "R", FS)), stdout = FALSE, stderr = FALSE); Sys.unsetenv("TRUNC_END")
  tr <- fread(file.path(proc, paste0("features_trunc_", gsub("[^0-9]", "", T), ".csv"))); tr[, time := to_utc(time)]
  tt <- as.POSIXct(T, tz = "UTC"); key <- function(x) as.character(x$time)
  f <- full[time <= tt]; setkey(f, time); setkey(tr, time)
  m <- merge(f[, c("time", feat), with = FALSE], tr[, c("time", feat), with = FALSE], by = "time", suffixes = c(".full", ".trunc"))
  mx <- vapply(feat, function(c) { a <- m[[paste0(c, ".full")]]; b <- m[[paste0(c, ".trunc")]]; ok <- is.finite(a) & is.finite(b); na_mismatch <- sum(is.finite(a) != is.finite(b)); c(max_abs_diff = if (any(ok)) max(abs(a[ok] - b[ok])) else 0, na_mismatch = na_mismatch) }, numeric(2))
  res[[T]] <- data.table(cut_time = T, rows_compared = nrow(m), n_features = length(feat), max_abs_diff_all = max(mx["max_abs_diff", ]), n_features_with_differences = sum(mx["max_abs_diff", ] > 1e-9 | mx["na_mismatch", ] > 0),
                         offending = paste(head(names(which(mx["max_abs_diff", ] > 1e-9 | mx["na_mismatch", ] > 0)), 8), collapse = ","))
  file.remove(file.path(proc, paste0("features_trunc_", gsub("[^0-9]", "", T), ".csv")))
}
R <- rbindlist(res); print(R); dir.create(file.path(outdir, "tables"), FALSE, TRUE); fwrite(R, file.path(outdir, "tables", "causality_test.csv"))
stopifnot(all(R$rows_compared > 1000))   # a vacuous comparison must never count as a pass
message(ifelse(all(R$n_features_with_differences == 0), "CAUSALITY TEST PASSED: no feature depends on future data", "CAUSALITY TEST FAILED: see table"))
