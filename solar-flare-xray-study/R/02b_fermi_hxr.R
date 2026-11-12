# 02b_fermi_hxr.R -- Fermi GBM CTIME (NaI) -> 1-minute hard X-ray count rates
# Output: data/processed/hxr_1min.csv  (time, hxr_lo = ~25-50 keV, hxr_hi = ~50-100 keV, in counts/s, n1+n2 combined)
suppressPackageStartupMessages({ library(FITSio); library(data.table) })
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
hxr_dir <- file.path(root, "data", "raw", "hxr"); proc <- file.path(root, "data", "processed")
MET0 <- as.POSIXct("2001-01-01 00:00:00", tz = "UTC")   # Fermi mission elapsed time origin (leap seconds ignored: <= ~5 s, irrelevant at 1 min)
BANDS <- list(hxr_lo = c(25, 50), hxr_hi = c(50, 100))   # keV

read_day <- function(f) {
  eb <- readFITS(f, hdu = 1); emin <- eb$col[[which(eb$colNames == "E_MIN")]]; emax <- eb$col[[which(eb$colNames == "E_MAX")]]
  sp <- readFITS(f, hdu = 2); g <- function(n) sp$col[[which(sp$colNames == n)]]
  cts <- g("COUNTS"); ex <- g("EXPOSURE"); q <- g("QUALITY"); tm <- g("TIME")
  # a CTIME channel belongs to a band if its centre energy lies inside it
  ctr <- (emin + emax) / 2
  d <- data.table(minute = floor(tm / 60), ex = ex, q = q)
  for (b in names(BANDS)) { ch <- which(ctr >= BANDS[[b]][1] & ctr < BANDS[[b]][2]); d[[b]] <- if (length(ch)) rowSums(cts[, ch, drop = FALSE]) else NA_real_ }
  d <- d[q == 0 & ex > 0]
  d[, .(ex = sum(ex), hxr_lo = sum(hxr_lo), hxr_hi = sum(hxr_hi)), by = minute][, `:=`(det = sub(".*ctime_(n[0-9ab])_.*", "\\1", basename(f)))][]
}

files <- list.files(hxr_dir, pattern = "ctime_n[0-9ab]_[0-9]{6}_v00\\.pha$", full.names = TRUE)
stopifnot(length(files) > 0)
eb <- readFITS(files[1], hdu = 1); message("Energy channels (keV): ",
  paste(sprintf("%.1f-%.1f", eb$col[[which(eb$colNames == "E_MIN")]], eb$col[[which(eb$colNames == "E_MAX")]]), collapse = " | "))
all <- rbindlist(lapply(files, function(f) tryCatch(read_day(f), error = function(e) { message("skip ", basename(f), ": ", e$message); NULL })))
# combine detectors: summed counts / summed exposure; keep only minutes with >= 30 s live exposure
m <- all[, .(ex = sum(ex), hxr_lo = sum(hxr_lo), hxr_hi = sum(hxr_hi), n_det = .N), by = minute][ex >= 30]
m[, `:=`(time = MET0 + minute * 60, hxr_lo = hxr_lo / ex, hxr_hi = hxr_hi / ex)]
out <- m[order(time), .(time, hxr_lo, hxr_hi, n_det)]
fwrite(out, file.path(proc, "hxr_1min.csv"))
message("HXR minutes: ", nrow(out), " covering ", format(min(out$time)), " to ", format(max(out$time)))
