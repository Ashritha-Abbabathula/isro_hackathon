# 20_features_v2.R -- (v2 of 10_features_adl1.R: adds sub-minute HXR features and a magnetic (SDO/HMI SHARP) block) Aditya-L1 study: merge SoLEXS (SXR), HEL1OS (HXR), GOES EUVS (EUV) and flare history; build causal features + labels
# Input : <SOLAR_PROC>/{solexs_1min,hel1os_1min,goes_1min,flare_events,sep_periods,sharp_fulldisk_12min}.csv
# Output: <SOLAR_PROC>/features_v2.csv
suppressPackageStartupMessages({ library(data.table) })
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
proc <- Sys.getenv("SOLAR_PROC", file.path(root, "data", "processed_v2"))
to_utc <- function(v) if (inherits(v, "POSIXct")) { attr(v, "tzone") <- "UTC"; v } else
  as.POSIXct(as.character(v), tz = "UTC", tryFormats = c("%Y-%m-%dT%H:%M:%SZ", "%Y-%m-%d %H:%M:%S", "%Y-%m-%d %H:%M:%OS", "%Y-%m-%d"))
rd <- function(f) { x <- fread(file.path(proc, f)); for (cn in intersect(c("time", "start", "peak", "end"), names(x))) x[[cn]] <- to_utc(x[[cn]]); x }
floor_min <- function(x) as.POSIXct(floor(as.numeric(x) / 60) * 60, origin = "1970-01-01", tz = "UTC")

sx <- rd("solexs_1min.csv"); hx <- rd("hel1os_1min.csv"); gs <- rd("goes_1min.csv"); ev <- rd("flare_events.csv"); sepdt <- rd("sep_periods.csv")
sx[, time := floor_min(time)]; hx[, time := floor_min(time)]; gs[, time := floor_min(time)]

# ---- grid over the span where SoLEXS and HEL1OS both exist ----
d <- merge(sx, hx, by = "time", all = FALSE); d <- merge(d, gs[, .(time, xrsb_flux, irr_256, irr_284, irr_304)], by = "time", all.x = TRUE)
setorder(d, time)
d <- merge(data.table(time = seq(min(d$time), max(d$time), by = "min")), d, by = "time", all.x = TRUE)
TRUNC <- Sys.getenv("TRUNC_END", ""); if (nzchar(TRUNC)) d <- d[time <= as.POSIXct(TRUNC, tz = "UTC")]   # causality test: delete everything after TRUNC_END
d[, hxr_hi_raw := czt_60_80 + czt_80_150]
setnames(d, c("czt_20_40", "czt_40_60"), c("hxr_lo_raw", "hxr_mid_raw"))

# ---- causal helpers ----
rm_  <- function(x, n) frollmean(x, n, na.rm = TRUE, align = "right")
rsd  <- function(x, n) sqrt(pmax(0, frollmean(x^2, n, na.rm = TRUE, align = "right") - rm_(x, n)^2))
rmax <- function(x, n) frollapply(x, n, function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE), align = "right")
lag_ <- function(x, k) shift(x, k)
chan_feats <- function(x, pre, bg_win = NULL) {
  out <- list(); lx <- log10(pmax(x, 1e-6))
  if (!is.null(bg_win)) { bg <- frollapply(lx, bg_win, median, na.rm = TRUE, align = "right"); out[[paste0(pre, "_excess")]] <- lx - bg }
  else out[[paste0(pre, "_log")]] <- lx
  for (k in c(5, 15, 30)) out[[paste0(pre, "_d", k)]] <- lx - lag_(lx, k)
  out[[paste0(pre, "_acc")]] <- (lx - lag_(lx, 5)) - (lag_(lx, 5) - lag_(lx, 10))
  out[[paste0(pre, "_sd15")]] <- rsd(lx, 15); out[[paste0(pre, "_sd60")]] <- rsd(lx, 60)
  out[[paste0(pre, "_max60")]] <- rmax(lx, 60) - lx
  as.data.table(out)
}
F <- list()
F$sxr <- chan_feats(d$sxr_cps, "sxr_sol")                                          # SoLEXS soft X-ray (Aditya-L1)
F$hxr <- cbind(chan_feats(d$hxr_lo_raw, "hxr_lo", 180), chan_feats(d$hxr_mid_raw, "hxr_mid", 180), chan_feats(d$hxr_hi_raw, "hxr_hi", 180))
F$hxr[, hxr_hardness := log10(pmax(d$hxr_hi_raw, 1e-6)) - log10(pmax(d$hxr_lo_raw, 1e-6))]
F$hxr[, hxr_hardness_d15 := hxr_hardness - lag_(hxr_hardness, 15)]
if ("lo_max8s" %in% names(d)) {   # sub-minute features built on 8 s blocks (the real HEL1OS readout cadence)
  F$hxr[, hxr_lo_peak := log10(pmax(d$lo_max8s, 1e-3)) - log10(pmax(d$hxr_lo_raw, 1e-3))]
  F$hxr[, hxr_hi_peak := log10(pmax(d$hi_max8s, 1e-3)) - log10(pmax(d$hxr_hi_raw, 1e-3))]
  F$hxr[, `:=`(hxr_lo_cv = d$lo_cv8s, hxr_hi_cv = d$hi_cv8s, hxr_lo_nb = d$lo_nburst, hxr_hi_nb = d$hi_nburst)]
  F$hxr[, `:=`(hxr_lo_nb15 = frollsum(d$lo_nburst, 15, na.rm = TRUE, align = "right"), hxr_hi_nb15 = frollsum(d$hi_nburst, 15, na.rm = TRUE, align = "right"))]
  F$hxr[, hxr_lo_peak15 := rmax(hxr_lo_peak, 15)]
}
F$euv <- cbind(chan_feats(d$irr_256, "euv_256"), chan_feats(d$irr_284, "euv_284"), chan_feats(d$irr_304, "euv_304"))
# ---- magnetic block: full-disk SDO/HMI SHARP aggregates (12 min cadence), usable 12 min after observation, carried forward <= 60 min ----
shf <- file.path(proc, "sharp_fulldisk_12min.csv")
if (file.exists(shf)) {
  sh <- fread(shf); sh[, time := as.POSIXct(sub("_TAI", "", T_REC), tz = "UTC", format = "%Y.%m.%d_%H:%M:%S") + 720]; setkey(sh, time)
  mm <- sh[data.table(time = d$time), on = "time", roll = 3600]
  lg <- function(x) log10(pmax(x, 1e-3))
  M <- data.table(mag_usflux = lg(mm$sum_usflux), mag_totusjh = lg(mm$sum_totusjh), mag_totpot = lg(mm$sum_totpot), mag_absnjzh = lg(mm$sum_absnjzh),
                  mag_area = lg(mm$sum_area), mag_rvalue = log10(1 + mm$max_r_value), mag_shr45 = mm$mean_shrgt45, mag_nreg = mm$n_regions)
  for (v in c("mag_usflux", "mag_totusjh", "mag_totpot", "mag_rvalue")) { M[[paste0(v, "_d6h")]] <- M[[v]] - lag_(M[[v]], 360); M[[paste0(v, "_d24h")]] <- M[[v]] - lag_(M[[v]], 1440) }
  F$mag <- M
}
# Neupert block: HXR vs d(SXR)/dt
dsxr <- (d$sxr_cps - lag_(d$sxr_cps, 1)) / pmax(d$sxr_cps, 1)
z <- function(v, n = 180) (v - rm_(v, n)) / (rsd(v, n) + 1e-12)
a <- z(dsxr); b <- z(F$hxr$hxr_lo_excess)
nep <- data.table(nep_dsxr = dsxr, nep_resid = b - a)
nep[, nep_corr15 := frollapply(seq_along(a), 15, function(i) { aa <- a[i]; bb <- b[i]; ok <- is.finite(aa) & is.finite(bb); if (sum(ok) < 8) NA_real_ else cor(aa[ok], bb[ok]) }, align = "right")]
nep[, nep_resid_sd15 := rsd(nep_resid, 15)]
F$nep <- nep

# ---- flare-history block (always available; causal: only flares already finished by t) ----
evC <- ev[class_letter %in% c("C", "M", "X") & !is.na(start) & !is.na(end) & !is.na(peak)]; setorder(evC, start)
tt <- as.numeric(d$time); st <- as.numeric(evC$start); en <- as.numeric(evC$end); lpk <- log10(evC$peak_flux); isM <- evC$class_letter %in% c("M", "X")
hist <- rbindlist(lapply(seq_along(tt), function(i) NULL))                       # placeholder to keep types simple
H <- matrix(NA_real_, length(tt), 5, dimnames = list(NULL, c("hist_n6h", "hist_n24h", "hist_hrs_since", "hist_maxlogflux24h", "hist_nM24h")))
ord <- order(en); en_s <- en[ord]; st_s <- st[ord]; lpk_s <- lpk[ord]; isM_s <- isM[ord]
for (i in seq_along(tt)) {
  t <- tt[i]; k <- findInterval(t, en_s)                                         # flares finished at or before t
  if (k == 0) { H[i, ] <- c(0, 0, 72, -7, 0); next }
  idx <- seq_len(k); e <- en_s[idx]; s <- st_s[idx]
  in6 <- s >= t - 6 * 3600; in24 <- s >= t - 24 * 3600
  H[i, "hist_n6h"] <- sum(in6); H[i, "hist_n24h"] <- sum(in24)
  H[i, "hist_hrs_since"] <- min(72, (t - max(e)) / 3600)
  H[i, "hist_maxlogflux24h"] <- if (any(in24)) max(lpk_s[idx][in24]) else -7
  H[i, "hist_nM24h"] <- sum(in24 & isM_s[idx])
}
F$hist <- as.data.table(H)

X <- do.call(cbind, c(list(d[, .(time)]), unname(F[intersect(c("sxr", "hxr", "euv", "nep", "hist", "mag"), names(F))])))   # unname: named args would prefix the column names

# ---- labels (new flare >= C / >= M starting within H minutes); rows inside an ongoing flare excluded ----
setorder(evC, start); st <- as.numeric(evC$start); en <- as.numeric(evC$end)
nxt <- findInterval(tt, st) + 1; nxt_start <- ifelse(nxt <= length(st), st[pmin(nxt, length(st))], Inf)
X[, ongoing := vapply(tt, function(t) any(st <= t & en >= t), logical(1))]
X[, mins_to_next := (nxt_start - tt) / 60]
for (H_ in c(15, 30, 60)) X[[paste0("y", H_)]] <- as.integer(X$mins_to_next > 0 & X$mins_to_next <= H_)
evM <- evC[class_letter %in% c("M", "X")]; stM <- as.numeric(evM$start); nM <- findInterval(tt, stM) + 1
X[, mins_to_next_M := ifelse(nM <= length(stM), stM[pmin(nM, length(stM))] - tt, Inf) / 60]
for (H_ in c(30, 60)) X[[paste0("yM", H_)]] <- as.integer(X$mins_to_next_M > 0 & X$mins_to_next_M <= H_)

# ---- SEP contamination flag (HEL1OS SEP catalog) ----
X[, sep := FALSE]; for (i in seq_len(nrow(sepdt))) { s0 <- sepdt$start[i]; s1 <- sepdt$end[i]; X[time >= s0 & time <= s1, sep := TRUE] }

feat_cols <- grep("^(sxr|hxr|euv|nep|hist|mag)_", names(X), value = TRUE); stopifnot(length(feat_cols) > 50)   # guard: never compute "complete" over zero columns
X[, complete := !ongoing & complete.cases(.SD), .SDcols = feat_cols]
fwrite(X, file.path(proc, if (nzchar(TRUNC)) paste0("features_trunc_", gsub("[^0-9]", "", TRUNC), ".csv") else "features_v2.csv"))
message("rows: ", nrow(X), " | complete & not-ongoing: ", sum(X$complete), " | of which outside SEP: ", sum(X$complete & !X$sep))
for (H_ in c(15, 30, 60)) message("positives y", H_, " (complete, no SEP): ", sum(X[complete & !sep][[paste0("y", H_)]]))
message("M+ positives yM30/yM60: ", sum(X[complete & !sep]$yM30), "/", sum(X[complete & !sep]$yM60))
message("flares >=C: ", nrow(evC), " (M+: ", nrow(evM), "); SEP minutes flagged: ", sum(X$sep), "; days: ", round(nrow(X) / 1440, 1))
