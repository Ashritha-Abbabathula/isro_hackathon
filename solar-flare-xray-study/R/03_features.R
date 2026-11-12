# 03_features.R -- merge bands, engineer causal features per block, build forecast labels
# Output: data/processed/features.csv  (one row per minute; columns prefixed sxr_ / hxr_ / euv_ / nep_ ; labels y15,y30,y60)
suppressPackageStartupMessages({ library(data.table) })
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
proc <- file.path(root, "data", "processed")
to_utc <- function(v) if (inherits(v, "POSIXct")) { attr(v, "tzone") <- "UTC"; v } else
  as.POSIXct(as.character(v), tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ")
rd <- function(f) { x <- fread(file.path(proc, f)); for (cn in intersect(c("time", "start", "peak", "end"), names(x))) x[[cn]] <- to_utc(x[[cn]]); x }
goes <- rd("goes_1min.csv"); hxr <- rd("hxr_1min.csv"); ev <- rd("flare_events.csv")
goes[, time := as.POSIXct(floor(as.numeric(time) / 60) * 60, origin = "1970-01-01", tz = "UTC")]
hxr[, time := as.POSIXct(floor(as.numeric(time) / 60) * 60, origin = "1970-01-01", tz = "UTC")]

# ---- common 1-min grid restricted to the span where all three bands exist ----
d <- merge(goes, hxr[, .(time, hxr_lo, hxr_hi)], by = "time", all = FALSE)
setorder(d, time)
full <- data.table(time = seq(min(d$time), max(d$time), by = "min"))
d <- merge(full, d, by = "time", all.x = TRUE)

# ---- causal helpers (right-aligned windows; only past data) ----
rm_  <- function(x, n) frollmean(x, n, na.rm = TRUE, align = "right")
rsd  <- function(x, n) sqrt(pmax(0, frollmean(x^2, n, na.rm = TRUE, align = "right") - rm_(x, n)^2))
rmax <- function(x, n) frollapply(x, n, function(v) if (all(is.na(v))) NA_real_ else max(v, na.rm = TRUE), align = "right")
lag_ <- function(x, k) shift(x, k)
# log-flux channels: SXR (GOES XRS), EUV lines (GOES EUVS) ; HXR: log count rate
chan_feats <- function(x, pre, bg_win = NULL) {
  out <- list(); lx <- log10(x)
  if (!is.null(bg_win)) { bg <- frollapply(lx, bg_win, median, na.rm = TRUE, align = "right"); out[[paste0(pre, "_excess")]] <- lx - bg }  # level above running background
  else out[[paste0(pre, "_log")]] <- lx
  for (k in c(5, 15, 30)) out[[paste0(pre, "_d", k)]] <- lx - lag_(lx, k)            # change over past k minutes
  out[[paste0(pre, "_acc")]] <- (lx - lag_(lx, 5)) - (lag_(lx, 5) - lag_(lx, 10))       # curvature
  out[[paste0(pre, "_sd15")]] <- rsd(lx, 15); out[[paste0(pre, "_sd60")]] <- rsd(lx, 60)
  out[[paste0(pre, "_max60")]] <- rmax(lx, 60) - lx                                   # distance below recent peak
  as.data.table(out)
}
F <- list()
# SXR block
F$sxr <- cbind(chan_feats(d$xrsb_flux, "sxr_b"), chan_feats(d$xrsa_flux, "sxr_a"))
F$sxr[, sxr_hardness := log10(d$xrsa_flux) - log10(d$xrsb_flux)]                      # spectral hardness / temperature proxy
F$sxr[, sxr_hardness_d15 := sxr_hardness - lag_(sxr_hardness, 15)]
# HXR block (orbit-modulated background removed by a 3-hour running median)
F$hxr <- cbind(chan_feats(d$hxr_lo, "hxr_lo", bg_win = 180), chan_feats(d$hxr_hi, "hxr_hi", bg_win = 180))
F$hxr[, hxr_ratio := log10(d$hxr_hi) - log10(d$hxr_lo)]
# EUV block (26, 28, 30 nm lines)
F$euv <- cbind(chan_feats(d$irr_256, "euv_256"), chan_feats(d$irr_284, "euv_284"), chan_feats(d$irr_304, "euv_304"))
# Neupert block: HXR should track d(SXR)/dt. Residual and rolling correlation, available only when SXR and HXR are both used.
dsxr <- d$xrsb_flux - lag_(d$xrsb_flux, 1); hxe <- F$hxr$hxr_lo_excess
nep <- data.table(nep_dsxr = dsxr / d$xrsb_flux)
z <- function(v, n = 180) (v - rm_(v, n)) / (rsd(v, n) + 1e-12)
a <- z(nep$nep_dsxr); b <- z(hxe)
nep[, nep_resid := b - a]                                                          # HXR deviation from what d(SXR)/dt predicts
nep[, nep_corr15 := frollapply(seq_along(a), 15, function(i) { aa <- a[i]; bb <- b[i]; ok <- is.finite(aa) & is.finite(bb); if (sum(ok) < 8) NA_real_ else cor(aa[ok], bb[ok]) }, align = "right")]
nep[, nep_resid_sd15 := rsd(nep_resid, 15)]
F$nep <- nep

X <- cbind(d[, .(time)], F$sxr, F$hxr, F$euv, F$nep)

# ---- labels: a NEW flare >= C starts within the next H minutes; rows inside an ongoing flare are excluded ----
ev <- ev[class_letter %in% c("C", "M", "X") & !is.na(start)]
setorder(ev, start)
tt <- as.numeric(X$time)
st <- as.numeric(ev$start); en <- as.numeric(ev$end)
nxt <- findInterval(tt, st) + 1                                       # index of next flare start strictly after t
nxt_start <- ifelse(nxt <= length(st), st[pmin(nxt, length(st))], Inf)
ongoing <- vapply(tt, function(t) any(st <= t & en >= t), logical(1))
X[, ongoing := ongoing]
X[, mins_to_next := (nxt_start - tt) / 60]
for (H in c(15, 30, 60)) X[[paste0("y", H)]] <- as.integer(X$mins_to_next > 0 & X$mins_to_next <= H)
# 'strong' target: M+ only
evM <- ev[class_letter %in% c("M", "X")]; stM <- as.numeric(evM$start); nM <- findInterval(tt, stM) + 1
X[, mins_to_next_M := ifelse(nM <= length(stM), stM[pmin(nM, length(stM))] - tt, Inf) / 60]
for (H in c(30, 60)) X[[paste0("yM", H)]] <- as.integer(X$mins_to_next_M > 0 & X$mins_to_next_M <= H)

# all rows must have every band so subset comparisons use identical samples
feat_cols <- setdiff(names(X), c("time", "ongoing", "mins_to_next", "mins_to_next_M", grep("^yM?[0-9]+$", names(X), value = TRUE)))
X[, complete := !ongoing & complete.cases(.SD), .SDcols = feat_cols]
fwrite(X, file.path(proc, "features.csv"))
message("rows: ", nrow(X), "; complete & not-ongoing: ", sum(X$complete))
for (H in c(15, 30, 60)) message("positives y", H, ": ", sum(X[complete == TRUE][[paste0("y", H)]]))
message("distinct flares (>=C): ", nrow(ev), "; M+: ", nrow(evM))
