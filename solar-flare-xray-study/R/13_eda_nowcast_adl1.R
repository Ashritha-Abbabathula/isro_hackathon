# 13_eda_nowcast_adl1.R -- exploratory figures, Neupert check and band-by-band nowcasting detection latency (Aditya-L1 data)
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
proc <- Sys.getenv("SOLAR_PROC", file.path(root, "data", "processed_adl1")); out <- Sys.getenv("SOLAR_OUT", file.path(root, "outputs_adl1"))
for (p in c("tables", "figures")) dir.create(file.path(out, p), FALSE, TRUE)
to_utc <- function(v) if (inherits(v, "POSIXct")) { attr(v, "tzone") <- "UTC"; v } else as.POSIXct(as.character(v), tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ")
rd <- function(f) { x <- fread(file.path(proc, f)); for (cn in intersect(c("time", "start", "peak", "end"), names(x))) x[[cn]] <- to_utc(x[[cn]]); x }
fl <- function(x) as.POSIXct(floor(as.numeric(x) / 60) * 60, origin = "1970-01-01", tz = "UTC")
sx <- rd("solexs_1min.csv"); hx <- rd("hel1os_1min.csv"); gs <- rd("goes_1min.csv"); ev <- rd("flare_events.csv"); sepdt <- rd("sep_periods.csv")
sx[, time := fl(time)]; hx[, time := fl(time)]; gs[, time := fl(time)]
d <- merge(merge(sx, hx, by = "time"), gs[, .(time, xrsb_flux, irr_304, irr_256)], by = "time", all.x = TRUE); setorder(d, time)
d <- merge(data.table(time = seq(min(d$time), max(d$time), by = "min")), d, by = "time", all.x = TRUE)
d[, `:=`(hxr_lo = czt_20_40, hxr_hi = czt_60_80 + czt_80_150)]
d[, sep := FALSE]; for (i in seq_len(nrow(sepdt))) { s0 <- sepdt$start[i]; s1 <- sepdt$end[i]; d[time >= s0 & time <= s1, sep := TRUE] }
ev <- ev[!is.na(start) & !is.na(peak) & !is.na(end) & class_letter %in% c("C", "M", "X")]
th <- theme_minimal(base_size = 9) + theme(panel.grid.minor = element_blank(), legend.position = "bottom")

# ---- 1. multi-band light curves of the strongest flare in the window ----
cand <- ev[peak >= min(d$time) & peak <= max(d$time) & !is.na(peak_flux)][order(-peak_flux)][1]
w <- d[time >= cand$start - 90 * 60 & time <= cand$end + 60 * 60]
long <- melt(w[, .(time, `SoLEXS 1-22 keV (cts/s)` = sxr_cps, `HEL1OS 20-40 keV (cts/s)` = hxr_lo, `HEL1OS 60-150 keV (cts/s)` = hxr_hi, `GOES EUV 30.4 nm` = irr_304)], id.vars = "time")
long[, band := fifelse(grepl("^SoLEXS", variable), "SXR", fifelse(grepl("^HEL1OS", variable), "HXR", "EUV"))]
g1 <- ggplot(long, aes(time, value, colour = variable)) + geom_line(linewidth = .35, na.rm = TRUE) + facet_grid(band ~ ., scales = "free_y") +
  geom_vline(xintercept = cand$peak, linewidth = .25) + geom_vline(xintercept = c(cand$start, cand$end), linetype = "dashed", linewidth = .25) + scale_y_log10() +
  labs(x = "UTC", y = NULL, colour = NULL) + th + guides(colour = guide_legend(nrow = 2))
ggsave(file.path(out, "figures", "fig_example_flare.png"), g1, width = 3.5, height = 4.4, dpi = 300)
writeLines(sprintf("%s|%s|%s", cand$class, format(cand$peak, "%d %B %Y"), format(cand$peak, "%H:%M")), file.path(out, "tables", "example_flare.txt"))

# ---- 2. Neupert check on SoLEXS vs HEL1OS during the rise phase (M/X flares outside SEP periods) ----
MX <- ev[class_letter %in% c("M", "X")]
per <- rbindlist(lapply(seq_len(nrow(MX)), function(i) { e <- MX[i]; x <- d[time >= e$start - 10 * 60 & time <= e$peak + 5 * 60]
  if (nrow(x) < 8 || any(x$sep, na.rm = TRUE)) return(NULL)
  x[, dsxr := c(NA, diff(sxr_cps))]; x[, hx := hxr_lo - median(hxr_lo[1:min(10, .N)], na.rm = TRUE)]
  ok <- x[complete.cases(dsxr, hx)]; if (nrow(ok) < 8) return(NULL)
  data.table(flare = e$flare_id, class = e$class, r = suppressWarnings(cor(ok$dsxr, ok$hx, method = "spearman")), n = nrow(ok)) }))
per <- per[is.finite(r)]; fwrite(per, file.path(out, "tables", "neupert_per_flare_spearman.csv"))
message("Neupert (M/X, SoLEXS vs HEL1OS 20-40 keV): median r = ", round(median(per$r), 2), " over ", nrow(per), " flares; r>0 in ", round(100 * mean(per$r > 0)), "%")
g2 <- ggplot(per, aes(r)) + geom_histogram(bins = 16, fill = "#c0392b", colour = "white") + geom_vline(xintercept = 0, linewidth = .3) +
  labs(x = "Spearman r (HEL1OS 20-40 keV vs. dSoLEXS/dt)", y = "M/X flares") + th
ggsave(file.path(out, "figures", "fig_neupert.png"), g2, width = 3.4, height = 2.3, dpi = 300)

# ---- 3. nowcasting detection latency per band relative to the GOES flare start ----
rz <- function(x, bgw = 180, lagm = 5) { lx <- log10(pmax(x, 1e-3))
  med <- shift(frollapply(lx, bgw, median, na.rm = TRUE, align = "right"), lagm)
  mad <- shift(frollapply(lx, bgw, function(v) median(abs(v - median(v, na.rm = TRUE)), na.rm = TRUE), align = "right"), lagm)
  (lx - med) / pmax(1.4826 * mad, 0.02) }
Z <- data.table(time = d$time, `SoLEXS (SXR)` = rz(d$sxr_cps), `HEL1OS 20-40 keV` = rz(d$hxr_lo), `HEL1OS 60-150 keV` = rz(d$hxr_hi), `GOES EUV 30.4 nm` = rz(d$irr_304))
Z0 <- 5; NCONS <- 2
detect <- function(z, t, s, pk) { i <- which(t >= s - 20 * 60 & t <= pk); if (length(i) < NCONS) return(NA_real_)
  zz <- z[i]; ok <- which(frollsum(as.numeric(zz > Z0), NCONS, align = "right") >= NCONS)
  if (!length(ok)) return(NA_real_); (as.numeric(t[i[ok[1]]]) - (NCONS - 1) * 60 - as.numeric(s)) / 60 }
tt <- as.numeric(Z$time); lat <- list()
ev_in <- ev[start >= min(d$time) + 4 * 3600 & peak <= max(d$time)]
for (j in seq_len(nrow(ev_in))) { e <- ev_in[j]; ix <- which(tt >= as.numeric(e$start) - 20 * 60 & tt <= as.numeric(e$peak)); if (length(ix) < 5) next
  if (any(d$sep[ix], na.rm = TRUE)) next
  for (b in setdiff(names(Z), "time")) lat[[length(lat) + 1]] <- data.table(flare = e$flare_id, cls = e$class_letter, class = e$class, band = b, latency = detect(Z[[b]][ix], Z$time[ix], e$start, e$peak), z_avail = mean(is.finite(Z[[b]][ix]))) }
L <- rbindlist(lat)[z_avail > 0.8]; fwrite(L, file.path(out, "tables", "nowcast_latency_per_flare.csv"))
Ls <- L[, .(n = .N, detected = mean(is.finite(latency)), median_latency_min = median(latency, na.rm = TRUE), q25 = quantile(latency, .25, na.rm = TRUE), q75 = quantile(latency, .75, na.rm = TRUE)), by = .(class_group = ifelse(cls == "C", "C", "M/X"), band)]
fwrite(Ls, file.path(out, "tables", "nowcast_latency_summary.csv")); print(Ls)
# paired comparison: on flares where both HEL1OS 20-40 keV and SoLEXS detect, who is earlier?
PW <- dcast(L[band %in% c("SoLEXS (SXR)", "HEL1OS 20-40 keV")], flare + cls ~ band, value.var = "latency")
PW <- PW[is.finite(`SoLEXS (SXR)`) & is.finite(`HEL1OS 20-40 keV`)][, diff := `HEL1OS 20-40 keV` - `SoLEXS (SXR)`]
PWs <- PW[, .(n_both = .N, median_diff_min = median(diff), hel_earlier = mean(diff < 0), same_minute = mean(diff == 0), hel_later = mean(diff > 0)), by = .(class_group = ifelse(cls == "C", "C", "M/X"))]
fwrite(PWs, file.path(out, "tables", "nowcast_paired.csv")); print(PWs)
g3 <- ggplot(L[is.finite(latency)][, grp := ifelse(cls == "C", "C-class", "M/X-class")], aes(band, latency, fill = band)) + geom_boxplot(outlier.size = .4, linewidth = .25) + facet_wrap(~grp) +
  geom_hline(yintercept = 0, linewidth = .3, linetype = "dashed") + coord_cartesian(ylim = c(-20, 15)) + scale_x_discrete(labels = c("GOES EUV 30.4 nm" = "EUV\n30.4", "HEL1OS 20-40 keV" = "HEL1OS\n20-40", "HEL1OS 60-150 keV" = "HEL1OS\n60-150", "SoLEXS (SXR)" = "SoLEXS")) + labs(x = NULL, y = "Detection time - GOES start (min)") + th + theme(axis.text.x = element_text(size = 6), legend.position = "none")
ggsave(file.path(out, "figures", "fig_nowcast_latency.png"), g3, width = 3.5, height = 2.8, dpi = 300)

# ---- dataset summary numbers ----
ev_w <- ev[peak >= min(d$time) & peak <= max(d$time)]
fwrite(ev_w[, .N, by = class_letter][order(class_letter)], file.path(out, "tables", "flare_counts_window.csv"))
message("EDA done; flares in analysis window: ", nrow(ev_w))
