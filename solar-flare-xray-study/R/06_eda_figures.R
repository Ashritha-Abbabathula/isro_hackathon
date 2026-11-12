# 06_eda_figures.R -- exploratory figures: multi-band light curves of the strongest flare, Neupert check, dataset summary
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
proc <- file.path(root, "data", "processed"); fig <- file.path(root, "outputs", "figures"); tab <- file.path(root, "outputs", "tables")
to_utc <- function(v) if (inherits(v, "POSIXct")) { attr(v, "tzone") <- "UTC"; v } else as.POSIXct(as.character(v), tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ")
g <- fread(file.path(proc, "goes_1min.csv")); g[, time := to_utc(time)]
h <- fread(file.path(proc, "hxr_1min.csv")); h[, time := to_utc(time)]
ev <- fread(file.path(proc, "flare_events.csv")); for (cn in c("start", "peak", "end")) ev[[cn]] <- to_utc(ev[[cn]]); ev <- ev[!is.na(start) & class_letter %in% c("C", "M", "X")]
d <- merge(g, h[, .(time, hxr_lo, hxr_hi)], by = "time"); th <- theme_minimal(base_size = 10) + theme(panel.grid.minor = element_blank())

# dataset summary table
sm <- ev[time_in <- peak >= min(d$time) & peak <= max(d$time), .N, by = class_letter][order(class_letter)]
fwrite(sm, file.path(tab, "flare_counts_hxr_window.csv"))
fwrite(ev[, .N, by = class_letter][order(class_letter)], file.path(tab, "flare_counts_goes_window.csv"))

# strongest flare inside the Fermi window: stacked panels
cand <- ev[peak >= min(d$time) & peak <= max(d$time)][order(-peak_flux)][1]
w <- d[time >= cand$start - 90 * 60 & time <= cand$end + 60 * 60]
long <- melt(w[, .(time, `SXR 1-8 A (W/m2)` = xrsb_flux, `SXR 0.5-4 A (W/m2)` = xrsa_flux, `HXR 26-50 keV (cts/s)` = hxr_lo, `HXR 50-101 keV (cts/s)` = hxr_hi, `EUV 30.4 nm (W/m2)` = irr_304, `EUV 25.6 nm (W/m2)` = irr_256)], id.vars = "time")
long[, band := fifelse(grepl("^SXR", variable), "SXR", fifelse(grepl("^HXR", variable), "HXR", "EUV"))]
g1 <- ggplot(long, aes(time, value, colour = variable)) + geom_line(linewidth = .4, na.rm = TRUE) + facet_grid(band ~ ., scales = "free_y") +
  geom_vline(xintercept = cand$peak, linetype = "solid", linewidth = .3) + geom_vline(xintercept = c(cand$start, cand$end), linetype = "dashed", linewidth = .3) +
  scale_y_log10() + labs(title = sprintf("Multi-band light curves of the strongest flare in the study window (%s, peak %s UTC)", cand$class, format(cand$peak, "%d %b %H:%M")),
                         subtitle = "Dashed: GOES start/end; solid: GOES peak", x = NULL, y = NULL, colour = NULL) + th + theme(legend.position = "bottom")
ggsave(file.path(fig, "fig_example_flare.png"), g1, width = 9, height = 7, dpi = 200)

# Neupert check: HXR (26-50 keV) vs d(SXR)/dt for the flares in the window, using the rise phase (start -> peak)
rise <- rbindlist(lapply(seq_len(nrow(ev[peak >= min(d$time) & peak <= max(d$time) & class_letter %in% c("M", "X")])), function(i) {
  e <- ev[peak >= min(d$time) & peak <= max(d$time) & class_letter %in% c("M", "X")][i]
  x <- d[time >= e$start - 10 * 60 & time <= e$peak + 5 * 60]; if (nrow(x) < 8) return(NULL)
  x[, `:=`(dsxr = c(NA, diff(xrsb_flux)), flare = e$flare_id)]
  x[, hx := hxr_lo - median(hxr_lo[1:min(10, .N)], na.rm = TRUE)]
  x[complete.cases(dsxr, hx)] }))
if (nrow(rise) > 20) {
  rise[, r := cor(dsxr, hx, method = "spearman"), by = flare]
  per <- unique(rise[, .(flare, r, n = .N), by = flare][, .(flare, r, n)])
  fwrite(per, file.path(tab, "neupert_per_flare_spearman.csv"))
  message("Neupert (M/X flares): median per-flare Spearman r = ", round(median(per$r, na.rm = TRUE), 2), " over ", nrow(per), " flares")
  g2 <- ggplot(per, aes(r)) + geom_histogram(bins = 15, fill = "#c0392b", colour = "white") + geom_vline(xintercept = 0, linewidth = .3) +
    labs(title = "Neupert effect: per-flare rank correlation of HXR (26-50 keV) with dSXR/dt", subtitle = "M/X flares, window start-10 min to peak+5 min", x = "Spearman r", y = "Flares") + th
  ggsave(file.path(fig, "fig_neupert.png"), g2, width = 6, height = 4, dpi = 200)
}
message("EDA figures done")
