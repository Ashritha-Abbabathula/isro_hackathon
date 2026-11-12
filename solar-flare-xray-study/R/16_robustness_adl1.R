# 16_robustness_adl1.R -- per-block stability, paired block analysis (captures refit variability), alert-episode metrics, forward-vs-LOO comparison
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
proc <- Sys.getenv("SOLAR_PROC", file.path(root, "data", "processed_adl1")); out <- Sys.getenv("SOLAR_OUT", file.path(root, "outputs_adl1")); outF <- Sys.getenv("SOLAR_OUT_FWD", file.path(root, "outputs_adl1_forward")); tb <- file.path(out, "tables")
to_utc <- function(v) if (inherits(v, "POSIXct")) { attr(v, "tzone") <- "UTC"; v } else as.POSIXct(as.character(v), tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ")
P <- readRDS(file.path(out, "models", "preds.rds")); P[, time := to_utc(time)]
ev <- fread(file.path(proc, "flare_events.csv")); ev[, start := to_utc(start)]; ev <- ev[!is.na(start) & class_letter %in% c("C", "M", "X")]
starts <- list(C = sort(as.numeric(ev$start)), M = sort(as.numeric(ev[class_letter %in% c("M", "X")]$start)))
HZ <- c(y15 = 15, y30 = 30, y60 = 60, yM30 = 30, yM60 = 60)
tss <- function(y, p, thr) { pr <- p >= thr; sum(pr & y == 1) / max(1, sum(y == 1)) - sum(pr & y == 0) / max(1, sum(y == 0)) }

# ---- 1. per-block TSS ----
PB <- P[, .(TSS = tss(y, p, thr[1]), n = .N, pos = sum(y)), by = .(target, model, subset, fold)]
fwrite(PB, file.path(tb, "per_block_tss.csv"))

# ---- 2. paired block analysis vs the history-only model (each block = an independent refit + disjoint test period) ----
H <- PB[subset == "NONE", .(target, model, fold, TSS_hist = TSS)]
PBd <- merge(PB[subset != "NONE" & subset != "level_only"], H, by = c("target", "model", "fold"))[, d := TSS - TSS_hist]
PA <- PBd[, .(n_blocks = .N, mean_delta = mean(d), se = sd(d) / sqrt(.N), blocks_positive = sum(d > 0)), by = .(target, model, subset)]
PA[, `:=`(lo = mean_delta - qt(.975, n_blocks - 1) * se, hi = mean_delta + qt(.975, n_blocks - 1) * se)]
fwrite(PA, file.path(tb, "paired_block_vs_history.csv"))
PA2 <- PA[, .(mean_delta = mean(mean_delta), ci_above0 = sum(lo > 0), ci_below0 = sum(hi < 0), n = .N, blocks_pos = mean(blocks_positive / n_blocks)), by = subset]
fwrite(PA2, file.path(tb, "paired_block_summary_by_subset.csv")); print(PA2)

# level-only baseline vs full model per block (M+ stability check)
LV <- dcast(PB[subset %in% c("level_only", "SXR+HXR+EUV", "NONE") & model %in% c("glmnet", "xgboost") & target %in% c("yM30", "yM60", "y30")], target + model + fold ~ subset, value.var = "TSS")
fwrite(LV, file.path(tb, "level_only_vs_full_per_block.csv"))

# ---- 3. alert-episode metrics ----
alert_metrics <- function(d, st, H) {
  d <- d[order(fold, time)]; al <- d[p >= thr]; ndays <- nrow(d) / 1440
  rng <- d[, .(a = min(as.numeric(time)), b = max(as.numeric(time))), by = fold]; fs <- st[vapply(st, function(s) any(s >= rng$a & s <= rng$b), logical(1))]
  if (!nrow(al)) return(data.table(n_episodes = 0L, episodes_per_day = 0, episode_precision = NA_real_, flare_detect_rate = 0, lead_median_min = NA_real_, alert_minute_share = 0, n_flares = length(fs)))
  tt <- as.numeric(al$time); new <- c(TRUE, diff(tt) > 300 | diff(al$fold) != 0); id <- cumsum(new)
  e0 <- tapply(tt, id, min); e1 <- tapply(tt, id, max)
  ver <- findInterval(e1 + H * 60, st) > findInterval(e0, st)                                  # >=1 flare starts after the episode begins and within H min of its end
  det <- vapply(fs, function(s) { i <- which(tt >= s - H * 60 & tt < s); if (length(i)) (s - tt[min(i)]) / 60 else NA_real_ }, numeric(1))
  data.table(n_episodes = length(e0), episodes_per_day = length(e0) / ndays, episode_precision = mean(ver), flare_detect_rate = mean(is.finite(det)), lead_median_min = median(det, na.rm = TRUE),
             alert_minute_share = nrow(al) / nrow(d), n_flares = length(fs)) }
AL <- P[, { st <- if (grepl("^yM", target[1])) starts$M else starts$C; alert_metrics(.SD, st, HZ[[target[1]]]) }, by = .(target, model, subset), .SDcols = c("fold", "time", "p", "thr")]
fwrite(AL, file.path(tb, "alert_episode_metrics.csv"))
print(AL[subset %in% c("NONE", "SXR+HXR+EUV") & model == "xgboost"][order(target, subset)])

# ---- 4. forward-in-time vs leave-one-block-out on the SAME test blocks (3..8) ----
if (file.exists(file.path(outF, "models", "preds.rds"))) {
  PF <- readRDS(file.path(outF, "models", "preds.rds")); PF[, time := to_utc(time)]
  a <- P[fold >= 3 & model %in% unique(PF$model), .(TSS_loo = tss(y, p, thr[1])), by = .(target, model, subset, fold)]
  b <- PF[, .(TSS_forward = tss(y, p, thr[1])), by = .(target, model, subset, fold)]
  FW <- merge(a, b, by = c("target", "model", "subset", "fold"))
  # pooled over blocks 3..8
  pa <- P[fold >= 3 & model %in% unique(PF$model), .(TSS_loo = tss(y, p, thr), n = .N), by = .(target, model, subset)]
  pa <- P[fold >= 3 & model %in% unique(PF$model)][, .(TSS_loo = { pr <- p >= thr; sum(pr & y == 1) / max(1, sum(y == 1)) - sum(pr & y == 0) / max(1, sum(y == 0)) }), by = .(target, model, subset)]
  pb <- PF[, .(TSS_forward = { pr <- p >= thr; sum(pr & y == 1) / max(1, sum(y == 1)) - sum(pr & y == 0) / max(1, sum(y == 0)) }), by = .(target, model, subset)]
  FWp <- merge(pa, pb, by = c("target", "model", "subset")); FWp[, diff := TSS_forward - TSS_loo]
  fwrite(FWp, file.path(tb, "forward_vs_loo_pooled.csv")); fwrite(FW, file.path(tb, "forward_vs_loo_per_block.csv"))
  print(FWp[, .(mean_loo = mean(TSS_loo), mean_forward = mean(TSS_forward), mean_diff = mean(diff)), by = subset][order(subset)])
  # forward-chaining deltas vs history (same paired block logic)
  HF <- PF[subset == "NONE", .(TSS_h = tss(y, p, thr[1])), by = .(target, model, fold)]
  DF <- merge(PF[subset %in% c("SXR", "HXR", "EUV", "SXR+HXR", "SXR+EUV", "HXR+EUV", "SXR+HXR+EUV"), .(TSS = tss(y, p, thr[1])), by = .(target, model, subset, fold)], HF, by = c("target", "model", "fold"))[, d := TSS - TSS_h]
  FA <- DF[, .(n_blocks = .N, mean_delta = mean(d), se = sd(d) / sqrt(.N), blocks_positive = sum(d > 0)), by = .(target, model, subset)]; FA[, `:=`(lo = mean_delta - qt(.975, n_blocks - 1) * se, hi = mean_delta + qt(.975, n_blocks - 1) * se)]
  fwrite(FA, file.path(tb, "forward_paired_block_vs_history.csv"))
}
# ---- figure: per-block TSS, full model vs history vs level-only ----
g <- ggplot(PB[subset %in% c("NONE", "SXR+HXR+EUV", "level_only") & model == "glmnet" & target %in% c("y30", "yM60")], aes(fold, TSS, colour = subset)) + geom_line() + geom_point() + facet_wrap(~target) +
  geom_hline(yintercept = 0, linewidth = .3) + labs(x = "Test block (chronological)", y = "TSS in block", colour = NULL) + theme_minimal(base_size = 9) + theme(legend.position = "bottom")
ggsave(file.path(out, "figures", "fig_per_block_tss.png"), g, width = 3.5, height = 2.6, dpi = 300)
message("robustness done")
