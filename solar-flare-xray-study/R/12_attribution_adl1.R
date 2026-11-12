# 12_attribution_adl1.R -- Shapley contribution of each band beyond flare history, with 6-hour block-bootstrap CIs.
# v(S) = pooled out-of-sample TSS of the model that uses flare history + the bands in S;  v(empty) = history-only model ("NONE").
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
set.seed(7)
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
out <- Sys.getenv("SOLAR_OUT", file.path(root, "outputs_adl1"))
P <- readRDS(file.path(out, "models", "preds.rds")); R <- fread(file.path(out, "tables", "results.csv"))
BL <- c("SXR", "HXR", "EUV")
nm <- function(S) if (length(S)) paste(S[order(match(S, BL))], collapse = "+") else "NONE"
all_sub <- c(list(character(0)), unlist(lapply(1:3, function(k) combn(BL, k, simplify = FALSE)), recursive = FALSE)); sub_names <- vapply(all_sub, nm, "")
tss <- function(y, p, thr) { pr <- p >= thr; sum(pr & y == 1) / max(1, sum(y == 1)) - sum(pr & y == 0) / max(1, sum(y == 0)) }
v_of <- function(tab, idx) vapply(sub_names, function(s) { d <- tab[[s]]; if (is.null(d)) NA_real_ else tss(d$y[idx], d$p[idx], d$thr[idx]) }, numeric(1))
shapley <- function(v) { names(v) <- sub_names
  vapply(BL, function(i) { oth <- setdiff(BL, i); tot <- 0
    for (k in 0:length(oth)) for (S in (if (k == 0) list(character(0)) else combn(oth, k, simplify = FALSE))) {
      w <- factorial(length(S)) * factorial(length(BL) - length(S) - 1) / factorial(length(BL)); tot <- tot + w * (v[[nm(c(S, i))]] - v[[nm(S)]]) }
    tot }, numeric(1)) }

rows <- list(); deltas <- list(); B <- 500
for (tg in unique(P$target)) for (md in unique(P$model)) {
  d <- P[target == tg & model == md & subset %in% sub_names]
  if (uniqueN(d$subset) < 8) next
  tab <- split(d[, .(time, y, p, thr)], d$subset)
  common <- Reduce(intersect, lapply(tab, function(x) as.numeric(x$time))); tab <- lapply(tab, function(x) { x <- x[as.numeric(time) %in% common]; setorder(x, time); x })
  n <- nrow(tab[[1]]); tnum <- as.numeric(tab[[1]]$time)
  blk <- floor((tnum - min(tnum)) / (6 * 3600)); bi <- split(seq_len(n), blk); ub <- names(bi)
  v0 <- v_of(tab, seq_len(n)); s0 <- shapley(v0)
  vb <- matrix(NA_real_, B, length(sub_names), dimnames = list(NULL, sub_names)); sb <- matrix(NA_real_, B, 3, dimnames = list(NULL, BL))
  for (b in seq_len(B)) { idx <- unlist(bi[sample(ub, length(ub), TRUE)], use.names = FALSE); vb[b, ] <- v_of(tab, idx); sb[b, ] <- shapley(vb[b, ]) }
  rows[[length(rows) + 1]] <- data.table(target = tg, model = md, band = BL, shapley_TSS = as.numeric(s0), lo = apply(sb, 2, quantile, .025, na.rm = TRUE), hi = apply(sb, 2, quantile, .975, na.rm = TRUE),
                                         TSS_history = v0[["NONE"]], TSS_full = v0[["SXR+HXR+EUV"]], gain = v0[["SXR+HXR+EUV"]] - v0[["NONE"]])
  for (s in setdiff(sub_names, "NONE")) { dd <- vb[, s] - vb[, "NONE"]
    deltas[[length(deltas) + 1]] <- data.table(target = tg, model = md, subset = s, TSS = v0[[s]], delta_vs_history = v0[[s]] - v0[["NONE"]], lo = quantile(dd, .025, na.rm = TRUE), hi = quantile(dd, .975, na.rm = TRUE)) }
}
SH <- rbindlist(rows); DL <- rbindlist(deltas)
fwrite(SH, file.path(out, "tables", "shapley_band_contribution.csv")); fwrite(DL, file.path(out, "tables", "delta_vs_history.csv"))
# consensus across the model families (mean of per-model Shapley values; spread = SD across models)
CS <- SH[, .(shapley_mean = mean(shapley_TSS), shapley_sd = sd(shapley_TSS), n_models = .N, n_models_ci_above0 = sum(lo > 0), n_models_ci_below0 = sum(hi < 0)), by = .(target, band)]
fwrite(CS, file.path(out, "tables", "shapley_consensus.csv"))

th <- theme_minimal(base_size = 10) + theme(panel.grid.minor = element_blank(), legend.position = "bottom")
pal <- c(SXR = "#1b6ca8", HXR = "#c0392b", EUV = "#2e8b57")
tg_lab <- c(y15 = "C+ 15 min", y30 = "C+ 30 min", y60 = "C+ 60 min", yM30 = "M+ 30 min", yM60 = "M+ 60 min")
CS[, tf := factor(tg_lab[target], levels = tg_lab)]; SH[, tf := factor(tg_lab[target], levels = tg_lab)]
g1 <- ggplot(CS, aes(band, shapley_mean, fill = band)) + geom_col(width = .6, alpha = .85) + geom_errorbar(aes(ymin = shapley_mean - shapley_sd, ymax = shapley_mean + shapley_sd), width = .15) +
  geom_point(data = SH, aes(band, shapley_TSS, shape = model), size = 1.1, inherit.aes = FALSE) +
  facet_wrap(~tf, nrow = 1) + scale_fill_manual(values = pal, guide = "none") + geom_hline(yintercept = 0, linewidth = .3) +
  labs(x = NULL, y = "Shapley contribution to TSS", shape = "Model family") + th
ggsave(file.path(out, "figures", "fig_shapley_bands.png"), g1, width = 7.2, height = 2.6, dpi = 300)
R[, subset := factor(subset, levels = c("level_only", "NONE", "SXR", "HXR", "EUV", "SXR+HXR", "SXR+EUV", "HXR+EUV", "SXR+HXR+EUV"))]
g2 <- ggplot(R[target %in% c("y15", "y30", "y60") & subset != "level_only"], aes(subset, TSS, fill = model)) + geom_col(position = "dodge", width = .75) + facet_wrap(~target, ncol = 1) +
  geom_hline(yintercept = 0, linewidth = .3) + labs(x = NULL, y = "Pooled out-of-sample TSS", fill = NULL) + th + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(out, "figures", "fig_ablation_tss.png"), g2, width = 3.5, height = 5.4, dpi = 300)
message("attribution done: ", nrow(SH), " shapley rows, ", nrow(DL), " delta rows")
