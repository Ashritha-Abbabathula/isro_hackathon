# 05_attribution.R -- how much does each band contribute? Shapley values over band subsets + block-bootstrap CIs
# Input: outputs/models/preds.rds, outputs/tables/results.csv ; Output: outputs/tables/*.csv, outputs/figures/*.png
suppressPackageStartupMessages({ library(data.table); library(ggplot2) })
set.seed(7)
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
out <- file.path(root, "outputs"); P <- readRDS(file.path(out, "models", "preds.rds")); R <- fread(file.path(out, "tables", "results.csv"))
BL <- c("SXR", "HXR", "EUV")
subs_of <- function(S) if (length(S)) paste(S, collapse = "+") else "NONE"
all_sub <- c(list(character(0)), unlist(lapply(1:3, function(k) combn(BL, k, simplify = FALSE)), recursive = FALSE))

tss <- function(y, p, thr) { pr <- p >= thr; pod <- sum(pr & y == 1) / max(1, sum(y == 1)); pofd <- sum(pr & y == 0) / max(1, sum(y == 0)); pod - pofd }

# characteristic function v(S): test TSS of the model using exactly the bands in S (v(empty) = 0, climatology)
v_of <- function(tab, idx) {   # tab: list name -> data.table(y,p,thr); idx: integer row index (possibly resampled)
  sapply(all_sub, function(S) if (!length(S)) 0 else { d <- tab[[subs_of(S)]]; if (is.null(d)) NA_real_ else tss(d$y[idx], d$p[idx], d$thr[idx]) }) }
shapley <- function(v) {  # exact Shapley for 3 players; v ordered as all_sub
  names(v) <- sapply(all_sub, subs_of); get <- function(S) v[[subs_of(S)]]
  sapply(BL, function(i) { oth <- setdiff(BL, i); tot <- 0
    for (k in 0:length(oth)) for (S in (if (k == 0) list(character(0)) else combn(oth, k, simplify = FALSE))) {
      Sc <- if (length(S)) S[order(match(S, BL))] else S; Si <- c(Sc, i); Si <- Si[order(match(Si, BL))]
      w <- factorial(length(S)) * factorial(length(BL) - length(S) - 1) / factorial(length(BL)); tot <- tot + w * (get(Si) - get(Sc)) }
    tot }) }

rows <- list(); ablate <- list()
for (tg in unique(P$target)) for (md in unique(P$model)) {
  d <- P[target == tg & model == md & subset %in% sapply(all_sub[-1], subs_of)]
  if (uniqueN(d$subset) < 7) next
  tab <- split(d[, .(y, p, thr, time)], d$subset); n <- nrow(tab[[1]])
  blk <- as.integer(floor((as.numeric(tab[[1]]$time) - min(as.numeric(tab[[1]]$time))) / (6 * 3600)))   # 6-hour blocks preserve autocorrelation
  ub <- unique(blk); bi <- split(seq_len(n), blk)
  v0 <- v_of(tab, seq_len(n)); s0 <- shapley(v0)
  B <- 400; sb <- matrix(NA_real_, B, 3, dimnames = list(NULL, BL)); vb <- matrix(NA_real_, B, length(all_sub))
  for (b in seq_len(B)) { idx <- unlist(bi[as.character(sample(ub, length(ub), TRUE))], use.names = FALSE); vb[b, ] <- v_of(tab, idx); sb[b, ] <- shapley(vb[b, ]) }
  names(v0) <- sapply(all_sub, subs_of)
  rows[[length(rows) + 1]] <- data.table(target = tg, model = md, band = BL, shapley_TSS = as.numeric(s0),
    lo = apply(sb, 2, quantile, .025, na.rm = TRUE), hi = apply(sb, 2, quantile, .975, na.rm = TRUE), TSS_full = v0[["SXR+HXR+EUV"]],
    share = as.numeric(s0) / sum(s0))
  for (b in BL) { oth <- setdiff(BL, b); full <- v0[["SXR+HXR+EUV"]]; drop <- v0[[subs_of(oth)]]; alone <- v0[[b]]
    j <- which(sapply(all_sub, subs_of) == subs_of(oth)); jb <- which(sapply(all_sub, subs_of) == "SXR+HXR+EUV")
    dd <- vb[, jb] - vb[, j]
    ablate[[length(ablate) + 1]] <- data.table(target = tg, model = md, band = b, TSS_alone = alone, TSS_drop_this_band = full - drop, lo = quantile(dd, .025, na.rm = TRUE), hi = quantile(dd, .975, na.rm = TRUE)) }
}
SH <- rbindlist(rows); AB <- rbindlist(ablate)
fwrite(SH, file.path(out, "tables", "shapley_band_contribution.csv")); fwrite(AB, file.path(out, "tables", "ablation_drop_one.csv"))

# ---- figures ----
th <- theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank(), legend.position = "bottom")
pal <- c(SXR = "#1b6ca8", HXR = "#c0392b", EUV = "#2e8b57")
if (nrow(SH)) {
  g <- ggplot(SH[model == "xgboost"], aes(band, shapley_TSS, fill = band)) + geom_col(width = .6) + geom_errorbar(aes(ymin = lo, ymax = hi), width = .15) +
    facet_wrap(~target, nrow = 1) + scale_fill_manual(values = pal) + geom_hline(yintercept = 0, linewidth = .3) +
    labs(title = "Shapley contribution of each band to test TSS (XGBoost)", subtitle = "Bars: exact Shapley value over the 7 band subsets; whiskers: 95% block-bootstrap CI", x = NULL, y = "Contribution to TSS") + th
  ggsave(file.path(out, "figures", "fig_shapley_bands.png"), g, width = 10, height = 3.8, dpi = 200)
}
R[, subset := factor(subset, levels = c("baseline_level", "SXR", "HXR", "EUV", "SXR+HXR", "SXR+EUV", "HXR+EUV", "SXR+HXR+EUV"))]
g2 <- ggplot(R[target %in% c("y15", "y30", "y60") & subset != "baseline_level"], aes(subset, TSS, fill = model)) + geom_col(position = "dodge", width = .75) + facet_wrap(~target, ncol = 1) +
  geom_hline(yintercept = 0, linewidth = .3) + labs(title = "Forecast skill (TSS) by band subset and model", x = NULL, y = "TSS (test)") + th + theme(axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(out, "figures", "fig_ablation_tss.png"), g2, width = 8, height = 8, dpi = 200)
g3 <- ggplot(R[subset == "SXR+HXR+EUV" & model == "xgboost" & !is.na(lead_median)], aes(horizon, lead_median, colour = target)) + geom_point(size = 3) + labs(title = "Median lead time of alerts before flare start", x = "Forecast horizon (min)", y = "Median lead (min)") + th
ggsave(file.path(out, "figures", "fig_lead_time.png"), g3, width = 6, height = 4, dpi = 200)
message("attribution done: ", nrow(SH), " shapley rows")
