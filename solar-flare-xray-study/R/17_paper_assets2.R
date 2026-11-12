# 17_paper_assets2.R -- extra tables/macros for the revised paper: forward-in-time results, paired blocks, alerts, causality
suppressPackageStartupMessages(library(data.table))
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
pap <- Sys.getenv("PAPER_DIR", file.path(root, "paper", "ieee")); OM <- Sys.getenv("OUT_MAIN", file.path(root, "outputs_adl1")); OF <- Sys.getenv("OUT_FWD", file.path(root, "outputs_adl1_forward")); tb <- file.path(OM, "tables"); tbF <- file.path(OF, "tables")
f2 <- function(x, d = 2) ifelse(is.na(x), "--", formatC(x, format = "f", digits = d))
digitless <- function(n) gsub("15", "Fifteen", gsub("30", "Thirty", gsub("60", "Sixty", n)))
m <- c(); add <- function(name, val) m <<- c(m, sprintf("\\newcommand{\\%s}{%s}", digitless(name), val))
tg_lab <- c(y15 = "C+ 15", y30 = "C+ 30", y60 = "C+ 60", yM30 = "M+ 30", yM60 = "M+ 60"); TG <- names(tg_lab)
sub_order <- c("NONE", "SXR", "HXR", "EUV", "SXR+HXR", "SXR+EUV", "HXR+EUV", "SXR+HXR+EUV")
sub_tex <- c(NONE = Sys.getenv("BASELINE_LABEL", "History only"), SXR = "+SXR", HXR = "+HXR", EUV = "+EUV", "SXR+HXR" = "+SXR+HXR", "SXR+EUV" = "+SXR+EUV", "HXR+EUV" = "+HXR+EUV", "SXR+HXR+EUV" = "+SXR+HXR+EUV")

RF <- fread(file.path(tbF, "results.csv")); FW <- fread(file.path(tb, "forward_vs_loo_pooled.csv")); PS <- fread(file.path(tb, "paired_block_summary_by_subset.csv"))
FA <- fread(file.path(tb, "forward_paired_block_vs_history.csv")); AL <- fread(file.path(tb, "alert_episode_metrics.csv")); CT <- fread(file.path(tb, "causality_test.csv"))
CS <- fread(file.path(tb, "shapley_consensus.csv")); CSF <- fread(file.path(tbF, "shapley_consensus.csv")); SHF <- fread(file.path(tbF, "shapley_band_contribution.csv"))

# ---- macros ----
add("FwdDiffMean", f2(abs(mean(FW$diff)))); add("FwdShareLower", round(100 * mean(FW$diff < 0)))
if ("HISTONLY" %in% FW$subset) { add("MagGainLoo", f2(mean(FW[subset == "NONE"]$TSS_loo) - mean(FW[subset == "HISTONLY"]$TSS_loo), 3)); add("MagGainFwd", f2(mean(FW[subset == "NONE"]$TSS_forward) - mean(FW[subset == "HISTONLY"]$TSS_forward), 3)) }
add("LevelLoo", f2(mean(FW[subset == "level_only"]$TSS_loo))); add("LevelFwd", f2(mean(FW[subset == "level_only"]$TSS_forward)))
add("HistLooBlocks", f2(mean(FW[subset == "NONE"]$TSS_loo))); add("HistFwd", f2(mean(FW[subset == "NONE"]$TSS_forward)))
FS <- FA[, .(mean_delta = mean(mean_delta), ci_above0 = sum(lo > 0), ci_below0 = sum(hi < 0), n = .N), by = subset]
for (s in FS$subset) add(paste0("FwdGain", gsub("+", "", s, fixed = TRUE)), f2(FS[subset == s]$mean_delta, 3))
for (s in PS$subset) add(paste0("PairGain", gsub("+", "", s, fixed = TRUE)), f2(PS[subset == s]$mean_delta, 3))
PSb <- PS[subset %in% sub_order[-1]]; add("PairSigTotal", sum(PSb$ci_above0)); add("PairSigAll", sum(PSb$n)); add("FwdSigTotal", sum(FS$ci_above0)); add("FwdSigAll", sum(FS$n))
add("CausalRowsMin", format(min(CT$rows_compared), big.mark = "{,}")); add("CausalFeatures", CT$n_features[1]); add("CausalCuts", nrow(CT)); add("CausalMaxDiff", format(max(CT$max_abs_diff_all), digits = 3))
ax <- AL[model == "xgboost" & subset == "SXR+HXR+EUV" & target %in% c("y15", "y30", "y60")]; ah <- AL[model == "xgboost" & subset == "NONE" & target %in% c("y15", "y30", "y60")]
add("EpDayFullLo", f2(min(ax$episodes_per_day), 0)); add("EpDayFullHi", f2(max(ax$episodes_per_day), 0)); add("EpPrecFullLo", round(100 * min(ax$episode_precision))); add("EpPrecFullHi", round(100 * max(ax$episode_precision)))
add("EpDayHistLo", f2(min(ah$episodes_per_day), 0)); add("EpDayHistHi", f2(max(ah$episodes_per_day), 0)); add("EpPrecHistLo", round(100 * min(ah$episode_precision))); add("EpPrecHistHi", round(100 * max(ah$episode_precision)))
add("DutyLo", round(100 * min(AL[model == "xgboost" & subset %in% c("NONE", "SXR+HXR+EUV")]$alert_minute_share))); add("DutyHi", round(100 * max(AL[model == "xgboost" & subset %in% c("NONE", "SXR+HXR+EUV")]$alert_minute_share)))
add("FwdTSSCLo", f2(min(RF[target %in% c("y15", "y30", "y60"), .(m = mean(TSS)), by = .(subset, target)]$m))); add("FwdTSSCHi", f2(max(RF[target %in% c("y15", "y30", "y60"), .(m = mean(TSS)), by = .(subset, target)]$m)))
add("FwdTSSMLo", f2(min(RF[target %in% c("yM30", "yM60"), .(m = mean(TSS)), by = .(subset, target)]$m))); add("FwdTSSMHi", f2(max(RF[target %in% c("yM30", "yM60"), .(m = mean(TSS)), by = .(subset, target)]$m)))
fc <- CSF[, .(lo = min(shapley_mean), hi = max(shapley_mean)), by = band]
for (b in c("SXR", "HXR", "EUV")) { add(paste0("FwdShapLo", b), f2(fc[band == b]$lo, 3)); add(paste0("FwdShapHi", b), f2(fc[band == b]$hi, 3)) }
cat(m, sep = "\n", file = file.path(pap, "numbers.tex"), append = TRUE)

# ---- Table: forward-in-time TSS ----
extra <- intersect("HISTONLY", unique(RF$subset)); ord <- c(sub_order, extra, "level_only")
sk <- RF[subset %in% ord & target %in% TG, .(TSS = mean(TSS)), by = .(subset, target)]
w <- dcast(sk, subset ~ target, value.var = "TSS"); w <- w[match(ord, subset)]
rows <- apply(w, 1, function(r) paste0(ifelse(r[["subset"]] == "level_only", "SoLEXS level only", ifelse(r[["subset"]] == "HISTONLY", "History only (no magnetic)", sub_tex[[r[["subset"]]]])), " & ", paste(f2(as.numeric(r[TG])), collapse = " & "), " \\\\"))
writeLines(c("\\begin{tabular}{lccccc}", "\\toprule", paste0("Features & ", paste(tg_lab[TG], collapse = " & "), " \\\\"), "\\midrule", rows[1:8], "\\midrule", rows[9:length(rows)], "\\bottomrule", "\\end{tabular}"), file.path(pap, "tab_forward.tex"))

# ---- Table: paired-block gain vs history, both designs ----
tp <- merge(PS[, .(subset, loo_gain = mean_delta, loo_sig = ci_above0, loo_n = n)], FS[, .(subset, fwd_gain = mean_delta, fwd_sig = ci_above0, fwd_n = n)], by = "subset")[match(sub_order[-1], subset)]
rows <- apply(tp, 1, function(r) sprintf("%s & %s & %s/%s & %s & %s/%s \\\\", sub_tex[[r[["subset"]]]], f2(as.numeric(r[["loo_gain"]]), 3), r[["loo_sig"]], r[["loo_n"]], f2(as.numeric(r[["fwd_gain"]]), 3), r[["fwd_sig"]], r[["fwd_n"]]))
writeLines(c("\\begin{tabular}{lcccc}", "\\toprule", "Features & LOO gain & CI$>0$ & Fwd gain & CI$>0$ \\\\", "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), file.path(pap, "tab_paired.tex"))

# ---- Table: Shapley consensus, leave-one-block-out vs forward ----
rows <- c(); for (tg in TG) { cells <- vapply(c("SXR", "HXR", "EUV"), function(b) sprintf("%s / %s", f2(CS[target == tg & band == b]$shapley_mean, 3), f2(CSF[target == tg & band == b]$shapley_mean, 3)), "")
  rows <- c(rows, paste0(tg_lab[[tg]], " & ", paste(cells, collapse = " & "), " \\\\")) }
writeLines(c("\\begin{tabular}{lccc}", "\\toprule", "Target & SXR & HXR & EUV \\\\", "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), file.path(pap, "tab_shapley2.tex"))

# ---- Table: alert episodes ----
a <- AL[model == "xgboost" & subset %in% c("NONE", "SXR+HXR+EUV")][order(match(target, TG), subset)]
rows <- apply(a, 1, function(r) sprintf("%s & %s & %s & %s\\%% & %s\\%% & %s\\%% \\\\", tg_lab[[r[["target"]]]], ifelse(r[["subset"]] == "NONE", "History", "All bands"), f2(as.numeric(r[["episodes_per_day"]]), 1), round(100 * as.numeric(r[["episode_precision"]])), round(100 * as.numeric(r[["flare_detect_rate"]])), round(100 * as.numeric(r[["alert_minute_share"]]))))
writeLines(c("\\begin{tabular}{llrrrr}", "\\toprule", "Target & Features & Episodes/day & Precision & Flares caught & Alert on \\\\", "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), file.path(pap, "tab_alert.tex"))
for (f in c("fig_per_block_tss.png")) if (file.exists(file.path(OM, "figures", f))) file.copy(file.path(OM, "figures", f), file.path(pap, "figs", f), overwrite = TRUE)
message("assets2 written")
