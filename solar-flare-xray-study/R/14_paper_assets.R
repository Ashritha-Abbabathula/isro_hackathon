# 14_paper_assets.R -- turn result tables into LaTeX tables + number macros for the IEEE paper (paper/ieee/)
suppressPackageStartupMessages({ library(data.table) })
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
out <- Sys.getenv("OUT_MAIN", file.path(root, "outputs_adl1")); outS <- Sys.getenv("OUT_SEP", file.path(root, "outputs_adl1_sepincl")); pap <- Sys.getenv("PAPER_DIR", file.path(root, "paper", "ieee")); fig <- file.path(pap, "figs")
dir.create(fig, FALSE, TRUE); tb <- file.path(out, "tables")
f2 <- function(x, d = 2) ifelse(is.na(x), "--", formatC(x, format = "f", digits = d)); pct <- function(x) paste0(round(100 * x), "\\%")
R <- fread(file.path(tb, "results.csv")); SH <- fread(file.path(tb, "shapley_band_contribution.csv")); DL <- fread(file.path(tb, "delta_vs_history.csv"))
CS <- fread(file.path(tb, "shapley_consensus.csv")); NL <- fread(file.path(tb, "nowcast_latency_summary.csv")); NP <- fread(file.path(tb, "neupert_per_flare_spearman.csv"))
FC <- fread(file.path(tb, "flare_counts_window.csv")); X <- fread(file.path(Sys.getenv("PROC_DIR", file.path(root, "data", "processed_adl1")), Sys.getenv("FEATURE_FILE", "features_adl1.csv")))
tg_lab <- c(y15 = "C+ 15", y30 = "C+ 30", y60 = "C+ 60", yM30 = "M+ 30", yM60 = "M+ 60"); TG <- names(tg_lab)
sub_order <- c("NONE", "SXR", "HXR", "EUV", "SXR+HXR", "SXR+EUV", "HXR+EUV", "SXR+HXR+EUV")
sub_tex <- c(NONE = Sys.getenv("BASELINE_LABEL", "History only"), SXR = "+SXR", HXR = "+HXR", EUV = "+EUV", "SXR+HXR" = "+SXR+HXR", "SXR+EUV" = "+SXR+EUV", "HXR+EUV" = "+HXR+EUV", "SXR+HXR+EUV" = "+SXR+HXR+EUV")

# ---- number macros ----
m <- c()
add <- function(name, val) { name <- gsub("15", "Fifteen", gsub("30", "Thirty", gsub("60", "Sixty", name))); m <<- c(m, sprintf("\\newcommand{\\%s}{%s}", name, val)) }  # LaTeX macro names cannot contain digits
add("NumDays", round(nrow(X) / 1440)); add("NumMinutes", format(nrow(X), big.mark = "{,}")); add("NumUsable", format(sum(X$complete & !X$sep), big.mark = "{,}"))
add("NumSepMin", format(sum(X$sep), big.mark = "{,}")); add("NumFlaresAll", sum(FC$N))
for (cl in c("C", "M", "X")) add(paste0("NumFlares", cl), ifelse(cl %in% FC$class_letter, FC[class_letter == cl]$N, 0))
add("NeuMedian", f2(median(NP$r))); add("NeuN", nrow(NP)); add("NeuPos", round(100 * mean(NP$r > 0)))
for (tg in TG) { r <- R[target == tg & subset == "SXR+HXR+EUV" & model == "xgboost"]; if (nrow(r)) { add(paste0("FullTSS", sub("y", "", tg)), f2(r$TSS)); add(paste0("FullFAR", sub("y", "", tg)), f2(r$FAR)); add(paste0("FullPOD", sub("y", "", tg)), f2(r$POD)) }
  h <- R[target == tg & subset == "NONE" & model == "xgboost"]; if (nrow(h)) add(paste0("HistTSS", sub("y", "", tg)), f2(h$TSS)) }
add("PosYfifteen", format(sum(X[complete & !sep]$y15), big.mark = "{,}")); add("PosYthirty", format(sum(X[complete & !sep]$y30), big.mark = "{,}")); add("PosYsixty", format(sum(X[complete & !sep]$y60), big.mark = "{,}"))
add("PosMthirty", format(sum(X[complete & !sep]$yM30), big.mark = "{,}")); add("PosMsixty", format(sum(X[complete & !sep]$yM60), big.mark = "{,}"))
add("MedFoldSD", f2(median(R$TSS_fold_sd, na.rm = TRUE)))
add("SigAll", nrow(DL)); add("SigTotal", sum(DL$lo > 0)); for (s in c("SXR", "HXR", "EUV", "SXR+HXR")) add(paste0("Sig", sub("+", "", s, fixed = TRUE)), sum(DL$lo > 0 & DL$subset == s)); add("SigPerSubset", nrow(DL[subset == "SXR"]))
cs <- CS[, .(lo = min(shapley_mean), hi = max(shapley_mean)), by = band]
for (b in c("SXR", "HXR", "EUV")) { add(paste0("ShapLo", b), f2(cs[band == b]$lo, 3)); add(paste0("ShapHi", b), f2(cs[band == b]$hi, 3)) }
hs <- R[subset == "NONE", .(m = mean(TSS)), by = target]; add("HistLoC", f2(min(hs[target %in% c("y15", "y30", "y60")]$m))); add("HistHiC", f2(max(hs[target %in% c("y15", "y30", "y60")]$m))); add("HistM", f2(mean(hs[target %in% c("yM30", "yM60")]$m)))
sk2 <- R[subset %in% sub_order, .(m = mean(TSS)), by = .(target, subset)]
add("RangeLoC", f2(min(sk2[target %in% c("y15", "y30", "y60")]$m))); add("RangeHiC", f2(max(sk2[target %in% c("y15", "y30", "y60")]$m)))
add("RangeLoM", f2(min(sk2[target %in% c("yM30", "yM60")]$m))); add("RangeHiM", f2(max(sk2[target %in% c("yM30", "yM60")]$m)))
lvm <- R[subset == "level_only" & target %in% c("yM30", "yM60")]; add("LevelM30", f2(lvm[target == "yM30"]$TSS)); add("LevelM60", f2(lvm[target == "yM60"]$TSS))
far <- R[subset == "SXR+HXR+EUV" & model == "xgboost"]; add("FARlo", f2(min(far$FAR))); add("FARhi", f2(max(far$FAR)))
add("NumFlaresEvalC", R[target == "y30"]$n_flares[1]); add("NumFlaresEvalM", R[target == "yM30"]$n_flares[1])
PW <- fread(file.path(tb, "nowcast_paired.csv")); for (g in c("C", "MX")) { r <- PW[class_group == ifelse(g == "C", "C", "M/X")]; add(paste0("PairN", g), r$n_both); add(paste0("PairMed", g), f2(r$median_diff_min, 1)); add(paste0("PairEarly", g), round(100 * r$hel_earlier)); add(paste0("PairSame", g), round(100 * r$same_minute)); add(paste0("PairLate", g), round(100 * r$hel_later)) }
writeLines(m, file.path(pap, "numbers.tex"))

# ---- Table: dataset ----
writeLines(c("\\begin{tabular}{lll}", "\\toprule", "Band & Instrument & Product used \\\\", "\\midrule",
  "SXR & Aditya-L1 SoLEXS (SDD2) & Level-1 light curve, 1 s $\\rightarrow$ 1 min \\\\",
  "HXR & Aditya-L1 HEL1OS CZT1+2 & 20--40, 40--60, 60--150 keV, $\\approx$8 s readout $\\rightarrow$ 1 min \\\\",
  "EUV & GOES-16 EUVS & 25.6, 28.4, 30.4 nm lines, 1 min \\\\",
  "Labels & GOES-16 XRS & NOAA flare summary (C, M, X) \\\\", "\\bottomrule", "\\end{tabular}"), file.path(pap, "tab_data.tex"))

# ---- Table: skill (mean TSS across the three model families) ----
sk <- R[subset %in% sub_order & target %in% TG, .(TSS = mean(TSS)), by = .(subset, target)]
sk_w <- dcast(sk, subset ~ target, value.var = "TSS"); sk_w <- sk_w[match(sub_order, subset)]
body <- apply(sk_w, 1, function(r) paste0(sub_tex[[r[["subset"]]]], " & ", paste(f2(as.numeric(r[TG])), collapse = " & "), " \\\\"))
if ("HISTONLY" %in% R$subset) { hh <- R[subset == "HISTONLY", .(TSS = mean(TSS)), by = target][match(TG, target)]; body <- c(body, paste0("History only (no magnetic) & ", paste(f2(hh$TSS), collapse = " & "), " \\\\")) }
lv <- R[subset == "level_only", .(TSS = mean(TSS)), by = target]; lv <- lv[match(TG, target)]
writeLines(c("\\begin{tabular}{lccccc}", "\\toprule", paste0("Features & ", paste(tg_lab[TG], collapse = " & "), " \\\\"), "\\midrule", body,
  "\\midrule", paste0("SoLEXS level only & ", paste(f2(lv$TSS), collapse = " & "), " \\\\"), "\\bottomrule", "\\end{tabular}"), file.path(pap, "tab_skill.tex"))

# ---- Table: Shapley consensus + gain over history ----
bb <- c("SXR", "HXR", "EUV"); rows <- c()
for (tg in TG) { cs <- CS[target == tg][match(bb, band)]
  cells <- vapply(seq_along(bb), function(i) sprintf("%s$\\pm$%s%s", f2(cs$shapley_mean[i], 3), f2(cs$shapley_sd[i], 3), ifelse(cs$n_models_ci_above0[i] >= 2, "$^{*}$", "")), "")
  g <- SH[target == tg, .(gain = mean(gain))]
  rows <- c(rows, paste0(tg_lab[[tg]], " & ", paste(cells, collapse = " & "), " & ", f2(g$gain, 3), " \\\\")) }
writeLines(c("\\begin{tabular}{lcccc}", "\\toprule", "Target & SXR & HXR & EUV & Gain \\\\", "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), file.path(pap, "tab_shapley.tex"))

# ---- Table: nowcasting latency ----
NL[, band := gsub("\\(SXR\\)", "", band)]
NL[, band := gsub("-", "--", band)]; NL[, bo := match(trimws(band), c("SoLEXS", "HEL1OS 20--40 keV", "HEL1OS 60--150 keV", "GOES EUV 30.4 nm"))]
rows <- apply(NL[order(match(class_group, c("C", "M/X")), bo)], 1, function(r) sprintf("%s & %s & %s & %s & %s \\\\", r[["class_group"]], gsub("_", "\\\\_", trimws(r[["band"]])), r[["n"]], pct(as.numeric(r[["detected"]])),
  ifelse(is.na(as.numeric(r[["median_latency_min"]])), "--", sprintf("%s [%s, %s]", f2(as.numeric(r[["median_latency_min"]]), 1), f2(as.numeric(r[["q25"]]), 1), f2(as.numeric(r[["q75"]]), 1)))))
writeLines(c("\\begin{tabular}{llrrr}", "\\toprule", "Class & Band & $n$ & Detected & Median [IQR] (min) \\\\", "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), file.path(pap, "tab_nowcast.tex"))

# ---- SEP sensitivity (if available) ----
if (file.exists(file.path(outS, "tables", "results.csv"))) {
  RS <- fread(file.path(outS, "tables", "results.csv")); a <- R[model %in% unique(RS$model) & subset %in% c("NONE", "SXR+HXR+EUV"), .(TSS_excl = mean(TSS)), by = .(target, subset)]
  b <- RS[subset %in% c("NONE", "SXR+HXR+EUV"), .(TSS_incl = mean(TSS)), by = .(target, subset)]; s <- merge(a, b, by = c("target", "subset")); fwrite(s, file.path(tb, "sensitivity_sep.csv"))
  rws <- apply(s[order(match(target, TG), subset)], 1, function(r) sprintf("%s & %s & %s & %s \\\\", tg_lab[[r[["target"]]]], ifelse(r[["subset"]] == "NONE", Sys.getenv("BASELINE_LABEL", "History only"), "All bands"), f2(as.numeric(r[["TSS_excl"]])), f2(as.numeric(r[["TSS_incl"]]))))
  writeLines(c("\\begin{tabular}{llcc}", "\\toprule", "Target & Features & SEP excluded & SEP included \\\\", "\\midrule", rws, "\\bottomrule", "\\end{tabular}"), file.path(pap, "tab_sep.tex")) }


# ---- SEP-sensitivity macros (appended to numbers.tex) ----
if (file.exists(file.path(outS, "tables", "shapley_consensus.csv"))) {
  m <- c(); RS2 <- fread(file.path(outS, "tables", "results.csv")); CS2 <- fread(file.path(outS, "tables", "shapley_consensus.csv"))
  a2 <- R[model %in% unique(RS2$model) & subset %in% sub_order, .(e = mean(TSS)), by = .(target, subset)]; b2 <- RS2[subset %in% sub_order, .(i = mean(TSS)), by = .(target, subset)]
  s2 <- merge(a2, b2, by = c("target", "subset")); add("SepMeanAbs", f2(mean(abs(s2$i - s2$e)), 3)); add("SepMaxAbs", f2(max(abs(s2$i - s2$e)), 3))
  for (b in c("SXR", "HXR", "EUV")) { add(paste0("SepShapLo", b), f2(min(CS2[band == b]$shapley_mean), 3)); add(paste0("SepShapHi", b), f2(max(CS2[band == b]$shapley_mean), 3)) }
  cat(m, sep = "
", file = file.path(pap, "numbers.tex"), append = TRUE) }
# ---- figures ----
for (f in c("fig_example_flare.png", "fig_neupert.png", "fig_nowcast_latency.png", "fig_shapley_bands.png", "fig_ablation_tss.png")) if (file.exists(file.path(out, "figures", f))) file.copy(file.path(out, "figures", f), file.path(fig, f), overwrite = TRUE)
message("paper assets written to ", pap)
