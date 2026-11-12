# 04_models.R -- band-subset ablation with three model families, chronological validation
# Output: outputs/tables/results.csv, outputs/models/preds.rds
suppressPackageStartupMessages({ library(data.table); library(glmnet); library(ranger); library(xgboost); library(pROC) })
set.seed(42)
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
proc <- file.path(root, "data", "processed"); out <- file.path(root, "outputs"); for (p in c("tables", "models", "figures")) dir.create(file.path(out, p), FALSE, TRUE)
to_utc <- function(v) if (inherits(v, "POSIXct")) { attr(v, "tzone") <- "UTC"; v } else as.POSIXct(as.character(v), tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ")
X <- fread(file.path(proc, "features.csv")); X[, time := to_utc(time)]
ev <- fread(file.path(proc, "flare_events.csv")); ev[, start := to_utc(start)]; ev <- ev[!is.na(start)]
ev <- ev[class_letter %in% c("C", "M", "X")]
D <- X[complete == TRUE]; setorder(D, time)

# ---- blocked (contiguous-in-time) cross-validation: K blocks; each block is the test set once, the next block
# is the validation set (threshold selection), the rest train. 60-minute purge gaps around test/val blocks so that
# labels and rolling windows cannot leak across a boundary. No random shuffling. ----
K <- 7; GAP <- 60 * 60
t0 <- as.numeric(min(D$time)); span <- as.numeric(max(D$time)) - t0
D[, fold := pmin(K, floor((as.numeric(time) - t0) / span * K) + 1)]
message("rows per fold: ", paste(D[, .N, by = fold][order(fold)]$N, collapse = "/"))
near <- function(tt, rng) tt >= rng[1] - GAP & tt <= rng[2] + GAP

# ---- feature subsets ----
BLOCKS <- c("SXR", "HXR", "EUV"); pref <- c(SXR = "^sxr_", HXR = "^hxr_", EUV = "^euv_")
subsets <- unlist(lapply(1:3, function(k) combn(BLOCKS, k, simplify = FALSE)), recursive = FALSE)
feat_cols <- function(S) { cols <- unlist(lapply(S, function(b) grep(pref[[b]], names(D), value = TRUE)))
  if (all(c("SXR", "HXR") %in% S)) cols <- c(cols, grep("^nep_", names(D), value = TRUE))  # Neupert features need both bands
  cols }
subset_name <- function(S) paste(S, collapse = "+")

# ---- metrics ----
confusion <- function(y, pred) { tp <- sum(y == 1 & pred == 1); fn <- sum(y == 1 & pred == 0); fp <- sum(y == 0 & pred == 1); tn <- sum(y == 0 & pred == 0); c(tp = tp, fn = fn, fp = fp, tn = tn) }
skill <- function(cm) { with(as.list(cm), { pod <- tp / max(1, tp + fn); pofd <- fp / max(1, fp + tn); far <- if (tp + fp == 0) NA_real_ else fp / (tp + fp)
  n <- tp + fn + fp + tn; ex <- ((tp + fn) * (tp + fp) + (tn + fn) * (tn + fp)) / n; hss <- if (n == ex) 0 else (tp + tn - ex) / (n - ex)
  c(POD = pod, POFD = pofd, FAR = far, TSS = pod - pofd, HSS = hss) }) }
best_thr <- function(y, p) { g <- sort(unique(quantile(p, probs = seq(0.50, 0.999, length.out = 300)))); s <- sapply(g, function(th) { cm <- confusion(y, as.integer(p >= th)); skill(cm)[["TSS"]] }); g[which.max(s)] }
pr_auc <- function(y, p) { o <- order(-p); y <- y[o]; tp <- cumsum(y); prec <- tp / seq_along(y); rec <- tp / max(1, sum(y)); sum(diff(c(0, rec)) * prec) }

# ---- event-level evaluation: for each flare start in the test period, was there an alert in [start-H, start)? lead time ----
event_eval <- function(times, alert, H, starts) {
  tt <- as.numeric(times); res <- rbindlist(lapply(starts, function(s) {
    idx <- which(tt >= s - H * 60 & tt < s); if (length(idx) < 0.5 * H) return(NULL)
    a <- idx[alert[idx] == 1]; data.table(hit = length(a) > 0, lead_min = if (length(a)) (s - tt[min(a)]) / 60 else NA_real_) }))
  if (!nrow(res)) res <- data.table(hit = logical(), lead_min = numeric())
  list(table = res) }

# ---- model fitters: return a prediction function ----
fit_model <- function(name, Xtr, ytr, Xva, yva) {
  wpos <- sum(ytr == 0) / max(1, sum(ytr == 1))
  if (name == "glmnet") {
    two <- function(M) { M <- as.matrix(M); if (ncol(M) == 1) cbind(M, M) else M }   # glmnet needs >= 2 columns (single-feature baseline)
    Xtr <- two(Xtr); mu <- colMeans(Xtr); sdv <- apply(Xtr, 2, sd); sdv[sdv == 0] <- 1; sc <- function(M) sweep(sweep(two(M), 2, mu), 2, sdv, "/")
    f <- cv.glmnet(sc(Xtr), ytr, family = "binomial", alpha = 0.5, weights = ifelse(ytr == 1, wpos, 1), nfolds = 5, type.measure = "auc", foldid = ceiling(seq_along(ytr) / (length(ytr) / 5)))  # contiguous folds
    return(function(M) as.numeric(predict(f, sc(M), s = "lambda.1se", type = "response")))
  }
  if (name == "ranger") {
    f <- ranger(x = as.data.frame(Xtr), y = factor(ytr, levels = c(0, 1)), probability = TRUE, num.trees = 300, min.node.size = 100, mtry = max(2, floor(sqrt(ncol(Xtr)))),
                class.weights = c(1, min(wpos, 20)), seed = 42, num.threads = 2)
    return(function(M) predict(f, as.data.frame(M))$predictions[, "1"])
  }
  if (name == "xgboost") {
    dtr <- xgb.DMatrix(as.matrix(Xtr), label = ytr); dva <- xgb.DMatrix(as.matrix(Xva), label = yva)
    par <- list(objective = "binary:logistic", eval_metric = "aucpr", max_depth = 3, eta = 0.05, subsample = 0.8, colsample_bytree = 0.8, min_child_weight = 20, scale_pos_weight = min(wpos, 20), nthread = 2)
    argn <- if ("evals" %in% names(formals(xgb.train))) "evals" else "watchlist"
    args <- list(params = par, data = dtr, nrounds = 400, early_stopping_rounds = 30, verbose = 0); args[[argn]] <- list(val = dva)
    f <- do.call(xgb.train, args)
    return(function(M) predict(f, xgb.DMatrix(as.matrix(M))))
  }
}

TARGETS <- c("y15", "y30", "y60", "yM30", "yM60"); HZ <- c(y15 = 15, y30 = 30, y60 = 60, yM30 = 30, yM60 = 60); MODELS <- c("glmnet", "ranger", "xgboost")
flare_starts <- function(tgt) { s <- if (grepl("^yM", tgt)) ev[class_letter %in% c("M", "X")]$start else ev$start; as.numeric(s) }
res <- list(); preds <- list()
configs <- c(lapply(subsets, function(S) list(name = subset_name(S), cols = feat_cols(S))), list(list(name = "baseline_level", cols = "sxr_b_log")))
for (tgt in TARGETS) {
  message("== ", tgt, " positives per fold: ", paste(D[, sum(get(tgt)), by = fold][order(fold)]$V1, collapse = "/"))
  if (sum(D[[tgt]]) < 100) { message("skip ", tgt); next }
  st_all <- flare_starts(tgt)
  for (cfg in configs) for (m in MODELS) {
    if (cfg$name == "baseline_level" && m != "glmnet") next
    P <- list(); E <- list(); failed <- FALSE
    for (k in seq_len(K)) {
      vf <- k %% K + 1; te <- D[fold == k]; va <- D[fold == vf]; tr <- D[!(fold %in% c(k, vf))]
      rte <- range(as.numeric(te$time)); rva <- range(as.numeric(va$time))
      tr <- tr[!near(as.numeric(time), rte) & !near(as.numeric(time), rva)]; va <- va[!near(as.numeric(time), rte)]   # purge
      if (sum(tr[[tgt]]) < 20 || sum(va[[tgt]]) < 5) next
      pf <- tryCatch(fit_model(m, tr[, cfg$cols, with = FALSE], tr[[tgt]], va[, cfg$cols, with = FALSE], va[[tgt]]), error = function(e) { message("fit fail ", m, " ", cfg$name, ": ", e$message); NULL })
      if (is.null(pf)) { failed <- TRUE; break }
      pv <- pf(va[, cfg$cols, with = FALSE]); pt <- pf(te[, cfg$cols, with = FALSE]); thr <- best_thr(va[[tgt]], pv)
      P[[k]] <- data.table(fold = k, time = te$time, y = te[[tgt]], p = pt, thr = thr)
      ee <- event_eval(te$time, as.integer(pt >= thr), HZ[[tgt]], st_all[st_all >= rte[1] & st_all <= rte[2]]); E[[k]] <- ee$table
    }
    if (failed || !length(P)) next
    pp <- rbindlist(P); et <- rbindlist(E); sk <- skill(confusion(pp$y, as.integer(pp$p >= pp$thr)))
    fold_tss <- sapply(split(pp, pp$fold), function(d) skill(confusion(d$y, as.integer(d$p >= d$thr)))[["TSS"]])
    res[[length(res) + 1]] <- data.table(target = tgt, horizon = HZ[[tgt]], model = m, subset = cfg$name, n_feat = length(cfg$cols), n_test = nrow(pp),
      t(sk), TSS_fold_sd = sd(fold_tss), AUC = as.numeric(auc(roc(pp$y, pp$p, quiet = TRUE))), PRAUC = pr_auc(pp$y, pp$p), prevalence = mean(pp$y),
      n_flares = nrow(et), event_POD = mean(et$hit), lead_median = median(et$lead_min, na.rm = TRUE), lead_mean = mean(et$lead_min, na.rm = TRUE))
    preds[[length(preds) + 1]] <- data.table(target = tgt, model = m, subset = cfg$name, pp)
    message(sprintf("  %-8s %-16s TSS=%.3f (fold sd %.3f) POD=%.2f FAR=%.2f lead=%.1f  [n flares %d]", m, cfg$name, sk[["TSS"]], sd(fold_tss), sk[["POD"]], sk[["FAR"]], median(et$lead_min, na.rm = TRUE), nrow(et)))
  }
}
R <- rbindlist(res); fwrite(R, file.path(out, "tables", "results.csv")); saveRDS(rbindlist(preds), file.path(out, "models", "preds.rds"))
message("saved ", nrow(R), " result rows")
