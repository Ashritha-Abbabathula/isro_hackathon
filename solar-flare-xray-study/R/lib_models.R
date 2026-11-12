# lib_models.R -- metrics, model fitters and one cross-validated task; sourced by the main script and by every parallel worker.
suppressPackageStartupMessages({ library(data.table); library(glmnet); library(ranger); library(xgboost); library(pROC) })

confusion <- function(y, pred) as.numeric(c(tp = sum(y == 1 & pred == 1), fn = sum(y == 1 & pred == 0), fp = sum(y == 0 & pred == 1), tn = sum(y == 0 & pred == 0))) |> setNames(c("tp", "fn", "fp", "tn"))  # double: avoids integer overflow in products
skill <- function(cm) { with(as.list(cm), { pod <- tp / max(1, tp + fn); pofd <- fp / max(1, fp + tn); far <- if (tp + fp == 0) NA_real_ else fp / (tp + fp)
  n <- tp + fn + fp + tn; ex <- ((tp + fn) * (tp + fp) + (tn + fn) * (tn + fp)) / n; hss <- if (n == ex) 0 else (tp + tn - ex) / (n - ex)
  c(POD = pod, POFD = pofd, FAR = far, TSS = pod - pofd, HSS = hss) }) }
best_thr <- function(y, p) { g <- sort(unique(quantile(p, probs = seq(0.50, 0.999, length.out = 300), na.rm = TRUE)))
  s <- vapply(g, function(th) skill(confusion(y, as.integer(p >= th)))[["TSS"]], numeric(1)); g[which.max(s)] }
pr_auc <- function(y, p) { o <- order(-p); y <- y[o]; tp <- cumsum(y); prec <- tp / seq_along(y); rec <- tp / max(1, sum(y)); sum(diff(c(0, rec)) * prec) }

# event level: for every flare start s in the test block, was there an alert in [s-H, s)?  lead = s - first alert
event_eval <- function(times, alert, H, starts) {
  tt <- as.numeric(times)
  res <- rbindlist(lapply(starts, function(s) { idx <- which(tt >= s - H * 60 & tt < s); if (length(idx) < 0.5 * H) return(NULL)
    a <- idx[alert[idx] == 1]; data.table(hit = length(a) > 0, lead_min = if (length(a)) (s - tt[min(a)]) / 60 else NA_real_) }))
  if (!nrow(res)) res <- data.table(hit = logical(), lead_min = numeric()); res }

fit_model <- function(name, Xtr, ytr, Xva, yva) {
  wpos <- sum(ytr == 0) / max(1, sum(ytr == 1))
  if (name == "glmnet") {
    two <- function(M) { M <- as.matrix(M); if (ncol(M) == 1) cbind(M, M) else M }
    Xtr <- two(Xtr); mu <- colMeans(Xtr); sdv <- apply(Xtr, 2, sd); sdv[!is.finite(sdv) | sdv == 0] <- 1
    sc <- function(M) sweep(sweep(two(M), 2, mu), 2, sdv, "/")
    f <- cv.glmnet(sc(Xtr), ytr, family = "binomial", alpha = 0.5, weights = ifelse(ytr == 1, min(wpos, 20), 1), nfolds = 5, type.measure = "auc", foldid = ceiling(seq_along(ytr) / (length(ytr) / 5)))
    return(function(M) as.numeric(predict(f, sc(M), s = "lambda.1se", type = "response")))
  }
  if (name == "ranger") {
    f <- ranger(x = as.data.frame(Xtr), y = factor(ytr, levels = c(0, 1)), probability = TRUE, num.trees = 200, min.node.size = 100, mtry = max(2, floor(sqrt(ncol(Xtr)))),
                class.weights = c(1, min(wpos, 20)), seed = 42, num.threads = 1)
    return(function(M) predict(f, as.data.frame(M), num.threads = 1)$predictions[, "1"])
  }
  if (name == "xgboost") {
    dtr <- xgb.DMatrix(as.matrix(Xtr), label = ytr); dva <- xgb.DMatrix(as.matrix(Xva), label = yva)
    par <- list(objective = "binary:logistic", eval_metric = "aucpr", max_depth = 3, eta = 0.05, subsample = 0.8, colsample_bytree = 0.8, min_child_weight = 20, scale_pos_weight = min(wpos, 20), nthread = 1)
    argn <- if ("evals" %in% names(formals(xgb.train))) "evals" else "watchlist"
    args <- list(params = par, data = dtr, nrounds = 300, early_stopping_rounds = 30, verbose = 0); args[[argn]] <- list(val = dva)
    f <- do.call(xgb.train, args)
    return(function(M) predict(f, xgb.DMatrix(as.matrix(M))))
  }
}

# one (target, configuration, model) experiment under K-fold blocked CV -> list(res, pred)
run_task <- function(task, D, starts_all, K, GAP, stride, HZ, mode = "loo") {  # mode "loo": leave-one-block-out; "forward": train only on EARLIER blocks (expanding window)
  tgt <- task$tgt; cols <- task$cols; m <- task$model
  near <- function(tt, rng) tt >= rng[1] - GAP & tt <= rng[2] + GAP
  P <- list(); E <- list()
  for (k in seq_len(K)) {
    if (mode == "forward") { if (k < 3) next; vf <- k - 1; te <- D[fold == k]; va <- D[fold == vf]; tr <- D[fold < vf] }
    else { vf <- k %% K + 1; te <- D[fold == k]; va <- D[fold == vf]; tr <- D[!(fold %in% c(k, vf))] }
    rte <- range(as.numeric(te$time)); rva <- range(as.numeric(va$time))
    tr <- tr[!near(as.numeric(time), rte) & !near(as.numeric(time), rva)]; va <- va[!near(as.numeric(time), rte)]
    tr <- tr[seq(1, .N, by = stride)]
    if (sum(tr[[tgt]]) < 20 || sum(va[[tgt]]) < 5 || sum(te[[tgt]]) < 1) next
    pf <- tryCatch(fit_model(m, tr[, cols, with = FALSE], tr[[tgt]], va[, cols, with = FALSE], va[[tgt]]), error = function(e) NULL)
    if (is.null(pf)) return(NULL)
    pv <- pf(va[, cols, with = FALSE]); pt <- pf(te[, cols, with = FALSE]); thr <- best_thr(va[[tgt]], pv)
    P[[k]] <- data.table(fold = k, time = te$time, y = te[[tgt]], p = pt, thr = thr)
    E[[k]] <- event_eval(te$time, as.integer(pt >= thr), HZ[[tgt]], starts_all[starts_all >= rte[1] & starts_all <= rte[2]])
  }
  if (!length(P)) return(NULL)
  pp <- rbindlist(P); et <- rbindlist(E); sk <- skill(confusion(pp$y, as.integer(pp$p >= pp$thr)))
  fold_tss <- vapply(split(pp, pp$fold), function(d) skill(confusion(d$y, as.integer(d$p >= d$thr)))[["TSS"]], numeric(1))
  res <- data.table(target = tgt, horizon = HZ[[tgt]], model = m, subset = task$name, n_feat = length(cols), n_test = nrow(pp), t(sk), TSS_fold_sd = sd(fold_tss),
    AUC = as.numeric(auc(roc(pp$y, pp$p, quiet = TRUE))), PRAUC = pr_auc(pp$y, pp$p), prevalence = mean(pp$y),
    n_flares = nrow(et), event_POD = mean(et$hit), lead_median = median(et$lead_min, na.rm = TRUE), lead_mean = mean(et$lead_min, na.rm = TRUE))
  list(res = res, pred = data.table(target = tgt, model = m, subset = task$name, pp))
}
