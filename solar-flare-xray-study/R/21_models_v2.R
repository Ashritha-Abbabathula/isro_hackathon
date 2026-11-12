# 21_models_v2.R -- (v2 of 11_models_adl1.R: baseline = flare history + magnetic state; extra HISTONLY config) band-subset ablation on Aditya-L1 (+GOES EUV) data; blocked CV; parallel over experiments.
# Env: SOLAR_PROC (input dir), SOLAR_OUT (output dir), SEP_MODE = exclude (default) | include, MODELS (comma list)
suppressPackageStartupMessages({ library(data.table); library(parallel) })
set.seed(42)
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
proc <- Sys.getenv("SOLAR_PROC", file.path(root, "data", "processed_v2")); out <- Sys.getenv("SOLAR_OUT", file.path(root, "outputs_v2"))
MODE <- Sys.getenv("MODE", "loo"); SEP_MODE <- Sys.getenv("SEP_MODE", "exclude"); MODELS <- strsplit(Sys.getenv("MODELS", "glmnet,ranger,xgboost"), ",")[[1]]
for (p in c("tables", "models", "figures")) dir.create(file.path(out, p), FALSE, TRUE)
source(file.path(root, "R", "lib_models.R"))
to_utc <- function(v) if (inherits(v, "POSIXct")) { attr(v, "tzone") <- "UTC"; v } else as.POSIXct(as.character(v), tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ")
X <- fread(file.path(proc, "features_v2.csv")); X[, time := to_utc(time)]
ev <- fread(file.path(proc, "flare_events.csv")); ev[, start := to_utc(start)]; ev <- ev[!is.na(start) & class_letter %in% c("C", "M", "X")]
D <- if (SEP_MODE == "exclude") X[complete == TRUE & sep == FALSE] else X[complete == TRUE]; setorder(D, time)

K <- 8; GAP <- 60 * 60; STRIDE <- 4
t0 <- as.numeric(min(D$time)); span <- as.numeric(max(D$time)) - t0
D[, fold := pmin(K, floor((as.numeric(time) - t0) / span * K) + 1)]
message("MODE=", MODE, " | SEP_MODE=", SEP_MODE, " | rows: ", nrow(D), " | per fold: ", paste(D[, .N, by = fold][order(fold)]$N, collapse = "/"))

BANDS <- c("SXR", "HXR", "EUV"); pref <- c(SXR = "^sxr_", HXR = "^hxr_", EUV = "^euv_")
hist_cols <- grep("^hist_", names(D), value = TRUE)
mag_cols <- grep("^mag_", names(D), value = TRUE)
feat_cols <- function(S) { cols <- c(hist_cols, mag_cols)                              # flare history AND magnetic state are in EVERY band configuration
  for (b in S) cols <- c(cols, grep(pref[[b]], names(D), value = TRUE))
  if (all(c("SXR", "HXR") %in% S)) cols <- c(cols, grep("^nep_", names(D), value = TRUE)); cols }
subs <- c(list(character(0)), unlist(lapply(1:3, function(k) combn(BANDS, k, simplify = FALSE)), recursive = FALSE))
configs <- c(lapply(subs, function(S) list(name = if (length(S)) paste(S, collapse = "+") else "NONE", cols = feat_cols(S))),
             list(list(name = "HISTONLY", cols = hist_cols), list(name = "level_only", cols = "sxr_sol_log")))

TARGETS <- c("y15", "y30", "y60", "yM30", "yM60"); HZ <- list(y15 = 15, y30 = 30, y60 = 60, yM30 = 30, yM60 = 60)
tasks <- list()
for (tg in TARGETS) { if (sum(D[[tg]]) < 100) { message("skip ", tg); next }
  for (cfg in configs) for (m in MODELS) { if (cfg$name == "level_only" && m != "glmnet") next
    tasks[[length(tasks) + 1]] <- list(tgt = tg, name = cfg$name, cols = cfg$cols, model = m) } }
starts <- list(C = as.numeric(ev$start), M = as.numeric(ev[class_letter %in% c("M", "X")]$start))
message("experiments: ", length(tasks))

nw <- max(1, min(7, detectCores() - 1)); cl <- makeCluster(nw)
clusterExport(cl, c("D", "starts", "K", "GAP", "STRIDE", "HZ", "root", "MODE"), envir = environment())
clusterEvalQ(cl, source(file.path(root, "R", "lib_models.R")))
t_start <- Sys.time()
out_list <- parLapplyLB(cl, tasks, function(tk) {
  r <- tryCatch(run_task(tk, D, if (grepl("^yM", tk$tgt)) starts$M else starts$C, K, GAP, STRIDE, HZ, MODE), error = function(e) list(error = paste(tk$tgt, tk$name, tk$model, conditionMessage(e))))
  r })
stopCluster(cl)
errs <- Filter(function(x) !is.null(x) && !is.null(x$error), out_list); if (length(errs)) { message("FAILED experiments: ", length(errs)); for (e in head(errs, 5)) message("  ", e$error) }
null_n <- sum(vapply(out_list, is.null, logical(1))); if (null_n) message("NULL (skipped) experiments: ", null_n)
ok <- Filter(function(x) !is.null(x) && is.null(x$error), out_list)
R <- rbindlist(lapply(ok, `[[`, "res")); Pd <- rbindlist(lapply(ok, `[[`, "pred"))
fwrite(R, file.path(out, "tables", "results.csv")); saveRDS(Pd, file.path(out, "models", "preds.rds"))
message("done ", nrow(R), "/", length(tasks), " experiments in ", round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), " min")
