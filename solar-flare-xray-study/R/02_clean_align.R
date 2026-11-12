# 02_clean_align.R -- read raw files, QC, align everything on a 1-minute UTC grid
# Outputs: data/processed/{goes_1min.csv, flare_events.csv, hxr_1min.csv, merged_1min.csv}
suppressPackageStartupMessages({ library(ncdf4); library(data.table) })
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
raw <- Sys.getenv("SOLAR_RAW", file.path(root, "data", "raw")); proc <- Sys.getenv("SOLAR_PROC", file.path(root, "data", "processed"))
dir.create(proc, FALSE, TRUE)
EPOCH <- as.POSIXct("2000-01-01 12:00:00", tz = "UTC")   # GOES L2 time origin

read_nc <- function(f, vars) {
  n <- nc_open(f); on.exit(nc_close(n))
  out <- list(time = EPOCH + ncvar_get(n, "time"))
  for (v in vars) out[[v]] <- as.numeric(ncvar_get(n, v))
  as.data.table(out)
}
bind_dir <- function(dir, pat, vars) rbindlist(lapply(list.files(dir, pattern = pat, full.names = TRUE),
                                                      function(f) tryCatch(read_nc(f, vars), error = function(e) { message("skip ", basename(f)); NULL })))

# ---- GOES XRS: SXR (1-8 A, xrsb) and harder SXR (0.5-4 A, xrsa) ----
xrs <- bind_dir(file.path(raw, "goes_xrs"), "\\.nc$", c("xrsa_flux", "xrsb_flux", "xrsa_flag", "xrsb_flag"))
xrs[xrsa_flag != 0 | xrsa_flux <= 0, xrsa_flux := NA_real_]      # flag != 0 means bad / excluded
xrs[xrsb_flag != 0 | xrsb_flux <= 0, xrsb_flux := NA_real_]
xrs[, c("xrsa_flag", "xrsb_flag") := NULL]

# ---- GOES EUVS: EUV lines 25.6, 28.4, 30.4 nm ----
euv <- bind_dir(file.path(raw, "goes_euvs"), "\\.nc$",
                c("irr_256", "irr_284", "irr_304", "irr_256_flag", "irr_284_flag", "irr_304_flag"))
for (l in c("256", "284", "304")) euv[get(paste0("irr_", l, "_flag")) != 0 | get(paste0("irr_", l)) <= 0, (paste0("irr_", l)) := NA_real_]
euv[, grep("_flag", names(euv), value = TRUE) := NULL]

goes <- merge(xrs, euv, by = "time", all.x = TRUE)
setorder(goes, time); goes <- unique(goes, by = "time")
fwrite(goes, file.path(proc, "goes_1min.csv"))

# ---- Flare events (ground truth) from NOAA flare summary ----
fl <- bind_dir(file.path(raw, "flare_catalog"), "\\.nc$", c("xrsb_flux", "background_flux", "integrated_flux", "flare_id"))
cls <- rbindlist(lapply(list.files(file.path(raw, "flare_catalog"), pattern = "\\.nc$", full.names = TRUE), function(f) {
  n <- nc_open(f); on.exit(nc_close(n))
  data.table(status = trimws(ncvar_get(n, "status")), flare_class = trimws(ncvar_get(n, "flare_class")))
}))
fl[, `:=`(status = cls$status, flare_class = cls$flare_class)]
ev <- fl[status == "EVENT_START", .(flare_id, start = time)][
  fl[status == "EVENT_PEAK", .(flare_id, peak = time, peak_flux = xrsb_flux, class = flare_class)], on = "flare_id"][
  fl[status == "EVENT_END", .(flare_id, end = time)], on = "flare_id"]
ev[, class_letter := substr(class, 1, 1)]
ev[, class_num := suppressWarnings(as.numeric(substring(class, 2)))]
setorder(ev, peak)
fwrite(ev, file.path(proc, "flare_events.csv"))
message("GOES minutes: ", nrow(goes), "; flares: ", nrow(ev), "; by class: ", paste(names(table(ev$class_letter)), table(ev$class_letter), collapse = ", "))

# ---- Fermi GBM hard X-ray (run once CTIME files are present) ----
# (implemented in R/02b_fermi_hxr.R)
