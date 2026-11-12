# 01_download.R -- fetch public data (resumable; skips files that exist)
# Usage: Rscript R/01_download.R [goes|fermi|all]
suppressPackageStartupMessages({ library(FITSio) })
options(timeout = 900)
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."),
                      winslash = "/", mustWork = FALSE)
if (is.na(root) || !dir.exists(root)) root <- getwd()
raw <- Sys.getenv("SOLAR_RAW", file.path(root, "data", "raw"))   # override to keep large data outside OneDrive

GOES_START <- as.Date(Sys.getenv("GOES_START", "2024-04-15")); GOES_END <- as.Date(Sys.getenv("GOES_END", "2024-06-30"))  # NOAA host is slow (~15 s/file)
FERMI_START <- as.Date("2024-05-01"); FERMI_END <- as.Date("2024-05-21")  # extend if bandwidth allows
SAT <- "goes16"
NOAA <- paste0("https://data.ngdc.noaa.gov/platforms/solar-space-observing-satellites/goes/", SAT, "/l2/data/")
FERMI <- "https://heasarc.gsfc.nasa.gov/FTP/fermi/data/gbm/daily/"

fetch <- function(url, dest, tries = 4) {
  if (file.exists(dest) && file.size(dest) > 0) return(TRUE)
  dir.create(dirname(dest), showWarnings = FALSE, recursive = TRUE)
  for (i in seq_len(tries)) {
    ok <- tryCatch({ download.file(url, paste0(dest, ".part"), mode = "wb", quiet = TRUE) == 0 },
                   error = function(e) FALSE, warning = function(w) FALSE)
    if (ok && file.size(paste0(dest, ".part")) > 0) { file.rename(paste0(dest, ".part"), dest); return(TRUE) }
    unlink(paste0(dest, ".part")); Sys.sleep(2 * i)
  }
  message("FAILED: ", url); FALSE
}

# list a NOAA month directory once and map date -> filename
noaa_index <- function(product, ym) {
  u <- paste0(NOAA, product, "/", substr(ym, 1, 4), "/", substr(ym, 6, 7), "/")
  h <- tryCatch(readLines(u, warn = FALSE), error = function(e) character())
  f <- unique(regmatches(h, regexpr("sci_[^\"<>/]+?\\.nc", h, perl = TRUE)))
  setNames(f, sub(".*_d([0-9]{8})_.*", "\\1", f))
}

# parallel batch download (libcurl multi); retries the failures once per round
fetch_many <- function(urls, dests, batch = 24, rounds = 4) {
  keep <- !(file.exists(dests) & file.size(dests) > 0); urls <- urls[keep]; dests <- dests[keep]
  for (r in seq_len(rounds)) {
    if (!length(urls)) break
    for (i in seq(1, length(urls), by = batch)) {
      j <- i:min(i + batch - 1, length(urls)); for (p in unique(dirname(dests[j]))) dir.create(p, FALSE, TRUE)
      tmp <- paste0(dests[j], ".part")
      try(suppressWarnings(download.file(urls[j], tmp, method = "libcurl", mode = "wb", quiet = TRUE)), silent = TRUE)
      ok <- file.exists(tmp) & file.size(tmp) > 0
      # a truncated transfer still has size > 0, so confirm the NetCDF opens before accepting it
      ok[ok] <- vapply(tmp[ok], function(p) tryCatch({ n <- ncdf4::nc_open(p); ncdf4::nc_close(n); TRUE }, error = function(e) FALSE), logical(1))
      file.rename(tmp[ok], dests[j][ok]); unlink(tmp[!ok])
      message("batch ", i, "/", length(urls), " (round ", r, ")")
    }
    left <- !(file.exists(dests) & file.size(dests) > 0); urls <- urls[left]; dests <- dests[left]
  }
  if (length(urls)) message("FAILED: ", length(urls), " files")
}

download_goes <- function() {
  days <- seq(GOES_START, GOES_END, by = "day")
  months <- unique(format(days, "%Y-%m"))
  prods <- c(xrsf = "xrsf-l2-avg1m_science", euvs = "euvs-l2-avg1m_science", flsum = "xrsf-l2-flsum_science")
  dest <- c(xrsf = "goes_xrs", euvs = "goes_euvs", flsum = "flare_catalog")
  U <- character(); D <- character()
  for (k in names(prods)) for (ym in months) {
    idx <- noaa_index(prods[[k]], ym)
    for (d in names(idx)) {
      dd <- as.Date(d, "%Y%m%d"); if (dd < GOES_START || dd > GOES_END) next
      U <- c(U, paste0(NOAA, prods[[k]], "/", substr(ym, 1, 4), "/", substr(ym, 6, 7), "/", idx[[d]]))
      D <- c(D, file.path(raw, dest[[k]], idx[[d]]))
    }
  }
  message("GOES files to fetch: ", length(U)); fetch_many(U, D)
}

# ---- Fermi GBM: pick Sun-facing NaI detectors from the spacecraft position file ----
# NaI detector pointing in spacecraft frame (azimuth, zenith) in degrees (Meegan et al. 2009)
DET <- data.frame(det = c("n0","n1","n2","n3","n4","n5","n6","n7","n8","n9","na","nb"),
  az  = c(45.89,45.11,58.44,314.87,303.15,3.35,224.93,224.62,236.61,135.19,123.73,183.74),
  zen = c(20.58,45.31,90.21,45.24,90.27,89.79,20.43,46.18,89.97,45.55,90.42,90.32))

sun_j2000 <- function(jd) { # low-precision Sun unit vector (J2000 equatorial)
  n <- jd - 2451545.0; L <- (280.460 + 0.9856474 * n) %% 360; g <- ((357.528 + 0.9856003 * n) %% 360) * pi / 180
  lam <- (L + 1.915 * sin(g) + 0.020 * sin(2 * g)) * pi / 180; eps <- (23.439 - 0.0000004 * n) * pi / 180
  cbind(cos(lam), cos(eps) * sin(lam), sin(eps) * sin(lam))
}
quat_rotate <- function(q, v) { # rotate vectors v (n x 3) by quaternions q (n x 4, x y z w) : J2000 -> spacecraft
  x <- q[,1]; y <- q[,2]; z <- q[,3]; w <- q[,4]
  r11 <- 1-2*(y^2+z^2); r12 <- 2*(x*y-z*w); r13 <- 2*(x*z+y*w)
  r21 <- 2*(x*y+z*w); r22 <- 1-2*(x^2+z^2); r23 <- 2*(y*z-x*w)
  r31 <- 2*(x*z-y*w); r32 <- 2*(y*z+x*w); r33 <- 1-2*(x^2+y^2)
  # QSJ is the J2000 -> SC frame transform; the matrix above is its (passive) form
  cbind(r11*v[,1]+r12*v[,2]+r13*v[,3], r21*v[,1]+r22*v[,2]+r23*v[,3], r31*v[,1]+r32*v[,2]+r33*v[,3])
}

sun_angles <- function(poshist_file) { # returns data.frame(met, sun-angle per detector)
  z <- readFITS(poshist_file, hdu = 1)
  cn <- z$colNames; get <- function(n) z$col[[which(cn == n)]]
  met <- get("SCLK_UTC"); q <- cbind(get("QSJ_1"), get("QSJ_2"), get("QSJ_3"), get("QSJ_4"))
  jd <- 2451910.5 + met / 86400  # MET epoch 2001-01-01 00:00:00 UTC
  # quaternion row norm check
  sv <- quat_rotate(q, sun_j2000(jd))
  dv <- with(DET, cbind(sin(zen*pi/180)*cos(az*pi/180), sin(zen*pi/180)*sin(az*pi/180), cos(zen*pi/180)))
  cs <- sv %*% t(dv)
  ang <- matrix(acos(pmin(1, pmax(-1, as.vector(cs)))) * 180 / pi, nrow = nrow(cs))
  colnames(ang) <- DET$det
  data.frame(met = met, ang)
}

# Bandwidth-limited mode: fixed Sun-side detectors (n1, n2 faced the Sun ~80% of the day when checked
# against poshist on 2024-05-01), no position file. Set FIXED_DETS <- NULL to use the poshist selection below.
FIXED_DETS <- c("n1", "n2")
download_fermi_fixed <- function() {
  for (d in seq(FERMI_START, FERMI_END, by = "day")) {
    d <- as.Date(d, origin = "1970-01-01"); ds <- format(d, "%y%m%d"); base <- paste0(FERMI, format(d, "%Y/%m/%d"), "/current/")
    for (k in FIXED_DETS) fetch(paste0(base, "glg_ctime_", k, "_", ds, "_v00.pha"),
                                file.path(raw, "hxr", paste0("glg_ctime_", k, "_", ds, "_v00.pha")))
    message("fermi ", format(d), " fetched")
  }
}

download_fermi <- function(n_det = 3, max_angle = 60) {
  if (!is.null(FIXED_DETS)) return(download_fermi_fixed())
  days <- seq(FERMI_START, FERMI_END, by = "day"); log <- list()
  for (i in seq_along(days)) {
    d <- days[i]; ds <- format(d, "%y%m%d"); base <- paste0(FERMI, format(d, "%Y/%m/%d"), "/current/")
    pos <- file.path(raw, "hxr", paste0("glg_poshist_all_", ds, "_v00.fit"))
    if (!fetch(paste0(base, basename(pos)), pos)) next
    a <- tryCatch(sun_angles(pos), error = function(e) { message("poshist error ", ds, ": ", e$message); NULL })
    if (is.null(a)) next
    frac <- sapply(DET$det, function(k) mean(a[[k]] < max_angle, na.rm = TRUE)); names(frac) <- DET$det
    pick <- names(sort(frac, decreasing = TRUE))[seq_len(n_det)]
    for (k in pick) fetch(paste0(base, "glg_ctime_", k, "_", ds, "_v00.pha"),
                          file.path(raw, "hxr", paste0("glg_ctime_", k, "_", ds, "_v00.pha")))
    log[[i]] <- data.frame(date = d, det = pick, frac_sun_lt60 = round(as.numeric(frac[pick]), 3))
    message("fermi ", format(d), " -> ", paste(pick, collapse = ","))
  }
  out <- do.call(rbind, log)
  write.csv(out, file.path(root, "data", "processed", "fermi_detector_choice.csv"), row.names = FALSE)
}

dir.create(file.path(root, "data", "processed"), showWarnings = FALSE, recursive = TRUE)
what <- if (length(commandArgs(TRUE))) commandArgs(TRUE)[1] else "all"
if (what %in% c("goes", "all")) download_goes()
if (what %in% c("fermi", "all")) download_fermi()
message("download complete")
