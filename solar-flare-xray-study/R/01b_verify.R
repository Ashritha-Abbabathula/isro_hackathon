# 01b_verify.R -- integrity check of downloaded files; deletes corrupt ones so 01_download.R can refetch them
suppressPackageStartupMessages(library(ncdf4)); suppressPackageStartupMessages(library(FITSio))
root <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."), winslash = "/")
raw <- file.path(root, "data", "raw"); bad <- character(); n_ok <- 0
for (f in list.files(raw, pattern = "\\.nc$", recursive = TRUE, full.names = TRUE)) {
  ok <- tryCatch({ n <- nc_open(f); t <- ncvar_get(n, "time"); nc_close(n); length(t) > 0 }, error = function(e) FALSE)
  if (ok) n_ok <- n_ok + 1 else bad <- c(bad, f)
}
for (f in list.files(file.path(raw, "hxr"), pattern = "ctime.*\\.pha$", full.names = TRUE)) {
  ok <- tryCatch({ z <- readFITS(f, hdu = 2); length(z$col[[which(z$colNames == "TIME")]]) > 1000 }, error = function(e) FALSE, warning = function(w) TRUE)
  if (ok) n_ok <- n_ok + 1 else bad <- c(bad, f)
}
message("valid: ", n_ok, "  corrupt: ", length(bad))
if (length(bad)) { print(table(basename(dirname(bad)))); unlink(bad) }
