# Install and check required packages (run once)
pkgs <- c("ncdf4", "FITSio", "data.table", "ranger", "xgboost", "glmnet",
          "pROC", "ggplot2", "rmarkdown", "zoo", "knitr")
miss <- setdiff(pkgs, rownames(installed.packages()))
if (length(miss)) {
  options(repos = c(CRAN = "https://cloud.r-project.org"))
  install.packages(miss, type = "binary")
}
ok <- sapply(pkgs, requireNamespace, quietly = TRUE)
print(ok)
