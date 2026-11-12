# Runs the whole pipeline in order. Usage: powershell -File R\run_pipeline.ps1
$root = Split-Path -Parent $PSScriptRoot
$rs = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
$steps = @("01_download.R fermi", "01b_verify.R", "01_download.R all", "02_clean_align.R", "02b_fermi_hxr.R", "03_features.R", "04_models.R", "05_attribution.R")
foreach ($s in $steps) {
  $parts = $s.Split(" "); $script = Join-Path $root ("R\" + $parts[0])
  if (-not (Test-Path $script)) { "skip (missing) $s" | Add-Content "$root\data\pipeline.log"; continue }
  "$(Get-Date -Format s) START $s" | Add-Content "$root\data\pipeline.log"
  $argList = @("`"$script`"") + $parts[1..($parts.Length)] | Where-Object { $_ }
  Start-Process -FilePath $rs -ArgumentList $argList -RedirectStandardError "$root\data\step_$($parts[0]).log" -RedirectStandardOutput "$root\data\step_$($parts[0]).out" -WindowStyle Hidden -Wait
  "$(Get-Date -Format s) END   $s" | Add-Content "$root\data\pipeline.log"
}
