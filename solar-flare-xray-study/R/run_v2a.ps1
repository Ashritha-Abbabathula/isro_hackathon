# run_v2a.ps1 -- full rerun on the existing 2-month window with magnetic + sub-minute features (processed_v2a)
$d = "C:\Users\jayag\OneDrive\Desktop\solar-flare-xray"; $rs = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"; $p = "C:/Users/jayag/OneDrive/Desktop/solar-flare-xray"
function Step($name, $script, $envs) {
  foreach ($k in $envs.Keys) { Set-Item "Env:$k" $envs[$k] }
  "$(Get-Date -Format s) START $name" | Add-Content "$d\data\chain_v2a.log"
  Start-Process -FilePath $rs -ArgumentList "`"$d\R\$script`"" -RedirectStandardError "$d\data\v2a_$name.log" -RedirectStandardOutput "$d\data\v2a_$name.out" -WindowStyle Hidden -Wait
  "$(Get-Date -Format s) END   $name" | Add-Content "$d\data\chain_v2a.log"
  foreach ($k in $envs.Keys) { Remove-Item "Env:$k" -ErrorAction SilentlyContinue } }
$base = @{ SOLAR_PROC = "$p/data/processed_v2a" }
Step "loo"      "21_models_v2.R"        ($base + @{ SOLAR_OUT = "$p/outputs_v2a"; MODE = "loo" })
Step "forward"  "21_models_v2.R"        ($base + @{ SOLAR_OUT = "$p/outputs_v2a_forward"; MODE = "forward"; MODELS = "glmnet,xgboost" })
Step "attr_loo" "12_attribution_adl1.R" @{ SOLAR_OUT = "$p/outputs_v2a" }
Step "attr_fwd" "12_attribution_adl1.R" @{ SOLAR_OUT = "$p/outputs_v2a_forward" }
Step "robust"   "16_robustness_adl1.R"  ($base + @{ SOLAR_OUT = "$p/outputs_v2a"; SOLAR_OUT_FWD = "$p/outputs_v2a_forward" })
Step "causal"   "15_causality_test.R"   ($base + @{ SOLAR_OUT = "$p/outputs_v2a"; FEATURE_SCRIPT = "20_features_v2.R"; FEATURE_FILE = "features_v2.csv" })
"$(Get-Date -Format s) ALL DONE" | Add-Content "$d\data\chain_v2a.log"
