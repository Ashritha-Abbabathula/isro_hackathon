# run_v2.ps1 -- FULL 6-month (1 Jul - 31 Dec 2024) pipeline: preprocessing -> features -> models (LOO + forward + SEP-included) -> attribution -> robustness -> causality -> nowcasting
# Start only after: (1) HEL1OS extension download finished, (2) GOES retry finished (EUV + flare summary for Jul-Sep), (3) sharp_fulldisk_12min.csv exists.
$p = "C:/Users/jayag/OneDrive/Desktop/solar-flare-xray"; $d = "C:\Users\jayag\OneDrive\Desktop\solar-flare-xray"; $rs = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
$proc = "$p/data/processed_v2"; $log = "$d\data\chain_v2.log"
function Stamp($m) { "$(Get-Date -Format s) $m" | Add-Content $log }
function Step($name, $script, $envs) {
  foreach ($k in $envs.Keys) { Set-Item "Env:$k" $envs[$k] }
  Stamp "START $name"
  Start-Process -FilePath $rs -ArgumentList "`"$d\R\$script`"" -RedirectStandardError "$d\data\v2_$name.log" -RedirectStandardOutput "$d\data\v2_$name.out" -WindowStyle Hidden -Wait
  Stamp "END   $name"
  foreach ($k in $envs.Keys) { Remove-Item "Env:$k" -ErrorAction SilentlyContinue } }
# 1. preprocessing (numpy only) for all 6 months
Stamp "START preprocess"
Start-Process -FilePath "python.exe" -ArgumentList "-u", "`"$d\scripts_py\preprocess_aditya_np.py`"", "`"$d\data\processed_v2`"", "20240701", "20241231" -RedirectStandardOutput "$d\data\v2_preprocess.log" -RedirectStandardError "$d\data\v2_preprocess.err" -WindowStyle Hidden -Wait
Stamp "END   preprocess"
# 2. GOES cleaning (all 6 months) -> goes_1min.csv, flare_events.csv
Step "goes" "02_clean_align.R" @{ SOLAR_RAW = "C:/solar-data/goes_aditya_period"; SOLAR_PROC = $proc }
# 3. static inputs
Copy-Item "$d\data\processed_adl1\sharp_fulldisk_12min.csv" "$d\data\processed_v2\" -Force
Copy-Item "$d\data\processed_v2\sep_periods_6mo.csv" "$d\data\processed_v2\sep_periods.csv" -Force
# 4. features, models, attribution, robustness, tests, nowcasting
$base = @{ SOLAR_PROC = $proc }
Step "features" "20_features_v2.R" $base
Step "loo"      "21_models_v2.R"  ($base + @{ SOLAR_OUT = "$p/outputs_v2"; MODE = "loo" })
Step "forward"  "21_models_v2.R"  ($base + @{ SOLAR_OUT = "$p/outputs_v2_forward"; MODE = "forward"; MODELS = "glmnet,xgboost" })
Step "attr_loo" "12_attribution_adl1.R" @{ SOLAR_OUT = "$p/outputs_v2" }
Step "attr_fwd" "12_attribution_adl1.R" @{ SOLAR_OUT = "$p/outputs_v2_forward" }
Step "nowcast"  "13_eda_nowcast_adl1.R" ($base + @{ SOLAR_OUT = "$p/outputs_v2" })
Step "robust"   "16_robustness_adl1.R"  ($base + @{ SOLAR_OUT = "$p/outputs_v2"; SOLAR_OUT_FWD = "$p/outputs_v2_forward" })
Step "causal"   "15_causality_test.R"   ($base + @{ SOLAR_OUT = "$p/outputs_v2"; FEATURE_SCRIPT = "20_features_v2.R"; FEATURE_FILE = "features_v2.csv" })
Step "sepincl"  "21_models_v2.R"  ($base + @{ SOLAR_OUT = "$p/outputs_v2_sepincl"; MODE = "loo"; SEP_MODE = "include"; MODELS = "glmnet,xgboost" })
Step "attr_sep" "12_attribution_adl1.R" @{ SOLAR_OUT = "$p/outputs_v2_sepincl" }
Stamp "ALL DONE"
