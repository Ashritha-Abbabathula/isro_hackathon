$d = "C:\Users\jayag\OneDrive\Desktop\solar-flare-xray"; $rs = "C:\Program Files\R\R-4.6.1\bin\Rscript.exe"
while (-not (Test-Path "$d\outputs_adl1\tables\results.csv")) { Start-Sleep -Seconds 20 }
Start-Sleep -Seconds 5
function Run($script, $log, $envs) { foreach ($k in $envs.Keys) { Set-Item "Env:$k" $envs[$k] }
  Start-Process -FilePath $rs -ArgumentList "`"$d\R\$script`"" -RedirectStandardError "$d\data\$log.log" -RedirectStandardOutput "$d\data\$log.out" -WindowStyle Hidden -Wait
  foreach ($k in $envs.Keys) { Remove-Item "Env:$k" -ErrorAction SilentlyContinue } }
"$(Get-Date -Format s) main models done" | Add-Content "$d\data\chain.log"
Run "12_attribution_adl1.R" "attr_adl1" @{}
"$(Get-Date -Format s) attribution done" | Add-Content "$d\data\chain.log"
Run "11_models_adl1.R" "models_sepincl" @{ SEP_MODE = "include"; SOLAR_OUT = "$d\outputs_adl1_sepincl"; MODELS = "glmnet,xgboost" }
"$(Get-Date -Format s) sep-include models done" | Add-Content "$d\data\chain.log"
Run "12_attribution_adl1.R" "attr_sepincl" @{ SOLAR_OUT = "$d\outputs_adl1_sepincl" }
"$(Get-Date -Format s) ALL DONE" | Add-Content "$d\data\chain.log"
