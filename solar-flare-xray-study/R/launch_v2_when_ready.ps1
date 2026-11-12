# launch_v2_when_ready.ps1 -- waits for the HEL1OS extension and the GOES retry to finish, then starts the full 6-month chain (run_v2.ps1).
$d = "C:\Users\jayag\OneDrive\Desktop\solar-flare-xray"; $log = "$d\data\chain_v2.log"
function Stamp($m) { "$(Get-Date -Format s) $m" | Add-Content $log }
Stamp "watcher started"
while ($true) {
  $hel = (Get-Content "C:\solar-data\aditya\hel1os\bulk_hel1os_ext.log" -ErrorAction SilentlyContinue | Select-String "Your downloads \(" -Quiet)
  $goes = (Get-Content "C:\solar-data\goes_aditya_period\dl_ext3.log" -ErrorAction SilentlyContinue | Select-String "download complete" -Quiet)
  $v2a = (Get-Content "$d\data\chain_v2a.log" -ErrorAction SilentlyContinue | Select-String "ALL DONE" -Quiet)   # do not compete for the CPU with the interim rerun
  if ($hel -and $goes -and $v2a) { break }
  Start-Sleep -Seconds 60
}
Stamp "downloads finished and interim rerun done - launching run_v2.ps1"
& powershell -NoProfile -ExecutionPolicy Bypass -File "$d\R\run_v2.ps1"

